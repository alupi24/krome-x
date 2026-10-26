#!/usr/bin/env python3
"""Checks of a network and of the code generated from it (python3 -m kromex
--check).

1. The generated code against an independent evaluation of the network
   (a Python interpreter of the parsed network, not the generated C++): the
   rate coefficients at sample temperatures, the right-hand side from the
   stoichiometry, the analytic Jacobian against finite differences of the
   generated right-hand side, the freeze directives. A mismatch is a
   failure (a code-generation problem: please report it).
2. The rates themselves, at the sample temperatures inside each reaction's
   T range (the solver evaluates them from T_CMB to 1e9 K): a rate that
   is not finite is a failure; a negative rate, a two-body rate above 1e-5
   cm^3 s^-1 or a three-body rate above 1e-24 cm^6 s^-1 is a warning
   (often a fit used beyond its validity: declare its T range; or a typo
   in an exponent or the units).
3. The network as a whole: species formed but never destroyed, destroyed
   but never formed, or in no reaction (warnings); the same reaction
   written more than once with the same rate (a warning: counted twice) or
   with different rates (a note: summed as channels); detailed balance
   (balance.py): collisional reactions that run one way although their
   reverse, missing from the network, reaches 1 % of the forward flux
   inside their temperature range (a warning, the first ones named), and
   endothermic ones whose exothermic reverse is missing (a note).

Parameters are set to a probe value (1.3e-17; Tdust 10 K): a rate that
depends on a parameter is checked for consistency, but its value is only
reported as a warning.

The development test suite runs the code checks (1) on every shipped
network through tests/selfcheck.py."""

import paths  # noqa: E402 (where the runtime is)
import math
import os
import random
import shutil
import subprocess
from pathlib import Path

import knet
import codegen
import constants

KB_EV = constants.boltzmann_eV


# ---------------------------------------------------------------------------
# independent Python evaluation of the IR
# ---------------------------------------------------------------------------

def _colion_v96(T, dE, P, A, X, K):   # legacy krome_subs.f90
    U = dE / (T * constants.boltzmann_eV)
    return A * (1.0 + P * math.sqrt(U)) * U ** K * math.exp(-U) / (X + U)


def _recNewV96(T, r1, r2, r3, r4):
    tt = math.sqrt(T / r3)
    return r1 / (tt * (tt + 1.0) ** (1.0 - r2)
                 * (1.0 + math.sqrt(T / r4)) ** (1.0 + r2))


def _recFeV96(T, r1, r2, r3):
    tt = math.sqrt(T * 1e-4)
    return r1 / tt ** (r2 + r3 + math.log10(tt))


# independent Python versions of the legacy rate-database helpers that
# resolved `auto` rates call (codegen AUTO_TFUNCS)
TFUNCS = {
    "colion_v96": _colion_v96,
    "recV96": lambda T, a, b: a * (1e4 / T) ** b,
    "recNewV96": _recNewV96,
    "recFeV96": _recFeV96,
}


def reverse_tfuncs(gen):
    """independent 1/K_c(T) of the generated reverse rates, from the same
    thermo data (NASA 7-coefficient, else NIST Shomate):
    exp(dG/RT) (kB T / 1 bar)^dn"""
    nasa, nist = gen.thermo_all()

    def gRT(sp, T):
        c = nasa.get(sp)
        if c is not None:
            T = min(max(T, c[0]), c[2])
            a = c[10:17] if T <= c[1] else c[3:10]
            H = (a[0] + a[1] * T / 2 + a[2] * T * T / 3 + a[3] * T ** 3 / 4
                 + a[4] * T ** 4 / 5 + a[5] / T)
            S = (a[0] * math.log(T) + a[1] * T + a[2] * T * T / 2
                 + a[3] * T ** 3 / 3 + a[4] * T ** 4 / 4 + a[6])
            return H - S
        r = nist[sp.upper()]
        T = min(max(T, r[0][0]), r[-1][1])
        a = next((q[2] for q in r if T <= q[1]), r[-1][2])
        t = T * 1e-3
        H = (a[0] * t + a[1] * t * t / 2 + a[2] * t ** 3 / 3
             + a[3] * t ** 4 / 4 - a[4] / t + a[5])
        S = (a[0] * math.log(t) + a[1] * t + a[2] * t * t / 2
             + a[3] * t ** 3 / 3 - a[4] / (2 * t * t) + a[6])
        return H / (constants.Rgas_kJ * T) - S / constants.Rgas_J

    out = {}
    for q, rx in enumerate(getattr(gen, "_reverse", [])):
        R = [r["sp"] for r in rx["reactants"]]
        P = [p0["sp"] for p0 in rx["products"]]

        def f(T, R=R, P=P):
            T = max(T, 1e-10)
            dg = sum(gRT(x, T) for x in P) - sum(gRT(x, T) for x in R)
            return math.exp(min(dg, 700.0)) * (constants.boltzmann_erg * T * 1e-6) ** (
                len(P) - len(R))
        out["revKc_%d" % q] = f
    return out


def _growth(args, env):
    """kinetic nucleation (Gail & Sedlmayr 2013): v_th pi (ra + rb)^2 S"""
    net = env["net"]
    m = args[0]["sp"]
    if len(args) == 3:   # cluster_growth_rate(M, N, stick)
        N1, N2, st = 1.0, eval_ast(args[1], env), eval_ast(args[2], env)
    else:                # general_cluster_growth_rate(M, N1, N2, stick)
        N1, N2, st = (eval_ast(args[1], env), eval_ast(args[2], env),
                      eval_ast(args[3], env))
    mass = net.species[m]["mass_amu"] * constants.amu
    r1 = codegen.NUCLEATION_RADII[m.upper()]
    v = math.sqrt(8 * constants.boltzmann_erg * env["Tgas"] * (1 / (N1 * mass)
                                                     + 1 / (N2 * mass))
                  / constants.pi)
    ra, rb = r1 * N1 ** (1 / 3), r1 * N2 ** (1 / 3)
    return v * constants.pi * (ra + rb) ** 2 * st


def _troe(args, env):
    k0, kinf, Fc, m = (eval_ast(x, env) for x in args)
    if not (k0 * m > 0 and kinf > 0):
        return 0.0
    rm = k0 * m / kinf
    return k0 * m / (1 + rm) * Fc ** (1 / (1 + math.log10(rm) ** 2))


# independent Python versions of generated extension functions
PY_EXTERNS = {"cluster_growth_rate": _growth,
              "general_cluster_growth_rate": _growth,
              "troe_falloff": _troe}


def gfe_tfuncs(gen):
    """independent K_c of the revKc_with_GFE sites from the .gfe tables
    (linear in T, clamped to the grid)"""
    out = {}
    data = getattr(gen, "_gfe_data", {})

    def G(sp, T):
        Ts, Gs = data[sp]
        if T <= Ts[0]:
            return Gs[0]
        if T >= Ts[-1]:
            return Gs[-1]
        i = min(int((T - Ts[0]) / (Ts[1] - Ts[0])), len(Ts) - 2)
        f = (T - Ts[i]) / (Ts[i + 1] - Ts[i])
        return Gs[i] + f * (Gs[i + 1] - Gs[i])
    for q, (R, P) in enumerate(getattr(gen, "_gfe", [])):
        def kc(T, R=R, P=P):
            dg = sum(G(x, T) for x in P) - sum(G(x, T) for x in R)
            return math.exp(min(-dg / (constants.Rgas_kJ * T), 700.0)) * (
                constants.boltzmann_erg * T * 1e-6) ** (len(R) - len(P))
        out["gfeKc_%d" % q] = kc
    return out


def eval_ast(a, env):
    if "n" in a:
        return a["n"]
    if "s" in a:
        k, s = a.get("kind"), a["s"]
        if k == "state":
            if a.get("unclamped"):  # 1/K_c of reverse rates: actual T
                return env["Tgas_actual"]
            return env[s]
        if k == "param":
            return env["params"][s]
        if k == "var":
            return eval_ast(env["net"].vars[s], env)
        if k == "const":
            return knet.CONSTANTS[s]
    if "u" in a:
        return -eval_ast(a["x"], env)
    if "c" in a:
        l, r = eval_ast(a["l"], env), eval_ast(a["r"], env)
        return {"<": l < r, "<=": l <= r, ">": l > r, ">=": l >= r,
                "==": l == r, "!=": l != r}[a["c"]]
    if "b" in a:
        l, r = eval_ast(a["l"], env), eval_ast(a["r"], env)
        op = a["b"]
        if op == "+":
            return l + r
        if op == "-":
            return l - r
        if op == "*":
            return l * r
        if op == "/":
            return l / r
        return l ** r
    if "f" in a:
        f, args = a["f"], a["a"]
        if a.get("kind") == "density":
            return env["n"][args[0]["sp"]]
        if f == "where":
            return (eval_ast(args[1], env) if eval_ast(args[0], env)
                    else eval_ast(args[2], env))
        if a.get("kind") == "extern" and f in PY_EXTERNS:
            return PY_EXTERNS[f](args, env)
        vals = [eval_ast(x, env) for x in args]
        if a.get("kind") == "tfunc":
            return TFUNCS[f](*vals)
        return {"exp": math.exp, "log": math.log, "log10": math.log10,
                "sqrt": math.sqrt, "abs": abs, "min": min, "max": max,
                "tanh": math.tanh}[f](*vals)
    raise ValueError("cannot evaluate %r" % a)


def py_rates(net, Tgas, ndens, params, supported, outside="asym"):
    """the rate coefficients; outside a reaction's T range (rates.outside):
    asym zero below and the value at the upper bound above, clamp the
    value at the nearer bound, zero zero"""
    env = {"net": net, "Tgas": Tgas, "Tgas_actual": Tgas,
           "Tgas_eV": Tgas * KB_EV,
           "Tdust": params.get("Tdust", 10.0),
           "ntot": sum(ndens.values()),
           "Hnuclei": sum(v * (net.species[s]["composition"] or {})
                          .get("H", 0) for s, v in ndens.items()),
           "n": ndens, "params": params}
    k = [0.0] * len(net.reactions)
    def rate(rx, e):
        # rate_scale applied here, on the unscaled AST, independently of
        # the folding codegen does
        r = rx["rate"]
        try:  # C++ gives NaN / Inf where Python raises
            v = eval_ast(r.get("ast_unscaled", r["ast"]), e)
            if net.rate_scale is not None:
                v *= eval_ast(net.rate_scale, e)
        except (ValueError, ZeroDivisionError):
            return math.nan
        except OverflowError:
            return math.inf
        return math.nan if isinstance(v, complex) else v
    groups = {}
    for rx in net.reactions:
        if rx["id"] - 1 not in supported:
            continue
        if rx.get("group") is None:
            k[rx["id"] - 1] = rate(rx, env)
        else:
            groups.setdefault(rx["group"], []).append(rx)
    for pieces in groups.values():
        pieces.sort(key=lambda r: r["trange"])
        glo, ghi = pieces[0]["trange"][0], pieces[-1]["trange"][1]
        if Tgas < glo and outside != "clamp":
            continue  # zero below the lower bound
        if Tgas > ghi and outside == "zero":
            continue
        T = min(max(Tgas, glo), ghi)
        e2 = dict(env, Tgas=T, Tgas_eV=T * KB_EV)
        for j, rx in enumerate(pieces):
            lo, hi = rx["trange"]
            last = j == len(pieces) - 1
            if lo <= T and (T < hi or (last and T <= hi)):
                k[rx["id"] - 1] = rate(rx, e2)
                break
    return k


def py_rhs(net, sp_sorted, ndens, k):
    dn = {s: 0.0 for s in sp_sorted}
    for rx in net.reactions:
        flux = k[rx["id"] - 1]
        for r in rx["reactants"]:
            if r["kinetic"]:
                flux *= ndens[r["sp"]]
        for r in rx["reactants"]:
            dn[r["sp"]] -= flux
        for pr in rx["products"]:
            dn[pr["sp"]] += flux
    return [dn[s] for s in sp_sorted]


# ---------------------------------------------------------------------------
# C++ probe
# ---------------------------------------------------------------------------

PROBE = r"""
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "@HEADER@"
#include "core/dust_rates.hpp"
#include "core/h2_shielding.hpp"
#include "core/special_rates.hpp"
using namespace krome;
int main(int argc, char** argv) {
  double Tgas = atof(argv[1]);
  double n[NSPEC], k[NREA], dn[NSPEC];
  static double J[NSPEC*NSPEC];
  for (int i = 0; i < NSPEC; ++i) n[i] = atof(argv[2 + i]);
  Params p;
  for (int i = 0; i < NPAR; ++i)  // Tdust 10: matches the Python side
    *param_ptr(p, i) = std::strcmp(param_names[i], "Tdust") ? 1.3e-17 : 10.0;
  rates(Tgas, n, p, k);
  rhs(n, k, dn);
  jacobian(Tgas, n, p, k, J);
  double sc[NSPEC];
  rhs_scale(n, k, sc);
  for (int i = 0; i < NREA; ++i) printf("k %d %.17e\n", i, k[i]);
  for (int i = 0; i < NSPEC; ++i) printf("d %d %.17e\n", i, dn[i]);
  for (int i = 0; i < NSPEC; ++i) printf("s %d %.17e\n", i, sc[i]);
  for (int i = 0; i < NSPEC*NSPEC; ++i) printf("J %d %.17e\n", i, J[i]);
  int fr[NFREEZE_MAX + 1];
  const int nf = frozen_species(p, fr);
  for (int i = 0; i < nf; ++i) printf("f %d %d\n", i, fr[i]);
  return 0;
}
"""


RUNTIME = paths.ROOT / "runtime"
T_SAMPLES = (12.0, 97.0, 314.0, 1.1e3, 5.5e3, 2.3e4, 1.7e5, 3.0e6, 8e7)
PROBE_PARAM = 1.3e-17
MAX_RATE = {2: 1e-5, 3: 1e-24}  # warning thresholds, cm^3/s and cm^6/s


def relerr(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def compiler():
    """the C++ compiler of the probe: $CXX, else c++, clang++ or g++"""
    for c in (os.environ.get("CXX"), "c++", "clang++", "g++"):
        if c and shutil.which(c):
            return c
    raise RuntimeError("no C++ compiler found (set CXX)")


def _order(rx):
    return sum(1 for r in rx["reactants"] if r["kinetic"])


def _uses_param(gen, rx):
    return any(nd.get("kind") == "param"
               for nd in gen.deep_walk(rx["rate"]["ast"]))


def _in_range(rx, T):
    tr = rx["trange"]
    return not tr or tr[0] <= T <= tr[1]


def budget_warnings(net):
    """species formed but never destroyed, destroyed but never formed, in
    no reaction"""
    formed, destroyed, seen = set(), set(), set()
    for rx in net.reactions:
        d = {}
        for r in rx["reactants"]:
            d[r["sp"]] = d.get(r["sp"], 0) - 1
        for pr in rx["products"]:
            d[pr["sp"]] = d.get(pr["sp"], 0) + 1
        for sp, v in d.items():
            seen.add(sp)
            if v > 0:
                formed.add(sp)
            elif v < 0:
                destroyed.add(sp)
    live = [s for s in sorted(net.species)
            if not net.species[s].get("opaque")]
    out = []
    for label, sel in (
            ("formed but never destroyed (they can only accumulate)",
             [s for s in live if s in formed and s not in destroyed]),
            ("destroyed but never formed (they can only be depleted)",
             [s for s in live if s in destroyed and s not in formed]),
            ("in no reaction, or only as catalysts (they never change)",
             [s for s in live if s not in formed and s not in destroyed])):
        if sel:
            out.append("species %s: %s" % (label, ", ".join(sel)))
    return out


def balance_warnings(net, gen, shown=8):
    """detailed balance (balance.analyse): a warning naming the one-way
    reactions whose reverse sets in at the lowest temperatures (the first
    `shown`), a note on the endothermic ones; ([warnings], [notes])"""
    import balance
    a = balance.analyse(net, balance.load_tables(gen.data_dir))
    warn, note = [], []
    fl = a["flagged"]
    if fl:
        def first(e):
            return min(v for v in e["cross"].values() if v is not None)
        fl = sorted(fl, key=first)
        warn.append(
            "detailed balance: %d collisional reaction%s run%s one way "
            "although %s reverse, not in the network, reaches %g %% of the "
            "forward flux inside %s temperature range (with every species "
            "at 1 to 1e12 cm^-3), e.g. %s%s: where the network is used "
            "there, add the reverse (A <-> B); kromex-doc lists them all "
            "(Reverse reactions)" % (
                len(fl), "" if len(fl) == 1 else "s",
                "s" if len(fl) == 1 else "", "its" if len(fl) == 1 else
                "their", 100 * balance.THRESHOLD,
                "its" if len(fl) == 1 else "their",
                "; ".join("%s (line %s) above %.0f K" % (
                    e["eq"], e["line"], first(e)) for e in fl[:shown]),
                "; ..." if len(fl) > shown else ""))
    en = a["endothermic"]
    if en:
        n3 = sum(1 for e in en if len(e["P"]) >= 3)
        note.append(
            "detailed balance: %d endothermic reaction%s without %s "
            "reverse, which is the dominant direction (an exothermic "
            "association or recombination; %d three-body, where the "
            "forward makes more particles: important at high density): "
            "%s%s" % (
                len(en), "" if len(en) == 1 else "s",
                "its" if len(en) == 1 else "their", n3,
                "; ".join(e["eq"] for e in en[:shown]),
                "; ..." if len(en) > shown else ""))
    return warn, note


def repeated(net):
    """the same reaction written more than once, applying at a common T
    (not the pieces of one piecewise rate): its rates are summed.
    (warnings: the same rate expression, notes: different ones)"""
    by = {}
    for rx in net.reactions:
        key = (tuple(sorted(r["sp"] for r in rx["reactants"])),
               tuple(sorted(p["sp"] for p in rx["products"])))
        by.setdefault(key, []).append(rx)
    warn, note = [], []
    for rxs in by.values():
        if len(rxs) < 2:
            continue
        clash = same = False
        for a in range(len(rxs)):
            for b in range(a + 1, len(rxs)):
                ra, rb = rxs[a]["trange"], rxs[b]["trange"]
                ga, gb = rxs[a].get("group"), rxs[b].get("group")
                if ga is not None and ga == gb:
                    continue  # pieces of one piecewise rate
                if not ra or not rb or (ra[0] < rb[1] and rb[0] < ra[1]):
                    clash = True
                    same |= (rxs[a]["rate"].get("ast") is not None
                             and rxs[a]["rate"].get("ast")
                             == rxs[b]["rate"].get("ast"))
        lines = ", ".join(str(r.get("line", "?")) for r in rxs)
        if same:
            warn.append("%s written %d times with the same rate (lines %s): "
                        "counted twice" % (rxs[0]["eq"], len(rxs), lines))
        elif clash:
            note.append("%s (lines %s)" % (rxs[0]["eq"], lines))
    if note:
        note = ["reactions written more than once with different rates, "
                "summed as channels: " + "; ".join(note)]
    return warn, note


def independent(net, gen):
    """the reactions whose rate the Python evaluation reproduces: analytic
    expressions without tables or external functions (photo, surface and
    tabulated rates are computed at run time)"""
    # with --dust-tables, Tdust is the table's (a function of the state),
    # not the parameter the independent evaluation uses
    tdust_table = bool(getattr(gen, "dust_tables", ""))

    def ok(rx):
        if rx["rate"]["kind"] != "expr":
            return False
        return not any((nd.get("kind") in ("extern", "table")
                        and nd.get("f") not in PY_EXTERNS)
                       or (tdust_table and nd.get("kind") == "state"
                           and nd.get("s") == "Tdust")
                       for nd in gen.deep_walk(rx["rate"]["ast"]))
    return [rx["id"] - 1 for rx in net.reactions if ok(rx)]


CHECK_T = (97.0, 5.5e3, 1.7e5)  # make check: the reference temperatures

CHECK_CPP = r"""// make check: the rate coefficients of this build against reference
// values computed when the project was generated, by an independent
// evaluation of the network (kromex/netcheck.py). A difference
// means the compiler or its flags changed the arithmetic (e.g.
// -ffast-math) -- regenerated with the project, do not edit.
#include <cmath>
#include <cstdio>
#include <cstring>
#include "krome.hpp"
using namespace krome;

static const double T_REF[@NT@] = {@T@};
static const double N_REF[NSPEC] = {@N@};
static const int I_REF[@NI@] = {@I@};
static const double K_REF[@NT@][@NI@] = {@K@};

int main() {
  Params p;
  for (int i = 0; i < NPAR; ++i)
    *param_ptr(p, i) = std::strcmp(param_names[i], "Tdust") ? @P@ : 10.0;
  double n[NSPEC], k[NREA];
  for (int i = 0; i < NSPEC; ++i) n[i] = N_REF[i];
  double worst = 0.0, wa = 0.0, wb = 0.0;
  int wt = 0, wi = 0;
  for (int t = 0; t < @NT@; ++t) {
    rates(T_REF[t], n, p, k);
    for (int j = 0; j < @NI@; ++j) {
      const double a = k[I_REF[j]], b = K_REF[t][j];
      const double d = std::fabs(a - b) /
                       std::fmax(std::fmax(std::fabs(a), std::fabs(b)), 1e-300);
      if (!(d <= worst)) { worst = d; wt = t; wi = j; wa = a; wb = b; }
    }
  }
  const bool ok = worst < 1e-8;
  printf("make check: %d rate coefficients at %d temperatures, max relative "
         "difference %.1e: %s\n", @NI@, @NT@, worst, ok ? "OK" : "FAILED");
  if (!ok)
    printf("  worst: reaction %d at T = %g K: %.9e, reference %.9e\n",
           I_REF[wi] + 1, T_REF[wt], wa, wb);
  return ok ? 0 : 1;
}
"""


def check_source(net, gen):
    """krome_check.cpp (make check): reference rates of the independently
    evaluated reactions at CHECK_T (a program that only says so when there
    are none)"""
    sp = sorted(net.species)
    rng = random.Random(7)
    ndens = {s: 10 ** rng.uniform(-6, 3) for s in sp}
    state = dict(ndens)
    for s_, g in getattr(gen, "pair_gas", {}).items():
        state[s_] = ndens[s_] + ndens[g]
    params = {q: PROBE_PARAM for q in gen.pnames}
    params["Tdust"] = 10.0
    TFUNCS.update(reverse_tfuncs(gen))
    TFUNCS.update(gfe_tfuncs(gen))
    sup = independent(net, gen)
    rows = []
    for T in CHECK_T:
        rows.append(py_rates(net, T, ndens, params, sup,
                             getattr(gen, "outside", "asym")))
    sup = [i for i in sup if all(math.isfinite(r[i]) for r in rows)]
    if not sup:
        return ('// make check: this network has no analytic rate to check\n'
                '#include <cstdio>\n#include "krome.hpp"\n'
                'int main() {\n  printf("make check: no analytic rate '
                'coefficient to check (photo, surface and tabulated rates '
                'only): OK\\n");\n  return 0;\n}\n')

    def fmt(v):
        return ",".join("%.17e" % x for x in v)
    return (CHECK_CPP.replace("@NT@", str(len(CHECK_T)))
            .replace("@NI@", str(len(sup)))
            .replace("@T@", fmt(CHECK_T))
            .replace("@N@", fmt(state[s] for s in sp))
            .replace("@I@", ",".join(str(i) for i in sup))
            .replace("@K@", ",".join("{%s}" % fmt(r[i] for i in sup)
                                     for r in rows))
            .replace("@P@", "%.17e" % PROBE_PARAM))


def build_probe(net, gen, tmp):
    hdr = tmp / ("net_%s.hpp" % net.name)
    hdr.write_text(gen.generate())
    TFUNCS.update(reverse_tfuncs(gen))
    TFUNCS.update(gfe_tfuncs(gen))
    src = tmp / ("probe_%s.cpp" % net.name)
    src.write_text(PROBE.replace("@HEADER@", hdr.name))
    exe = tmp / ("probe_" + net.name)
    r = subprocess.run([compiler(), "-O2", "-std=c++17", "-ffp-contract=off",
                        "-I", str(RUNTIME), str(src), "-o", str(exe)],
                       cwd=tmp, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode:
        raise RuntimeError("the probe did not compile:\n" + r.stderr[-3000:])
    return exe


def run_probe(exe, Tgas, sp_sorted, ndens):
    out = subprocess.run(
        [str(exe), "%.17e" % Tgas] + ["%.17e" % ndens[s] for s in sp_sorted],
        stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout
    res = {"k": {}, "d": {}, "s": {}, "J": {}, "f": {}}
    for line in out.splitlines():
        tag, i, v = line.split()
        res[tag][int(i)] = float(v)
    return res


def check(net, gen, tmp):
    """run the checks on a parsed network and its Codegen; returns a dict:
    ok (no failure), the worst code differences, n_stub, failures and
    warnings (lists of messages)"""
    exe = build_probe(net, gen, Path(tmp))
    sp_sorted = sorted(net.species)
    params = {p: PROBE_PARAM for p in gen.pnames}
    params["Tdust"] = 10.0
    rng = random.Random(42)

    pair = getattr(gen, "pair_gas", {})  # --gas-total: surface -> gas
    supported = independent(net, gen)
    expr = [rx for rx in net.reactions if rx["rate"]["kind"] == "expr"]
    with_param = {rx["id"] - 1 for rx in expr if _uses_param(gen, rx)}

    # FD columns are expensive (2 process runs each); sample for big nets
    N = len(sp_sorted)
    fd_cols = list(range(N)) if N <= 40 else rng.sample(range(N), 24)
    fd_ts = T_SAMPLES if N <= 40 else T_SAMPLES[1::3]

    worst_k = worst_d = worst_j = worst_f = 0.0
    worst_j_info = worst_k_info = None
    bad_rate, big_rate = {}, {}  # reaction index -> first offending (T, k)
    for Tgas in T_SAMPLES:
        ndens = {s: 10 ** rng.uniform(-8, 4) for s in sp_sorted}
        # the solver state: with --gas-total the surface slot of a gas /
        # surface pair holds the total (gas + surface)
        state = dict(ndens)
        for s_, g in pair.items():
            state[s_] = ndens[s_] + ndens[g]
        res = run_probe(exe, Tgas, sp_sorted, state)
        ck, cd, cs, cj, cf = (res[t] for t in "kdsJf")
        # freeze directives: the frozen rows the solver will zero
        pf = sorted({sp_sorted.index(s) for fr in net.freeze
                     if fr["cond"] is None
                     or eval_ast(fr["cond"], {"params": params})
                     for s in fr["species"]})
        if sorted({int(v) for v in cf.values()}) != pf:
            worst_f = 1.0

        # the rates themselves
        for rx in expr:
            i = rx["id"] - 1
            if not _in_range(rx, Tgas):
                continue
            v = ck[i]
            if not math.isfinite(v) or v < 0:
                bad_rate.setdefault(i, []).append((Tgas, v))
            o = _order(rx)
            if o in MAX_RATE and math.isfinite(v) and v > MAX_RATE[o] \
                    and i not in with_param:
                big_rate.setdefault(i, []).append((Tgas, v))

        pk = py_rates(net, Tgas, ndens, params, supported,
                      getattr(gen, "outside", "asym"))
        for i in supported:
            if math.isfinite(ck[i]) or math.isfinite(pk[i]):
                e = relerr(ck[i], pk[i])
                if e > worst_k:
                    worst_k, worst_k_info = e, (Tgas, i, ck[i], pk[i])
        # rhs stoichiometry test uses the C++ k values, so it is meaningful
        # even for reactions whose rate is stubbed
        pd = py_rhs(net, sp_sorted, ndens,
                    [ck[i] for i in range(len(net.reactions))])
        if pair:  # d(total)/dt = d(surface)/dt + d(gas)/dt
            ix = {s_: q for q, s_ in enumerate(sp_sorted)}
            pd = [pd[q] + (pd[ix[pair[s_]]] if s_ in pair else 0.0)
                  for q, s_ in enumerate(sp_sorted)]
        scale = max((abs(x) for x in pd if math.isfinite(x)), default=1.0) \
            or 1.0
        for i, s in enumerate(sp_sorted):
            if math.isfinite(cd[i]) and math.isfinite(pd[i]):
                worst_d = max(worst_d, abs(cd[i] - pd[i]) / scale)

        # analytic J vs central finite differences of the C++ rhs itself
        if Tgas not in fd_ts:
            continue
        for j in fd_cols:
            sj = sp_sorted[j]
            h = max(state[sj] * 1e-7, 1e-30)
            fd = []
            for sign in (+1, -1):
                nd2 = dict(state)
                nd2[sj] += sign * h
                fd.append(run_probe(exe, Tgas, sp_sorted, nd2)["d"])
            for i in range(N):
                fdv = (fd[0][i] - fd[1][i]) / (2 * h)
                an = cj[i * N + j]
                if not (math.isfinite(fdv) and math.isfinite(an)):
                    continue  # reported by the rate check
                # central differences cannot resolve entries whose effect on
                # dn_i is below roundoff of the |flux| scale over the step h
                # (dn_i itself may be a near-cancelled sum of huge fluxes)
                noise = 1e-10 * max(abs(cd[i]), cs[i]) / h
                denom = max(abs(an), abs(fdv), noise, 1e-30)
                rel = abs(an - fdv) / denom
                if rel > worst_j:
                    worst_j = rel
                    worst_j_info = (Tgas, sp_sorted[i], sj, an, fdv)

    failures, warnings = [], []
    if worst_k >= 1e-9:
        T, i, a, b = worst_k_info
        rx = net.reactions[i]
        failures.append("generated rates differ from the network's "
                        "expressions by %.1e (line %s, %s at T = %g K: "
                        "%.6e vs %.6e)" % (worst_k, rx.get("line", "?"),
                                           rx["eq"], T, a, b))
    if worst_d >= 1e-12:
        failures.append("generated rhs differs from the stoichiometry by "
                        "%.1e" % worst_d)
    if worst_j >= 1e-4:
        T, si, sj2, an, fdv = worst_j_info
        failures.append("analytic Jacobian differs from finite differences "
                        "by %.1e (T=%g, d(dn_%s)/dn_%s: %.6e vs %.6e)"
                        % (worst_j, T, si, sj2, an, fdv))
    if worst_f:
        failures.append("the frozen species differ from the freeze "
                        "directives")
    def where(rx, pts):
        ts = ", ".join("%g" % T for T, _ in pts)
        hint = ""
        if not rx["trange"] or math.isinf(rx["trange"][1]):
            hint = " (no upper T limit is declared: a fit used beyond its " \
                   "validity? declare its T range)"
        return "at T = %s K%s" % (ts, hint)
    for i, pts in sorted(bad_rate.items()):
        rx = net.reactions[i]
        head = "line %s, %s: " % (rx.get("line", "?"), rx["eq"])
        nonfin = [(T, v) for T, v in pts if not math.isfinite(v)]
        if nonfin and i not in with_param:
            failures.append(head + "rate %g %s" % (nonfin[0][1],
                                                  where(rx, nonfin)))
        elif nonfin:
            warnings.append(head + "rate %g %s (with the parameters at the "
                            "probe value %g)" % (nonfin[0][1],
                                                 where(rx, nonfin),
                                                 PROBE_PARAM))
        neg = [(T, v) for T, v in pts if math.isfinite(v)]
        if neg:
            warnings.append(head + "negative rate (%.3g) %s" % (
                neg[0][1], where(rx, neg)))
    for i, pts in sorted(big_rate.items()):
        rx = net.reactions[i]
        o = _order(rx)
        warnings.append("line %s, %s: rate up to %.3g cm^%d s^-1, above %g, "
                        "%s; or a typo in an exponent or the units?" % (
                            rx.get("line", "?"), rx["eq"],
                            max(v for _, v in pts), 3 * (o - 1),
                            MAX_RATE[o], where(rx, pts)))
    rw, rn = repeated(net)
    bw, bn = balance_warnings(net, gen)
    warnings += rw + budget_warnings(net) + bw
    rn = rn + bn
    return {"ok": not failures, "k": worst_k, "rhs": worst_d,
            "jac": worst_j, "freeze": worst_f,
            "n_stub": len(net.reactions) - len(supported),
            "failures": failures, "warnings": warnings, "notes": rn}


def _say(tag, msg):
    import textwrap
    print(textwrap.fill(msg, 79, initial_indent="  %-9s" % tag,
                        subsequent_indent=" " * 11, break_long_words=False,
                        break_on_hyphens=False))


def report(net, res, gen_warnings=()):
    """the user-facing report of check() (python3 -m kromex --check)"""
    print("network %s: %d species, %d reactions"
          % (net.name, len(net.species), len(net.reactions)))
    print("  generated code vs an independent evaluation of the network "
          "(%d temperatures, %g-%g K):" % (len(T_SAMPLES), T_SAMPLES[0],
                                           T_SAMPLES[-1]))
    print("    rates     %.1e  (max relative difference; %d reactions not "
          "evaluated independently:\n               photo, tables and "
          "external functions)" % (res["k"], res["n_stub"]))
    print("    rhs       %.1e  (from the stoichiometry)" % res["rhs"])
    print("    Jacobian  %.1e  (analytic vs finite differences)" % res["jac"])
    for f in res["failures"]:
        _say("FAIL", f)
    for w in list(gen_warnings) + res["warnings"]:
        _say("WARNING", w)
    for n in res["notes"] + [n for n in getattr(net, "notices", [])
                             if "additive channels" not in n]:
        _say("note", n)
    print("RESULT: %s (%d failures, %d warnings)" % (
        "OK" if res["ok"] else "FAILED", len(res["failures"]),
        len(res["warnings"]) + len(gen_warnings)))
