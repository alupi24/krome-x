"""Detailed balance of a network: the reverse reactions it may need.

A reaction that runs one way only cannot reach chemical equilibrium: where
its reverse is fast, the network pumps the species one way (a hot gas that
keeps its molecules, ions that never recombine three-body at high
density). In equilibrium the reverse process of the same collision has
the rate k_rev = k_fwd / K_c(T) (detailed balance), so the ratio of the
reverse to the forward flux, with every species at the density n, is

    B / F = n^dn / K_c(T),   1 / K_c = exp(dG / RT) (kB T / 1 bar)^dn,

dG = sum G(products) - sum G(reactants), dn = n_products - n_reactants:
thermodynamics alone, the rate coefficients cancel. Where B / F reaches a
threshold (1 % by default) inside the reaction's temperature range, the
reverse is not negligible there; for dn != 0 that temperature depends on
the density (dissociation and ionization win at low density).

The reverse is a collision of the same kind only for collisional
processes: two-body exchanges (dn = 0), three-body associations (a third
body among the reactants, or a rate proportional to ntot: the reverse is
collisional dissociation) and fragmentations (dn > 0: the reverse is a
three-body association or recombination). The reverse of a radiative
association or recombination (two reactants, one product, no third body)
is a photo-process set by the radiation field, not by detailed balance
with the gas: those are listed apart and never flagged; photo, surface,
tabulated and nuclear reactions are not analysed.

G/RT of a species: Burcat's NASA 7-coefficient polynomials
(data/thermo_burcat.toml), else NIST-JANAF Shomate (data/thermoNIST.dat),
as the code generator's `<->` reverse rates, extended beyond the data
with a constant heat capacity (gRT); the electron, which the tables leave
out, from the Sackur-Tetrode entropy of an ideal gas of degeneracy 2
(H(298.15 K) = 0, the tables' electron convention); ortho and para H2
from H2 and their rotational partition sums (the tables' records of the
spin-resolved species are the normal species': the other spin-resolved
species have no data here).

Used by python3 -m kromex --check (netcheck: a warning per flagged
reaction) and kromex_doc (the page of the reverse reactions)."""

import math
from pathlib import Path

import constants

P0 = 1e6                     # the standard state, 1 bar [dyn cm^-2]
THRESHOLD = 1e-2             # B / F that makes a reverse not negligible
DENSITIES = (1.0, 1e4, 1e8, 1e12)   # [cm^-3]
TGRID = (10.0, 1e5, 81)      # the temperatures examined: lo, hi, points


# ---- thermochemistry ----------------------------------------------------

def _find(data_dir, name):
    for base in (data_dir.parent, data_dir, data_dir.parent.parent):
        if (base / name).exists():
            return base / name
    return None


def load_tables(data_dir):
    """(nasa, nist): Burcat's polynomials, species name -> [Tlow, Tcommon,
    Thigh, 7 upper-range, 7 lower-range coefficients]; NIST Shomate, NAME
    (upper case) -> [(Tmin, Tmax, [A..G])], sorted (the code generator's
    load_thermo / load_thermo_nist)"""
    import tomlio
    data_dir = Path(data_dir)
    nasa, nist = {}, {}
    path = _find(data_dir, "thermo_burcat.toml")
    if path is not None:
        sp = tomlio.load(path).get("species", {})
        nasa = {k: list(v["T"]) + list(v["high"]) + list(v["low"])
                for k, v in sp.items()}
    path = _find(data_dir, "thermoNIST.dat")
    if path is not None:
        for row in path.read_text().splitlines():
            s = row.strip()
            if not s or s[0] == "#":
                continue
            a = s.split(":")
            v = [float(x) for x in a[1:]]
            nist.setdefault(a[0].upper(), []).append((v[0], v[1], v[2:9]))
        for k in nist:
            nist[k].sort()
    return nasa, nist


def electron_S0():
    """the electron's S/R - 5/2 ln T at 1 bar (Sackur-Tetrode, g = 2): the
    a7 of its NASA record (-11.72)"""
    m = constants.e_mass
    kB = constants.boltzmann_erg
    h = constants.planck_erg
    return (math.log(2.0 * (2.0 * math.pi * m * kB / h ** 2) ** 1.5
                     * kB / P0) + 2.5)


_E_S0 = electron_S0()


# spin-resolved H2: the tables' H2_ortho and H2_para records are normal H2
# (the same numbers), so their G/RT is H2's plus ln(q_H2 / q_spin), the
# rotational partition sums of a rigid rotor (B, as the code generator's
# adiabatic index): para the even J, ortho 3 x the odd J (nuclear spin)
H2_BROT = 85.4    # K
_SPIN_H2 = {"H2_ortho": 1, "H2_para": 0}


def _h2_spin_lnratio(sp, T):
    """ln(q_H2 / q_spin) of H2_ortho or H2_para at T"""
    jmax = int(30 + 3 * math.sqrt(T / H2_BROT))
    qe = sum((2 * j + 1) * math.exp(-H2_BROT * j * (j + 1) / T)
             for j in range(0, jmax, 2))
    qo = 3.0 * sum((2 * j + 1) * math.exp(-H2_BROT * j * (j + 1) / T)
                   for j in range(1, jmax, 2))
    return math.log((qe + qo) / (qo if _SPIN_H2[sp] else qe))


def _spin_resolved(sp):
    return "_ortho" in sp or "_para" in sp


def has_data(tables, sp):
    """thermochemical data of sp: the electron, a record, or H2_ortho /
    H2_para (from H2); the other spin-resolved species (H3_ortho+, ...)
    have none (their records are those of the normal species)"""
    nasa, nist = tables
    if sp in _SPIN_H2:
        return "H2" in nasa or "H2" in nist
    if _spin_resolved(sp):
        return False
    return sp == "e-" or sp in nasa or sp.upper() in nist


def _nasa_HS(a, T):
    """H/RT, S/R, cp/R of the NASA polynomial a at T"""
    H = (a[0] + a[1] * T / 2 + a[2] * T * T / 3 + a[3] * T ** 3 / 4
         + a[4] * T ** 4 / 5 + a[5] / T)
    S = (a[0] * math.log(T) + a[1] * T + a[2] * T * T / 2
         + a[3] * T ** 3 / 3 + a[4] * T ** 4 / 4 + a[6])
    cp = a[0] + a[1] * T + a[2] * T * T + a[3] * T ** 3 + a[4] * T ** 4
    return H, S, cp


def _shomate_HS(a, T):
    """H/RT, S/R, cp/R of the Shomate coefficients a at T"""
    t = T * 1e-3
    H = (a[0] * t + a[1] * t * t / 2 + a[2] * t ** 3 / 3
         + a[3] * t ** 4 / 4 - a[4] / t + a[5])
    S = (a[0] * math.log(t) + a[1] * t + a[2] * t * t / 2
         + a[3] * t ** 3 / 3 - a[4] / (2 * t * t) + a[6])
    cp = a[0] + a[1] * t + a[2] * t * t + a[3] * t ** 3 + a[4] / (t * t)
    return (H / (constants.Rgas_kJ * T), S / constants.Rgas_J,
            cp / constants.Rgas_J)


def gRT(tables, sp, T):
    """G/RT = H/RT - S/R of species sp at T (1 bar); None without data.
    Outside the data, with the heat capacity of its nearest end held
    constant (H + cp (T - T0), S + cp ln(T / T0): exact for atoms and
    ions, cp = 5/2 R): the tables of most species stop at a few thousand
    K, and holding G/RT at the end instead would miss the entropy that
    makes ionization and dissociation win at high T"""
    nasa, nist = tables
    if sp == "e-":     # H = 5/2 R (T - 298.15), S/R = 5/2 ln T + S0
        return 2.5 * (1.0 - 298.15 / T) - (2.5 * math.log(T) + _E_S0)
    if sp in _SPIN_H2:
        g = gRT(tables, "H2", T)
        return None if g is None else g + _h2_spin_lnratio(sp, T)
    if _spin_resolved(sp):
        return None
    c = nasa.get(sp)
    if c is not None:
        T0 = min(max(T, c[0]), c[2])

        def f(x):
            return _nasa_HS(c[10:17] if x <= c[1] else c[3:10], x)
    else:
        r = nist.get(sp.upper())
        if r is None:
            return None
        T0 = min(max(T, r[0][0]), r[-1][1])

        def f(x):
            return _shomate_HS(next((q[2] for q in r if x <= q[1]),
                                    r[-1][2]), x)
    H0, S0, cp = f(T0)
    if T0 == T:
        return H0 - S0
    # H/RT = (H0 T0 + cp (T - T0)) / T in units of R
    return (H0 * T0 + cp * (T - T0)) / T - (S0 + cp * math.log(T / T0))


def inv_kc(tables, reactants, products, T):
    """1 / K_c(T) = exp(dG/RT) (kB T / 1 bar)^dn [cm^3^dn] of reactants ->
    products, the factor of the reverse rate k_rev = k_fwd / K_c; None
    without data"""
    g = [gRT(tables, x, T) for x in list(reactants) + list(products)]
    if None in g:
        return None
    nr = len(reactants)
    dg = sum(g[nr:]) - sum(g[:nr])
    dn = len(products) - nr
    return (math.exp(min(dg, 700.0))
            * (constants.boltzmann_erg * T / P0) ** dn)


def flux_ratio(tables, reactants, products, T, n):
    """B / F: the reverse over the forward flux of the reaction in
    equilibrium, every species at the density n [cm^-3]; None without
    data"""
    ik = inv_kc(tables, reactants, products, T)
    if ik is None:
        return None
    return ik * n ** (len(products) - len(reactants))


# ---- the analysis -------------------------------------------------------

def _uses_ntot(ast):
    stack = [ast]
    while stack:
        a = stack.pop()
        if isinstance(a, dict):
            if a.get("s") == "ntot":
                return True
            stack.extend(a.values())
        elif isinstance(a, list):
            stack.extend(a)
    return False


def kind_of(rx):
    """the process of a reaction, for its reverse: "exchange" (two-body,
    dn = 0), "association" (three-body: a third body among the reactants
    or a rate proportional to ntot), "fragmentation" (two-body, dn > 0),
    "radiative" (two reactants, one product, no third body: the reverse is
    a photo-process), "driven" (one reactant: cosmic-ray, photo or decay
    processes written as expressions, not collisions of the gas), or None
    (not analysed: photo, surface, tabulated, nuclear)"""
    if rx["rate"]["kind"] != "expr" or rx.get("nuclear"):
        return None
    nr, np_ = len(rx["reactants"]), len(rx["products"])
    if nr == 1:
        return "driven"
    if nr >= 3 or _uses_ntot(rx["rate"]["ast"]):
        return "association"
    if np_ > nr:
        return "fragmentation"
    if np_ == nr:
        return "exchange"
    return "radiative"


def _key(side):
    return tuple(sorted(x["sp"] for x in side))


def _collisional(rx):
    """a collision of the gas: two reactants or more, or one with a rate
    proportional to ntot (a partner implied)"""
    if len(rx["reactants"]) >= 2:
        return True
    return rx["rate"]["kind"] == "expr" and _uses_ntot(rx["rate"]["ast"])


def analyse(net, tables, tgrid=TGRID, densities=DENSITIES,
            threshold=THRESHOLD):
    """the reverse reactions of net, a reaction's temperature pieces taken
    together: {"flagged", "endothermic", "radiative", "present", "nodata":
    [entry], "driven", "negligible": counts}. Entries {"eq", "line", "ids"
    (its pieces), "R", "P" (the species), "kind", "dn", "range": (Tlo,
    Thi)}; a reaction is in
      present      its reverse is in the network: `<->`, or written out
                   ("reverse_ids": that reaction's pieces);
      nodata       a species has no thermochemical data ("missing");
      radiative    the reverse is a photo-process (not flagged);
      endothermic  B / F >= threshold from the lowest temperature of its
                   range at every density: the reverse is the dominant
                   direction throughout (an exothermic association or
                   recombination the network lacks; with dn > 0 a
                   three-body one, important at high density);
      flagged      the reverse grows to B / F >= threshold inside its
                   range: "cross" {n: the lowest such T [K], or None};
    driven (one reactant) and negligible (B / F below threshold
    throughout) are counted"""
    lo, hi, npt = tgrid
    grid = [lo * (hi / lo) ** (i / (npt - 1)) for i in range(npt)]
    # the reverses that count: collisions (a one-reactant reaction is a
    # cosmic-ray, photo or decay process, not the collisional reverse,
    # unless its rate is proportional to ntot: dissociation by a partner)
    by_key = {}
    for rx in net.reactions:
        by_key.setdefault((_key(rx["reactants"]), _key(rx["products"])),
                          []).append(rx)
    keys = {k for k, rxs in by_key.items()
            if any(_collisional(rx) for rx in rxs)}
    out = {"flagged": [], "endothermic": [], "radiative": [],
           "present": [], "nodata": [], "driven": 0, "negligible": 0}
    seen = set()
    for rx in net.reactions:
        k = kind_of(rx)
        if k is None:
            continue
        if k == "driven":
            out["driven"] += 1
            continue
        R, P = _key(rx["reactants"]), _key(rx["products"])
        if (R, P) in seen:     # the pieces of a piecewise rate, together
            continue
        seen.add((R, P))
        pieces = [r for r in by_key[(R, P)] if kind_of(r) is not None]
        tlo = min((r["trange"] or (0.0, math.inf))[0] for r in pieces)
        thi = max((r["trange"] or (0.0, math.inf))[1] for r in pieces)
        entry = {"eq": rx["eq"], "line": rx["line"],
                 "ids": [r["id"] for r in pieces], "R": R, "P": P,
                 "kind": k, "dn": len(P) - len(R), "range": (tlo, thi)}
        if any(r.get("reverse") for r in pieces) or (P, R) in keys:
            entry["reverse_ids"] = [r["id"] for r in by_key.get((P, R), [])
                                    if _collisional(r)]
            out["present"].append(entry)
            continue
        missing = [x for x in R + P if not has_data(tables, x)]
        if missing:
            entry["missing"] = sorted(set(missing))
            out["nodata"].append(entry)
            continue
        if k == "radiative":
            out["radiative"].append(entry)
            continue
        Ts = [T for T in grid if tlo <= T <= thi]
        if not Ts:
            out["negligible"] += 1
            continue
        # B/F at n = B/F at 1 cm^-3 times n^dn: one evaluation per T
        r1 = [flux_ratio(tables, R, P, T, 1.0) for T in Ts]
        dn = entry["dn"]
        cross = {n: next((T for T, r in zip(Ts, r1)
                          if r * n ** dn >= threshold), None)
                 for n in densities}
        entry["cross"] = cross
        if all(v == Ts[0] for v in cross.values()):
            out["endothermic"].append(entry)
        elif any(v is not None for v in cross.values()):
            out["flagged"].append(entry)
        else:
            out["negligible"] += 1
    return out
