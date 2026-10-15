#!/usr/bin/env python3
"""A KIDA network (kida.uva: gas_reactions_*.in) as a KROME-X network
(.knet):

    kida2krome gas_reactions_kida.uva.2024.in kida2024.knet
    (from a checkout: python3 -m kromex_convert.kida2krome ...)

The file format (kida.uva.2024_readme.txt): 3(a11),1x,5(a11),1x,
3(e10.3,1x),e8.2,1x,e8.2,1x,a4,1x,i2,1x,i6,1x,i6,1x,i2,1x,i5,i2,i3 =
reactants, products, alpha, beta, gamma, F, g, type of uncertainty,
itype, Tmin, Tmax, formula, number, number of (alpha, beta, gamma),
recommendation. The formulas (Wakelam et al. 2012, ApJS 199, 21; 2024,
A&A 689, A63):
 1 cosmic rays (CR, CRP: direct, and the UV they induce): alpha zeta
   [s^-1], zeta
   the H2 ionisation rate (the parameter crflux, default 1.3e-17 s^-1);
 2 photodissociation in the Draine (1978) field: alpha exp(-gamma Av)
   [s^-1], times the UV scale (--uv-param, the field in Draine units);
 3 Kooij: alpha (T/300)^beta exp(-gamma/T) [cm^3 s^-1];
 4 ionpol1: alpha beta (0.62 + 0.4767 gamma (300/T)^0.5);
 5 ionpol2: alpha beta (1 + 0.0967 gamma (300/T)^0.5 + gamma^2/10.526
   (300/T)), alpha the branching ratio, beta the Langevin rate, gamma the
   x factor at 300 K.
Other formulas (0: grain charging, 10-11: Nautilus's surface hydrogen)
are skipped with a comment. For one-body reactions the pseudo-reactant
(CR, CRP, Photon) decides between 1 and 2: in KIDA 2024 the formula
column disagrees with it in 25 lines (a warning for each).

Temperature ranges: -9999, 9999 (not relevant) is none; 0, 800 is KIDA's
default for no information, kept as given; a range of one temperature is
none. Lines of one reaction whose ranges tile are one piecewise rate;
overlapping ones are KIDA's alternative values for the same range (its
readme: the user chooses one): the best recommended is kept (then the
first), the others are comments. Pseudo-species (CR, CRP, Photon) leave
the reaction."""

import argparse

from . import ratedb

PSEUDO = {"CR", "CRP", "PHOTON", ""}
# the Fortran format's columns
FIELDS = ([("r%d" % i, 11) for i in range(3)] + [("", 1)]
          + [("p%d" % i, 11) for i in range(5)] + [("", 1)]
          + [("a", 10), ("", 1), ("b", 10), ("", 1), ("c", 10), ("", 1),
             ("F", 8), ("", 1), ("g", 8), ("", 1), ("unc", 4), ("", 1),
             ("itype", 2), ("", 1), ("tmin", 6), ("", 1), ("tmax", 6),
             ("", 1), ("formula", 2), ("", 1), ("num", 5), ("nabg", 2),
             ("recom", 3)])


def parse_line(line):
    out, pos = {}, 0
    for name, wd in FIELDS:
        if name:
            out[name] = line[pos:pos + wd].strip()
        pos += wd
    return out


def kind_of(d):
    """the formula to use: the pseudo-reactant decides for one-body
    reactions (CR, CRP: 1; Photon: 2), the formula column for the others.
    KIDA 2024 has 25 one-body lines whose formula disagrees: CRP with
    formula 2 (alpha 1e3-1e4, gamma 0: cosmic-ray factors), Photon with
    formula 1 (alpha ~1e-9, gamma ~2: UV photorates). Returns (formula,
    note or None)"""
    pseudo = {d["r%d" % i].upper() for i in range(3)} & (PSEUDO - {""})
    f = d["formula"]
    want = 1 if pseudo & {"CR", "CRP"} else 2 if "PHOTON" in pseudo else f
    if want != f:
        return want, ("formula %d with %s: used as formula %d"
                      % (f, "/".join(sorted(pseudo)), want))
    return f, None


def rate_of(d, cr_param, uv_param):
    f = kind_of(d)[0]
    a, b, c = d["a"], d["b"], d["c"]

    def fac(x):     # a factor: a negative one in parentheses
        return ratedb.num(x) if float(x) >= 0 else "(%s)" % ratedb.num(x)
    if f == 1:
        return "%s*%s" % (fac(a), cr_param)
    if f == 2:
        s = "%s*%s" % (fac(a), uv_param)
        return s + ("*exp(-(%s)*Av)" % ratedb.num(c) if float(c) else "")
    if f == 3:
        return ratedb.kooij(a, b, c)
    if f == 4:
        return "%s*%s*(0.62 + 0.4767*%s*sqrt(300/Tgas))" % (
            fac(a), fac(b), fac(c))
    if f == 5:
        return ("%s*%s*(1 + 0.0967*%s*sqrt(300/Tgas) + %s**2/10.526*"
                "300/Tgas)" % (fac(a), fac(b), fac(c), fac(c)))
    return None


def suspicious(d):
    """parameters that do not fit their formula (KIDA 2024 has a few:
    ionpol parameters under Kooij, Kooij ones under ionpol)"""
    f, a, b = kind_of(d)[0], float(d["a"]), float(d["b"])
    if f == 3 and a >= 0.1 and 0.0 < b < 1e-6:
        return ("formula 3 (Kooij) with ion-polar-like parameters (alpha "
                "%g, beta %g): k = %g cm^3 s^-1" % (a, b, a))
    if f in (4, 5) and (b <= 0.0 or a < 1e-3):
        return ("formula %d (ionpol) with Kooij-like parameters (alpha %g, "
                "beta %g): k %s" % (f, a, b, "< 0" if a * b < 0 else
                                    "~ alpha beta = %g" % (a * b)))
    return None


def convert(src, cr_param="crflux", uv_param="uv_scale", keep_ranges=True):
    w = ratedb.Writer(src, "KIDA kida.uva format", [
        (cr_param, 1.3e-17, "cosmic-ray ionisation rate of H2 zeta [s^-1] "
         "(KIDA: k = alpha zeta)"),
        (uv_param, 1.0, "the UV field in units of Draine's (1978) (KIDA "
         "photorates)")])
    groups = {}
    order = []
    raw_lines = open(src).read().splitlines()
    for line in raw_lines:
        if not line.strip() or line.lstrip().startswith(("!", "#")):
            continue
        d = parse_line(line)
        try:
            d["formula"] = int(d["formula"])
            d["tmin_v"], d["tmax_v"] = float(d["tmin"]), float(d["tmax"])
            d["recom_v"] = int(d["recom"] or 1)
            float(d["a"]), float(d["b"]), float(d["c"])
        except ValueError:
            w.skip(line, "unreadable line")
            continue
        d["line"] = line
        reac = [ratedb.species(d["r%d" % i]) for i in range(3)
                if d["r%d" % i].upper() not in PSEUDO]
        prod = [ratedb.species(d["p%d" % i]) for i in range(5)
                if d["p%d" % i].upper() not in PSEUDO]
        # the process too: the cosmic-ray-induced and the UV dissociation
        # of a species (CRP, Photon) are separate channels that add up
        process = tuple(sorted({d["r%d" % i].upper() for i in range(3)}
                               & (PSEUDO - {""})))
        key = (tuple(sorted(reac)), tuple(sorted(prod)), process)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((reac, prod, d))
    for key in order:
        entries = groups[key]
        reac, prod = entries[0][0], entries[0][1]
        ok = []
        for r, p, d in entries:
            if rate_of(d, cr_param, uv_param) is None:
                w.skip(d["line"], "formula %d not supported" % d["formula"])
            else:
                ok.append(d)
        if not ok:
            continue
        # the ranges: none (-9999, 9999), one temperature (none), or as
        # given
        def rng(d):
            """the range: none for -9999, 9999 (not relevant); else the
            bounds as given (9999 alone is a temperature: H + Cl has 3000-
            9999 then 10000-20000 K; a one-temperature range is a piece
            here, and no range for a reaction of one line)"""
            lo, hi = d["tmin_v"], d["tmax_v"]
            if (lo <= -9999 and hi >= 9999) or not keep_ranges:
                return None
            return (max(lo, 0.0), hi)
        # alternatives: entries whose ranges overlap; the pieces of one
        # rate: ranges that do not overlap
        chosen = []
        for d in sorted(ok, key=lambda d: (-d["recom_v"],
                                           int(d["num"] or 0))):
            r = rng(d)
            def overlap(p, q):
                # a one-temperature piece [T, T] overlaps only a range
                # strictly around T
                return p[0] < q[1] and q[0] < p[1] if p[0] < p[1] and \
                    q[0] < q[1] else (q[0] < p[0] < q[1] if q[0] < q[1]
                                      else p[0] < q[0] < p[1] if p[0] < p[1]
                                      else p[0] == q[0])
            clash = [c for c in chosen if r is None or rng(c) is None or
                     overlap(r, rng(c))]
            if clash:
                w.skip(d["line"], "an alternative to KIDA %s (recommendation "
                       "%d, kept: %d)" % (clash[0]["num"], d["recom_v"],
                                          clash[0]["recom_v"]))
                continue
            chosen.append(d)
        chosen.sort(key=lambda d: (rng(d) or (0.0, 0.0))[0])
        notes = []
        if len(chosen) > 1:
            pieces, notes = ratedb.tile([rng(d) + (d,) for d in chosen])
            chosen = [p[2] for p in pieces]
            ranges = [(p[0], p[1]) for p in pieces]
        else:
            r0 = rng(chosen[0])
            if r0 is not None and r0[0] == r0[1]:
                notes.append("measured at %s K: no range" % ratedb.num(r0[0]))
                r0 = None
            ranges = [r0]
        for i, (d, r) in enumerate(zip(chosen, ranges)):
            odd = suspicious(d) or kind_of(d)[1]
            if odd:
                w.suspicious.append("KIDA %s %s: %s" % (
                    d["num"], " + ".join(reac), odd))
            note = "KIDA %s (itype %s, formula %d%s)%s" % (
                d["num"], d["itype"], d["formula"],
                "; " + "; ".join(notes) if notes and i == 0 else "",
                "\nWARNING: %s; check the KIDA entry" % odd if odd else "")
            w.reaction(reac, prod, ratedb.trange(*r) if r else "[]",
                       rate_of(d, cr_param, uv_param), note=note)
    return w


def main(argv=None):
    ap = argparse.ArgumentParser(prog="kida2krome", description=__doc__.split("\n\n")[0])
    ap.add_argument("reactions", help="the kida.uva reaction file")
    ap.add_argument("output", help="the .knet to write")
    ap.add_argument("--cr-param", default="crflux",
                    help="the cosmic-ray rate parameter (default crflux)")
    ap.add_argument("--uv-param", default="uv_scale",
                    help="the UV scale parameter (default uv_scale)")
    ap.add_argument("--no-ranges", action="store_true",
                    help="no temperature ranges")
    a = ap.parse_args(argv)
    w = convert(a.reactions, a.cr_param, a.uv_param, not a.no_ranges)
    for x in w.suspicious:
        print("WARNING: %s" % x)
    net = w.write(a.output, [
        "KIDA: Wakelam et al. 2012, ApJS 199, 21; the kida.uva networks: "
        "Wakelam et al. 2015, ApJS 217, 20; 2024, A&A 689, A63",
        "rates outside their T range: KROME-X's default (rates.outside = "
        '"asym") gives 0 below the range;',
        'use rates.outside = "clamp" for gas colder than a range\'s lower '
        "end"])
    return 1 if net.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
