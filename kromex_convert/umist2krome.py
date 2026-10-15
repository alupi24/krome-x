#!/usr/bin/env python3
"""A UMIST Database for Astrochemistry ratefile (RATE22, RATE12: the
colon-separated format) as a KROME-X network (.knet):

    umist2krome rate22_final.rates rate22.knet
    (from a checkout: python3 -m kromex_convert.umist2krome ...)

The rates (Millar et al. 2024, A&A 682, A109, section 2.2):
- two-body (every type but CP, CR, PH): alpha (T/300)^beta exp(-gamma/T)
  [cm^3 s^-1];
- CP, cosmic-ray ionisation: alpha zeta/zeta0 [s^-1];
- CR, cosmic-ray-induced photoreactions: alpha (T/300)^beta gamma /
  (1 - omega) zeta/zeta0 [s^-1], omega the grain albedo in the far UV
  (0.5, as the database's models; --albedo);
- PH, photoreactions: alpha exp(-gamma Av) [s^-1], alpha the rate in the
  unshielded interstellar field; times the UV scale (--uv-param, default
  1: the database's field).
zeta0 = 1.36e-17 s^-1, the database's normalisation; the parameter
crflux is zeta (default zeta0: the database's rates).

Each line holds 18 fields (number, type, 2 reactants, 4 products, NTR,
alpha, beta, gamma, Tl, Tu, ST, ACC, REF, notes), and 9 more for each
further temperature range (NTR > 1): a piecewise rate. Separate lines of
one reaction are a sum of fits (RATE22, section 2.2): additive channels in
KROME-X. The DOI of each rate (REF) becomes its `ref`. Pseudo-species
(PHOTON, CRP, CRPHOT) leave the reaction."""

import argparse
import csv

from . import ratedb

ZETA0 = 1.36e-17
PSEUDO = {"PHOTON", "CRP", "CRPHOT", ""}


def convert(src, cr_param="crflux", uv_param="uv_scale", albedo=0.5,
            keep_ranges=True):
    w = ratedb.Writer(src, "UMIST RATE format", [
        (cr_param, ZETA0, "cosmic-ray ionisation rate zeta [s^-1]; the "
         "CP and CR rates scale as zeta / 1.36e-17"),
        (uv_param, 1.0, "the UV field in units of the database's "
         "unshielded interstellar field (PH rates)")])
    with open(src, newline="") as f:
        rows = list(csv.reader(f, delimiter=":", quotechar='"'))
    seen = {}
    for row in rows:
        raw = ":".join(row)
        if len(row) < 17 or not row[0].strip().isdigit():
            if raw.strip():
                w.skip(raw, "not a rate line")
            continue
        typ = row[1].strip()
        reac = [ratedb.species(x) for x in row[2:4]
                if x.strip() not in PSEUDO]
        prod = [ratedb.species(x) for x in row[4:8]
                if x.strip() not in PSEUDO]
        ntr = int(row[8])
        pieces = []
        for k in range(ntr):
            b = 9 + 9 * k
            if len(row) < b + 8:
                break
            a, be, g = row[b], row[b + 1], row[b + 2]
            tl, tu = float(row[b + 3]), float(row[b + 4])
            ref = row[b + 7].strip()
            pieces.append((tl, tu, a, be, g, ref))
        if len(pieces) != ntr:
            w.skip(raw, "fewer ranges than NTR = %d" % ntr)
            continue
        # a range repeated in one line (UMIST 7246) is one range
        uniq = []
        for p in pieces:
            if not any(p[:2] == q[:2] and [float(x) for x in p[2:5]] ==
                       [float(x) for x in q[2:5]] for q in uniq):
                uniq.append(p)
        pieces, ntr = uniq, len(uniq)
        # the same rate written twice (same reaction and parameters,
        # whatever the ranges: UMIST 1452 and 2061) would be summed as a
        # second channel: the first is kept
        key = (tuple(sorted(reac)), tuple(sorted(prod)), typ,
               tuple(tuple(float(x) for x in p[2:5]) for p in pieces))
        if key in seen:
            w.skip(raw, "the rate of UMIST %s again" % seen[key])
            continue
        seen[key] = row[0].strip()
        pieces, notes = ratedb.tile(pieces)
        if ntr == 1:
            rng1, n1 = ratedb.single_t(pieces[0][0], pieces[0][1])
            if n1:
                notes.append(n1)
        for i, (tl, tu, a, be, g, ref) in enumerate(pieces):
            if typ == "CP":
                rate = "%s*%s/%s" % (ratedb.num(a), cr_param,
                                     ratedb.num(ZETA0))
            elif typ == "CR":
                rate = "%s*%s/(1-%s)*%s/%s" % (
                    ratedb.kooij(a, be, 0.0), ratedb.num(g),
                    ratedb.num(albedo), cr_param, ratedb.num(ZETA0))
            elif typ == "PH":
                rate = "%s*%s" % (ratedb.num(a), uv_param)
                if float(g) != 0.0:
                    rate += "*exp(-(%s)*Av)" % ratedb.num(g)
            else:
                rate = ratedb.kooij(a, be, g)
            rng = ratedb.trange(tl, tu) if keep_ranges and (
                ntr > 1 or typ not in ("CP", "CR", "PH")) else "[]"
            if ntr == 1 and rng1:
                rng = rng1
            note = "UMIST %s %s (%s)" % (row[0].strip(), typ,
                                         "; ".join(notes) if notes and i == 0
                                         else "range %d of %d" % (i + 1, ntr)
                                         if ntr > 1 else "")
            w.reaction(reac, prod, rng, rate,
                       ref="doi:" + ref if ref and not ref.startswith(
                           "http") else ref or None,
                       note=note.replace(" ()", ""))
    return w


def main(argv=None):
    ap = argparse.ArgumentParser(prog="umist2krome", description=__doc__.split("\n\n")[0])
    ap.add_argument("ratefile")
    ap.add_argument("output", help="the .knet to write")
    ap.add_argument("--cr-param", default="crflux",
                    help="the cosmic-ray rate parameter (default crflux)")
    ap.add_argument("--uv-param", default="uv_scale",
                    help="the UV scale parameter (default uv_scale)")
    ap.add_argument("--albedo", type=float, default=0.5,
                    help="the grain albedo omega of the CR rates (default "
                         "0.5, the database's models)")
    ap.add_argument("--no-ranges", action="store_true",
                    help="no temperature ranges (KROME-X: k = 0 below a "
                         "range unless rates.outside = clamp)")
    a = ap.parse_args(argv)
    w = convert(a.ratefile, a.cr_param, a.uv_param, a.albedo,
                not a.no_ranges)
    net = w.write(a.output, [
        "UMIST Database for Astrochemistry: Millar et al. 2024, A&A 682, "
        "A109 (RATE22); McElroy et al. 2013, A&A 550, A36 (RATE12)",
        "rates outside their T range: KROME-X's default (rates.outside = "
        '"asym") gives 0 below the range;',
        'use rates.outside = "clamp" for gas colder than a range\'s lower '
        "end (the database's Tl is often 10 K)"])
    return 1 if net.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
