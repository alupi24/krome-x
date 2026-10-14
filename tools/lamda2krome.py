#!/usr/bin/env python3
"""LAMDA molecular data file -> KROME coolant file (the coolZ.dat format
read by codegen thermal.coolant_files), a port of legacy tools/lamda2.py.

    python3 tools/lamda2krome.py o2.dat coolO2.dat [--name O2]

Levels (energies [K], degeneracies), radiative transitions (up -> low,
Aij) and, per collision partner, the de-excitation rates as cubic-spline
tables in Tgas (fspline, legacy's Godunov spline, the end values outside
the table). Differences from legacy lamda2.py:
- energies: E [cm^-1] x hc/k = 1.4387769 K cm (legacy 1.42879, 0.7 % low;
  the legacy data/coolO2.dat carries that error);
- the file is parsed by its '!' section headers instead of fixed line
  positions, and LAMDA partner 7 (H+) is accepted;
- a partner's ortho/para note is written only for the ortho/para partners
  (legacy never reset its flag between partners).
LAMDA partner ids: 1 H2, 2 para-H2, 3 ortho-H2, 4 e-, 5 H, 6 He, 7 H+.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "kromex"))  # constants

import constants  # noqa: E402

HC_K = constants.hc_k  # K cm: second radiation constant hc/k
PARTNERS = {1: "H2", 2: "H2pa", 3: "H2or", 4: "e", 5: "H", 6: "He",
            7: "H+"}


def sections(lines):
    """the non-comment lines grouped under their '!' headers"""
    out, cur = [], None
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("!"):
            cur = [s, []]
            out.append(cur)
        elif cur is not None:
            cur[1].append(s)
    return out


def f90(x):
    x = x.lower()
    if "e" in x:
        return x.replace("e", "d")
    return x if "d" in x else x + "d0"


def convert(src, name=None):
    sec = sections(open(src).read().splitlines())
    it = iter(sec)

    def nxt(*need, without=()):
        """the next section whose header has every word of need and none
        of without (LAMDA files word their headers differently)"""
        for h, body in it:
            H = h.upper()
            if all(k in H for k in need) and not any(k in H for k in
                                                     without):
                return body
        raise ValueError("%s: no section with %s" % (src, need))

    mol = name or nxt("MOLECULE")[0].split()[0]
    nxt("WEIGHT")
    nlev = int(nxt("NUMBER", "LEVELS")[0])
    levels = [ln.split() for ln in nxt("LEVEL", without=("NUMBER",))][:nlev]
    nrad = int(nxt("NUMBER", "RADIATIVE")[0])
    rads = [ln.split() for ln in nxt("TRANS", without=("NUMBER",))][:nrad]
    npart = int(nxt("NUMBER", "PARTNER")[0])
    out = ["# converted from LAMDA %s by tools/lamda2krome.py" % src,
           "metal:%s" % mol,
           "# level n: energy (K), degeneracy g"]
    for lv in levels:
        out.append("level %d: %r, %s" % (int(lv[0]) - 1,
                                        float(lv[1]) * HC_K, lv[2]))
    out.append("# up -> low, Aij (1/s)")
    for r in rads:
        out.append("%d -> %d, %s" % (int(r[1]) - 1, int(r[2]) - 1, r[3]))
    for _ in range(npart):
        pid = int(nxt("BETWEEN")[0].split()[0])
        if pid not in PARTNERS:
            raise ValueError("%s: unknown LAMDA partner id %d" % (src, pid))
        p = PARTNERS[pid]
        ncoll = int(nxt("NUMBER", "TRANS")[0])
        nxt("NUMBER", "TEMP")
        temps = nxt("TEMP", without=("NUMBER",))[0].split()
        rows = nxt("TRANS", without=("NUMBER",))[:ncoll]
        out.append("# collider, level_up, level_down, de-excitation rate "
                   "(cm^3/s), LAMDA partner %d" % pid)
        if pid in (2, 3):
            out.append("ortho/para: 3/1")
        xv = "(/" + ", ".join(f90(t) for t in temps) + "/)"
        for r in rows:
            f = r.split()
            yv = "(/" + ", ".join(f90(v) for v in f[3:]) + "/)"
            out.append("%s, %d, %d, fspline(%s, %s, Tgas)"
                       % (p, int(f[1]) - 1, int(f[2]) - 1, xv, yv))
    out.append("end metal")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--name", help="the coolant species (default: the "
                    "LAMDA molecule name)")
    a = ap.parse_args()
    if a.input == a.output:
        sys.exit("input and output are the same file")
    txt = convert(a.input, a.name)
    open(a.output, "w").write(txt)
    print("wrote %s" % a.output)


if __name__ == "__main__":
    main()
