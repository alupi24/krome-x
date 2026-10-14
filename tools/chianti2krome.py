#!/usr/bin/env python3
"""CHIANTI atomic database -> KROME coolant file (the coolZ.dat format read
by codegen thermal.coolant_files), a port of legacy tools/chianti.py to Python 3 and
to the current CHIANTI formats (v8 and later).

    python3 tools/chianti2krome.py $XUVTOP coolChianti.dat c o ne fe
    python3 tools/chianti2krome.py $XUVTOP coolC.dat c_2 c_3 --max-level 10

The ions are CHIANTI names (c_2 = C+) or elements (c = every ion of C in
the database). Per ion:
- levels from .elvlc: the observed energy, the theoretical one when there
  is no observed value [cm^-1 -> K with hc/k = 1.4387769 K cm], g = 2J + 1;
- radiative transitions from .wgfa: up -> low, Aij [s^-1], wavelength
  [Angstrom] (CHIANTI's theoretical wavelengths are negative, written
  as their absolute value; 0, two-photon and autoionization lines, kept 0);
  repeated pairs are summed;
- electron de-excitation rates from .scups (Burgess & Tully 1992 scaled
  effective collision strengths, the six transition types, CHIANTI's
  natural cubic spline in the scaled temperature, descale_scups.pro):
  q_ul = 8.62913e-6 Upsilon / (sqrt(T) g_u) [cm^3 s^-1];
- proton de-excitation rates from .psplups: 5- or 9-point spline fits on
  uniform knots of the scaled EXCITATION rate coefficients q_lu
  [cm^3 s^-1] (the files' own units line), so by detailed balance
  q_ul = q_lu (g_l / g_u) exp(dE / kT). The fits hold over a limited
  range, which the files state in a comment ('valid for temperatures A
  to B K', read when present): outside it q_ul is held at its value at
  the nearer end. Without a stated range the lower end is dE / k (the
  Boltzmann factor at most e) and there is no upper end;
- a transition listed twice in a collision file is taken once; a pair
  listed upper level first (by the level energies) is swapped for the
  electrons (Upsilon is symmetric), and for the protons its rate is taken
  as the de-excitation one.
The rates are tabulated as log10 q on a log10 T grid (--tmin, --tmax,
--per-decade) and written as 1d1**fspline(logT, logq, logTgas).

Differences from legacy tools/chianti.py:
- CHIANTI v8+ formats: .elvlc energies in cm^-1 (legacy read Rydberg
  columns of the v7 files), .scups (legacy read the v7 .splups);
- the proton rates: legacy converted the .psplups rate coefficients as
  collision strengths (x 8.629e-6 / (sqrt(T) g_u), about 1e-7 of the
  right value at 1e4 K) and without detailed balance, so its H+
  collisions did nothing (ChiantiPy 0.16.1's Ion.upsilonDescale does
  the same, against the files' units);
- the table: 4 points per decade over 10 K - 1e9 K by default (legacy:
  10 points over 1e3 - 1e8 K), with CHIANTI's spline (legacy: scipy's
  interp1d on the knots). The natural spline is that of CHIANTI's IDL
  descale_scups; ChiantiPy uses scipy splrep (not-a-knot): over each
  ion's formation range (0.1 - 10 x the peak-abundance temperature) the
  two agree to a median 1e-5 - 1e-3 per collision strength, and differ
  more only between sparse, strongly curved knots, far from it;
- the wavelengths of theoretical lines kept (legacy wrote them negative),
  levels above --max-level dropped with their lines (legacy 20, the same
  default), the ion list and the database path are arguments.
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "kromex"))  # constants
import constants  # noqa: E402

RY_K = constants.ry_to_eV / constants.boltzmann_eV  # K per Rydberg
CM_K = constants.hc_k  # K per cm^-1 (hc / kB)
# h^2 / ((2 pi m_e)^1.5 kB^0.5) [cm^3 s^-1 K^0.5]
QCONST = constants.planck_erg ** 2 / (
    (2.0 * constants.pi * constants.e_mass) ** 1.5
    * constants.boltzmann_erg ** 0.5)
ELEMENTS = ("h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v "
            "cr mn fe co ni cu zn").split()


def natural_spline(xs, ys):
    """second derivatives of the natural cubic spline through (xs, ys)
    (IDL spl_init with zero end derivatives, as CHIANTI)"""
    n = len(xs)
    y2, u = [0.0] * n, [0.0] * n
    for i in range(1, n - 1):
        sig = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
                       - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]))
                / (xs[i + 1] - xs[i - 1]) - sig * u[i - 1]) / p
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


def spline_eval(xs, ys, y2, x):
    """the spline at x (the end intervals' cubics outside the knots, as
    IDL spl_interp)"""
    lo, hi = 0, len(xs) - 1
    if x <= xs[0]:
        hi = 1
    elif x >= xs[-1]:
        lo = hi - 1
    else:
        while hi - lo > 1:
            k = (hi + lo) // 2
            if xs[k] > x:
                hi = k
            else:
                lo = k
    h = xs[hi] - xs[lo]
    a = (xs[hi] - x) / h
    b = (x - xs[lo]) / h
    return (a * ys[lo] + b * ys[hi]
            + ((a ** 3 - a) * y2[lo] + (b ** 3 - b) * y2[hi]) * h * h / 6.0)


def descale(ttype, C, dE_ry, T, xs, ys, y2):
    """Burgess & Tully (1992) descaling at T [K] of a transition of energy
    dE_ry [Ry], type ttype and scaling parameter C"""
    kte = T / (dE_ry * RY_K)
    if ttype in (1, 4):
        xt = 1.0 - math.log(C) / math.log(kte + C)
    elif ttype in (2, 3, 5, 6):
        xt = kte / (kte + C)
    else:
        raise ValueError("unknown Burgess-Tully type %d" % ttype)
    s = spline_eval(xs, ys, y2, xt)
    if ttype == 1:
        return s * math.log(kte + math.e)
    if ttype == 2:
        return s
    if ttype == 3:
        return s / (kte + 1.0)
    if ttype == 4:
        return s * math.log(kte + C)
    if ttype == 5:
        return s / kte
    return 10.0 ** s


def data_lines(path):
    """the lines of a CHIANTI file up to its ' -1' terminator"""
    out = []
    for ln in Path(path).read_text(errors="replace").splitlines():
        if ln.strip() in ("-1", "-1.0") or ln.strip().startswith("%"):
            break
        if ln.strip():
            out.append(ln)
    return out


def read_elvlc(path):
    """{level: (E [K], g)} of a v8+ .elvlc (i7, a30, a5, i5, a5, f5.1,
    2f15.3: index, configuration, label, 2S+1, L, J, E_obs, E_theory
    [cm^-1]; E_obs -1 when unknown; some files add columns after these)"""
    lev = {}
    for ln in data_lines(path):
        i, J = int(ln[0:7]), float(ln[52:57])
        eobs, eth = float(ln[57:72]), float(ln[72:87])
        e = eobs if (eobs > 0.0 or (eobs == 0.0 and i == 1)) else eth
        if e < 0.0:
            raise ValueError("%s: level %d has no energy" % (path, i))
        lev[i] = (e * CM_K, 2.0 * J + 1.0)
    return lev


def read_wgfa(path):
    """[(low, up, wavelength [A], A [s^-1])] of a .wgfa (2i5, f15.3,
    2e15.3, labels)"""
    out = []
    for ln in data_lines(path):
        lo, up = int(ln[0:5]), int(ln[5:10])
        wvl, A = float(ln[10:25]), float(ln[40:55])
        out.append((lo, up, wvl, A))
    return out


def read_scups(path):
    """[(low, up, dE [Ry], type, C, xs, ys)] of a v8+ .scups (per
    transition: 'l1 l2 de gf lim nt type C', the nt scaled temperatures,
    the nt scaled collision strengths)"""
    rows = data_lines(path)
    out = []
    for k in range(0, len(rows) - 2, 3):
        h = rows[k].split()
        nt = int(h[5])
        xs = [float(v) for v in rows[k + 1].split()]
        ys = [float(v) for v in rows[k + 2].split()]
        if len(xs) != nt or len(ys) != nt:
            raise ValueError("%s: transition %s %s: %d points expected"
                             % (path, h[0], h[1], nt))
        out.append((int(h[0]), int(h[1]), float(h[2]), int(h[6]),
                    float(h[7]), xs, ys))
    return out


def psplups_range(path):
    """the validity range (Tlo, Thi) [K] a .psplups comment states ('Fits
    valid for temperatures 2e3 to 4e5 K'), or None"""
    num = r"([0-9.]+(?:[eE][+-]?[0-9]+)?)"
    for ln in Path(path).read_text(errors="replace").splitlines():
        m = re.search(r"valid[a-z ]*?" + num + r"\s*(?:K\s*)?to\s*" + num,
                      ln, re.I)
        if m:
            lo, hi = float(m.group(1)), float(m.group(2))
            if 0.0 < lo < hi and hi > 100.0:
                return lo, hi
    return None


def read_psplups(path):
    """[(low, up, dE [Ry], type, C, xs, ys)] of a .psplups (3i3, then
    e10.3 fields: gf, dE, C and the 5 or 9 spline values on uniform
    knots in [0, 1])"""
    out = []
    for ln in data_lines(path):
        lo, up, ty = int(ln[0:3]), int(ln[3:6]), int(ln[6:9])
        body = ln[9:].rstrip()
        f = [float(body[j:j + 10]) for j in range(0, len(body), 10)
             if body[j:j + 10].strip()]
        dE, C, ys = f[1], f[2], f[3:]
        if len(ys) not in (5, 9):
            raise ValueError("%s: %d %d: %d spline points" % (path, lo, up,
                                                               len(ys)))
        xs = [j / (len(ys) - 1.0) for j in range(len(ys))]
        out.append((lo, up, dE, ty, C, xs, ys))
    return out


def f90(x):
    return ("%.6e" % x).replace("e", "d")


def f90log(x):
    """a log10 value in the tables: 6 decimals"""
    return "%.6fd0" % x


def ion_name(ion):
    """CHIANTI c_2 -> the KROME species C+"""
    el, stage = ion.split("_")
    sym = el.capitalize()
    return sym + "+" * (int(stage) - 1)


def convert_ion(base, maxlev, temps):
    """the coolant block of one ion (base: the path without extension);
    None when the ion has no levels or no radiative data"""
    ion = base.name
    elv, wgf = base.with_suffix(".elvlc"), base.with_suffix(".wgfa")
    if not elv.exists() or not wgf.exists():
        return None, "no .elvlc / .wgfa"
    lev = {i: v for i, v in read_elvlc(elv).items() if i <= maxlev}
    out = ["", "#############", "# CHIANTI %s" % ion,
           "metal:%s" % ion_name(ion),
           "# level n: energy (K), degeneracy g"]
    for i in sorted(lev):
        out.append("level %d: %s, %g" % (i - 1, f90(lev[i][0]), lev[i][1]))
    rad = {}
    for lo, up, wvl, A in read_wgfa(wgf):
        if lo in lev and up in lev and A > 0.0:
            lo, up = (lo, up) if lev[up][0] >= lev[lo][0] else (up, lo)
            a0, w0 = rad.get((up, lo), (0.0, abs(wvl)))
            rad[(up, lo)] = (a0 + A, w0)
    out.append("")
    out.append("# up -> low, Aij (1/s), wavelength (angstrom)")
    for (up, lo) in sorted(rad, key=lambda k: (k[1], k[0])):
        A, wvl = rad[(up, lo)]
        out.append("%d -> %d, %s, %s" % (up - 1, lo - 1, f90(A), f90(wvl)))
    logT = ", ".join(f90log(math.log10(T)) for T in temps)
    ncoll = ndup = 0
    for coll, suffix in (("e", ".scups"), ("H+", ".psplups")):
        path = base.with_suffix(suffix)
        if not path.exists():
            continue
        data = read_scups(path) if coll == "e" else read_psplups(path)
        rng = psplups_range(path) if coll == "H+" else None
        out.append("")
        out.append("# collider, level_up, level_down, de-excitation rate "
                   "(cm^3/s), %s%s" % (path.name, "" if coll == "e" else
                                       ", fits valid %g - %g K" % rng
                                       if rng else ", no stated range"))
        seen = set()
        for lo, up, dE, ty, C, xs, ys in data:
            if lo not in lev or up not in lev or dE <= 0.0:
                continue
            if (lo, up) in seen:
                ndup += 1
                continue
            seen.add((lo, up))
            inverted = lev[up][0] < lev[lo][0]
            if inverted:
                lo, up = up, lo
            gl, gu = lev[lo][1], lev[up][1]
            y2 = natural_spline(xs, ys)
            lq = []
            for T in temps:
                if coll == "e":
                    u = max(descale(ty, C, dE, T, xs, ys, y2), 0.0)
                    q = QCONST * u / (math.sqrt(T) * gu)
                    lq.append(math.log10(q) if q > 0.0 else -40.0)
                    continue
                # protons: the excitation rate at T clamped to the fits'
                # range, then detailed balance with the level energies
                dEK = lev[up][0] - lev[lo][0]
                Tlo, Thi = rng if rng else (dEK, math.inf)
                Tc = min(max(T, Tlo), Thi)
                u = descale(ty, C, dE, Tc, xs, ys, y2)
                if inverted:  # the file's rate is already downward
                    lq.append(math.log10(u) if u > 0.0 else -40.0)
                    continue
                lq.append(math.log10(u * gl / gu) + dEK / Tc / math.log(10.0)
                          if u > 0.0 else -40.0)
            lq = [min(max(v, -40.0), 0.0) for v in lq]
            out.append("%s, %d, %d, 1d1**fspline((/%s/), (/%s/), logTgas)"
                       % (coll, up - 1, lo - 1, logT,
                          ", ".join(f90log(v) for v in lq)))
            ncoll += 1
    out.append("")
    out.append("endmetal")
    return "\n".join(out) + "\n", "%d levels, %d lines, %d rates%s" % (
        len(lev), len(rad), ncoll,
        ", %d repeated collision entries skipped" % ndup if ndup else "")


def ion_dirs(root, names):
    """the ion directories of the CHIANTI database root for the names
    (ions c_2, or elements c)"""
    out = []
    for nm in names:
        nm = nm.lower()
        if "_" in nm:
            el = nm.split("_")[0]
            d = root / el / nm
            if not d.is_dir():
                raise ValueError("%s: no ion %s" % (root, nm))
            out.append(d)
        else:
            if nm not in ELEMENTS:
                raise ValueError("unknown element %r" % nm)
            ions = [d for d in (root / nm).glob(nm + "_*")
                    if d.is_dir() and d.name.split("_")[1].isdigit()]
            out += sorted(ions, key=lambda d: int(d.name.split("_")[1]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("database", type=Path,
                    help="the CHIANTI database root (XUVTOP: the "
                    "directory with the element directories)")
    ap.add_argument("output", type=Path)
    ap.add_argument("ions", nargs="*", default=["c", "o", "ne", "fe"],
                    help="CHIANTI ions (c_2) or elements (c: all its ions)"
                    "; default c o ne fe, as legacy")
    ap.add_argument("--max-level", type=int, default=20,
                    help="keep the CHIANTI levels 1..N (default 20)")
    ap.add_argument("--tmin", type=float, default=10.0)
    ap.add_argument("--tmax", type=float, default=1e9)
    ap.add_argument("--per-decade", type=int, default=4,
                    help="rate-table points per decade of T (default 4: "
                    "the chianti test's cooling curve within 1e-3 of 8)")
    a = ap.parse_args()
    if not (0.0 < a.tmin < a.tmax) or a.per_decade < 1 or a.max_level < 2:
        sys.exit("need 0 < tmin < tmax, per-decade >= 1, max-level >= 2")
    n = int(round(math.log10(a.tmax / a.tmin) * a.per_decade)) + 1
    temps = [a.tmin * (a.tmax / a.tmin) ** (k / (n - 1.0)) for k in range(n)]
    ver = a.database / "VERSION"
    head = ["# converted from CHIANTI%s by tools/chianti2krome.py" % (
        " " + ver.read_text().split()[0] if ver.exists() else ""),
        "# rates: log10 q on %d log10 T points, %g - %g K" % (n, a.tmin,
                                                             a.tmax)]
    blocks = []
    try:
        dirs = ion_dirs(a.database, a.ions)
        for d in dirs:
            txt, info = convert_ion(d / d.name, a.max_level, temps)
            print("%-8s %-10s %s" % (d.name, ion_name(d.name),
                                     info if txt else "skipped: " + info))
            if txt:
                blocks.append(txt)
    except (ValueError, OSError) as exc:
        sys.exit("ERROR: %s" % exc)
    if not blocks:
        sys.exit("ERROR: no ion converted")
    a.output.write_text("\n".join(head) + "\n" + "".join(blocks))
    print("wrote %s (%d coolants): codegen thermal.coolant_files=%s" % (
        a.output, len(blocks), a.output))


if __name__ == "__main__":
    main()
