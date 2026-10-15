"""Shared by kida2krome and umist2krome: astrochemical rate
databases (KIDA, UMIST) as KROME-X networks (.knet).

Species: an isomer prefix becomes a label (l-C3H -> C3H_l, c-C3H2+ ->
C3H2_c+: a .knet name is a formula with _labels); the database's
pseudo-species (cosmic rays, photons) leave the reaction, its rate
carries them. Rates: written with the parameters crflux (the cosmic-ray
ionisation rate, in the database's convention), Av (visual extinction,
KROME-X's implicit parameter) and a UV scale (the field in units of the
database's, default 1).

A rate defined over several temperature ranges becomes a piecewise rate:
the ranges must tile, so a gap between two (UMIST: 10-100 K and 101-3000
K) is closed at the start of the next one, with a comment."""

import math
import re
from pathlib import Path

import kromex  # noqa: E402,F401 (puts kromex/ on the path)
import knet  # noqa: E402


def species(name):
    """a database name as a .knet species"""
    n = name.strip()
    if n in ("e-", "E-", "E"):
        return "e-"
    m = re.match(r"^([a-z])-(.+?)([+\-]*)$", n)
    if m:     # isomer prefix: a label
        return "%s_%s%s" % (m.group(2), m.group(1), m.group(3))
    return n


def num(x):
    """a number as written in a .knet (shortest round-trip form)"""
    v = float(x)
    if v == int(v) and abs(v) < 1e6:
        return "%d" % v if v != 0 else "0"
    return repr(v).replace("e+", "e")


def kooij(a, b, c):
    """alpha (T/300)^beta exp(-gamma/T), the factors that are not 1 left
    out"""
    a, b, c = float(a), float(b), float(c)
    if a == 0.0:
        return "0"
    s = num(a)
    if b != 0.0:
        s += "*(Tgas/300)**(%s)" % num(b)
    if c != 0.0:
        s += "*exp(-(%s)/Tgas)" % num(c)
    return s


def trange(lo, hi):
    """a .knet T range ("" for none)"""
    if lo is None or hi is None:
        return "[]"
    return "[%s, %s]" % (num(lo), "inf" if math.isinf(hi) else num(hi))


def single_t(lo, hi):
    """a rate measured at one temperature (Tl = Tu): written without a
    range (with KROME-X's default, rates.outside = "asym", the range
    [T, T] would make it zero below T); returns (range, note)"""
    if lo is not None and lo == hi:
        return "[]", "measured at %s K: no range" % num(lo)
    return None, None


def tile(pieces):
    """the pieces [(lo, hi, ...)] of one rate, sorted, with each gap
    closed at the start of the next piece; returns (pieces, notes)"""
    pieces = sorted(pieces, key=lambda p: p[0])
    out, notes = [], []
    for p in pieces:
        if out and out[-1][1] < p[0]:
            notes.append("gap %s-%s K closed: the piece below extends to %s"
                         % (num(out[-1][1]), num(p[0]), num(p[0])))
            out[-1] = (out[-1][0], p[0]) + tuple(out[-1][2:])
        out.append(tuple(p))
    return out, notes


class Writer:
    def __init__(self, source, database, params):
        self.source = source
        self.database = database
        self.params = params        # [(name, default, comment)]
        self.lines = []
        self.skipped = []
        self.suspicious = []
        self.nrx = 0

    def comment(self, text):
        for t in text.splitlines():
            self.lines.append("# " + t if t else "#")

    def reaction(self, reac, prod, rng, rate, ref=None, note=None):
        if note:
            self.comment(note)
        if ref:
            self.lines.append("ref %s" % ref)
        self.lines.append("%-40s %-18s %s" % (
            "%s -> %s" % (" + ".join(reac), " + ".join(prod)), rng, rate))
        self.nrx += 1

    def skip(self, raw, why):
        self.skipped.append(why)
        self.comment("SKIPPED (%s): %s" % (why, raw.strip()))

    def text(self, header):
        h = ["krome_format 1",
             "# converted from %s (%s) by kromex_convert" % (
                 Path(self.source).name, self.database)]
        h += ["# " + x for x in header]
        h.append("")
        for name, default, comment in self.params:
            h.append("param %-12s = %-10s # %s" % (name, num(default),
                                                    comment))
        h.append("")
        return "\n".join(h + self.lines) + "\n"

    def write(self, path, header):
        path = Path(path)
        path.write_text(self.text(header))
        net = knet.parse_network(path)
        print("wrote %s: %d reactions, %d species%s" % (
            path, len(net.reactions), len(net.species),
            "; %d entries skipped (comments in the file)" % len(
                self.skipped) if self.skipped else ""))
        opaque = sorted(s for s, d in net.species.items()
                        if d["composition"] is None)
        if opaque:
            print("species of unknown composition (no conservation check): "
                  + ", ".join(opaque))
        for e in net.errors:
            print("ERROR: %s" % e)
        return net
