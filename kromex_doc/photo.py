"""The photo reactions: their cross sections sigma(E), from the sources the
code generator reads (the expression of photo(xsec=...), the
photoionization database for xsec=auto, the Leiden and SWRI files), and
their unshielded rate in the Draine (1978) interstellar field."""

import math
import re

import numpy as np

import codegen
import paths

from . import evaluate
from .evaluate import sigma_v96, sigma_vy95

FITS = {"sigma_v96": sigma_v96, "sigma_vy95": sigma_vy95}

DB = paths.ROOT / "data" / "database"
DRAINE = (5.0, 13.6)   # eV


def draine_F(E):
    """Draine 1978 photon intensity [cm^-2 s^-1 sr^-1 eV^-1], 5-13.6 eV
    (runtime/core/spectra.hpp)"""
    F = 1.658e6 * E - 2.152e5 * E ** 2 + 6.919e3 * E ** 3
    return np.where((E >= DRAINE[0]) & (E <= DRAINE[1]), F, 0.0)


class Xsec:
    """sigma(E) of one photo reaction: kind, source (text), the energy
    window [Emin, Emax], and sigma on an array of energies (zero outside
    the window); ok False with a reason when the data are missing"""

    def __init__(self, rx, net, photo_db, warnings):
        r = rx["rate"]
        x = r["xsec"]
        self.rx, self.net = rx, net
        self.kind = x["kind"]
        self.ok, self.reason = True, ""
        self.Emin = r.get("Emin", 0.0)
        self.Emax = r.get("Emax", math.inf)
        self.table = None
        self.segs = None
        self.shield = r.get("shield")
        if self.kind == "expr":
            self.source = "the expression in the network"
            try:     # a function only the generated code computes
                self(np.array([13.6]))
            except evaluate.NotEvaluable as e:
                self.ok, self.reason = False, str(e)
        elif self.kind == "h2pd":
            self.ok = False
            self.source = "the H2 line data (data/H2pdB.dat)"
            self.reason = ("H2 photodissociation from its lines: computed "
                           "on the photon bins by the generated code")
        elif self.kind == "auto":
            reacts = sorted(q["sp"] for q in rx["reactants"])
            prods = sorted(p["sp"] for p in rx["products"])
            segs = [e for e in photo_db if e.get("type") == "photoion"
                    and e.get("reacts") == reacts
                    and e.get("prods") == prods]
            segs.sort(key=lambda e: e["limits"][0])
            self.source = ("Verner et al. 1996 fits "
                           "(data/database/photoionization.dat)")
            if not segs:
                self.ok = False
                self.reason = "no entry in the photoionization database"
            else:
                self.segs = []
                for e in segs:
                    # sigma_v96(...), or a sum of sigma_vy95(...) (one per
                    # shell, above the inner-shell edge)
                    rate = codegen.d2e(e["rate"]).strip()
                    terms = re.findall(r"(sigma_v96|sigma_vy95)\(\s*"
                                       r"energy_ev\s*,([^)]*)\)", rate, re.I)
                    self.segs.append((e["limits"][0], [
                        (FITS[f.lower()], [float(v) for v in a.split(",")])
                        for f, a in terms]))
                if "Emin" not in r:
                    self.Emin = segs[0]["limits"][0]
        else:
            ref = (x.get("ref") or "").lower()
            self.kind = ref
            self.source = {"leiden": "the Leiden database",
                           "swri": "the SWRI (phidrates) database"}.get(
                               ref, "%r (not implemented)" % x.get("ref"))
            try:
                if ref == "leiden":
                    self.table = codegen.leiden_table(DB, rx)
                elif ref == "swri":
                    self.table = codegen.swri_table(DB, rx, warnings)
                else:
                    self.ok = False
                    self.reason = "source %r not implemented" % x.get("ref")
            except FileNotFoundError as exc:
                self.ok, self.reason = False, str(exc)
            if self.table:
                if "Emin" not in r:
                    self.Emin = self.table[0][0]
                if "Emax" not in r:
                    self.Emax = self.table[0][-1]
            elif self.ok:
                self.ok, self.reason = False, "an empty table"

    def __call__(self, E):
        E = np.asarray(E, dtype=float)
        if not self.ok:
            raise evaluate.NotEvaluable(self.reason)
        if self.kind == "expr":
            ev = evaluate.Evaluator(self.net, {"E_eV": E})
            with np.errstate(all="ignore"):
                s = np.broadcast_to(ev.ev(self.rx["rate"]["xsec"]["ast"]),
                                    E.shape).astype(float)
        elif self.segs is not None:
            s = np.zeros_like(E)
            for lo, terms in self.segs:    # each from its own threshold
                with np.errstate(all="ignore"):
                    v = sum(fit(E, *args) for fit, args in terms)
                s = np.where(E >= lo, v, s)
        else:
            tE, tS = self.table
            s = np.interp(E, tE, tS, left=0.0, right=0.0)
        inside = (E >= self.Emin) & (E <= self.Emax)
        return np.where(inside & np.isfinite(s), s, 0.0)

    def grid(self, lo, hi, n=2000):
        """energies for plots and integrals: log-spaced, with the table's
        points and the window edges inside [lo, hi]"""
        E = np.geomspace(lo, hi, n)
        extra = [self.Emin, self.Emax]
        if self.table:
            extra += list(self.table[0])
        if self.segs:
            extra += [s[0] for s in self.segs]
        extra = [e for e in extra if lo <= e <= hi and math.isfinite(e)]
        for e in list(extra):     # the steps: a point on each side
            extra += [e * (1 - 1e-9), e * (1 + 1e-9)]
        E = np.unique(np.concatenate([E, extra]))
        return E[(E >= lo) & (E <= hi)]

    def draine_rate(self):
        """the unshielded rate in the Draine field [s^-1]: 4 pi int F
        sigma dE"""
        E = self.grid(*DRAINE, n=4000)
        return float(4.0 * math.pi * np.trapezoid(draine_F(E) * self(E), E)
                     if hasattr(np, "trapezoid") else
                     4.0 * math.pi * np.trapz(draine_F(E) * self(E), E))


def xsecs(net):
    """{reaction id: Xsec} for the photo reactions, and the generator's
    warnings on the data files"""
    db = codegen.load_photo_db(DB / "photoionization.dat")
    warnings = []
    out = {rx["id"]: Xsec(rx, net, db, warnings) for rx in net.reactions
           if rx["rate"]["kind"] == "photo"}
    return out, warnings
