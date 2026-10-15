"""Species data for the documentation: formation enthalpies from the
thermochemical tables KROME-X ships, as the code generator reads them
(data/thermo_burcat.toml: Burcat's database, Goos, Burcat & Ruscic, for
the species of the networks; data/thermoNIST.dat, NIST-JANAF, for a few
others)."""

import paths
import tomlio

T298 = 298.15
R_KJ = 8.314462618e-3          # kJ / (mol K)
KJMOL_TO_K = 1.0 / R_KJ        # an enthalpy [kJ/mol] as a temperature [K]


def load_nist(path):
    """NAME (upper case) -> the formation enthalpy at 298.15 K [kJ/mol],
    the records' last column (H in H(T) - H(298.15) = A t + ... + F - H)"""
    out = {}
    for row in path.read_text().splitlines():
        s = row.strip()
        if not s or s[0] == "#":
            continue
        a = s.split(":")
        out[a[0].upper()] = float(a[-1])
    return out


class Thermo:
    def __init__(self):
        d = paths.ROOT / "data"
        self.burcat = {k: v["Hf298"] for k, v in tomlio.load(
            d / "thermo_burcat.toml").get("species", {}).items()}
        self.nist = load_nist(d / "thermoNIST.dat")

    def h298(self, sp):
        """the standard formation enthalpy at 298.15 K [kJ/mol] (the
        electron: 0), or None without data"""
        if sp == "e-":
            return 0.0
        if sp in self.burcat:
            return self.burcat[sp]
        key = sp.upper()
        return self.nist.get(key)

    def dH(self, reactants, products):
        """the reaction enthalpy at 298.15 K [K] (negative: exothermic), or
        None"""
        hs = [self.h298(x) for x in reactants + products]
        if None in hs:
            return None
        nr = len(reactants)
        return (sum(hs[nr:]) - sum(hs[:nr])) * KJMOL_TO_K
