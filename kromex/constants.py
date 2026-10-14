"""The physical constants of KROME-X, in one place: CODATA 2018 (the
exact SI values of the 2019 redefinition where they apply) and the IAU
2012/2015 astronomical units, in cgs unless the name says otherwise.

The codegen emits them as krome::constants in every generated network
(the C++ code, the runtime headers and the network expressions use those
names); runtime/interfaces/C/krome.h (KROME_CONST_*), runtime/interfaces/python/krome.py and
runtime/interfaces/fortran/krome.f90 carry copies that tests/test_constants.py
checks against this table. Legacy KROME used older values (CODATA 2006-
2010, pi = 3.14159265359, 1 Mpc = 3.08e24 cm, a 365-day year).
"""

import math

_h = 6.62607015e-34      # J s, exact
_kB = 1.380649e-23       # J / K, exact
_e = 1.602176634e-19     # C, exact
_NA = 6.02214076e23      # 1 / mol, exact
_c = 2.99792458e8        # m / s, exact
_au = 1.495978707e13     # cm, IAU 2012 (exact)
_Ry = 13.605693122994    # eV, CODATA 2018

# name -> (value, unit, note); the order is the order of the listings
CONSTANTS = {
    "boltzmann_eV": (_kB / _e, "eV / K", "exact"),
    "boltzmann_J": (_kB, "J / K", "exact"),
    "boltzmann_erg": (_kB * 1e7, "erg / K", "exact"),
    "planck_eV": (_h / _e, "eV s", "exact"),
    "planck_J": (_h, "J s", "exact"),
    "planck_erg": (_h * 1e7, "erg s", "exact"),
    "hbar_erg": (_h * 1e7 / (2.0 * math.pi), "erg s", "exact"),
    "hc_k": (_h * _c * 1e2 / _kB, "cm K",
             "second radiation constant: E/k [K] = hc_k x E [cm^-1]"),
    "gravity": (6.67430e-8, "cm^3 g^-1 s^-2", "CODATA 2018"),
    "e_mass": (9.1093837015e-28, "g", "CODATA 2018"),
    "p_mass": (1.67262192369e-24, "g", "CODATA 2018"),
    "n_mass": (1.67492749804e-24, "g", "CODATA 2018"),
    "amu": (1.66053906660e-24, "g", "CODATA 2018"),
    "clight": (_c * 1e2, "cm s^-1", "exact"),
    "pi": (math.pi, "", ""),
    "eV_to_erg": (_e * 1e7, "eV -> erg", "exact"),
    "ry_to_eV": (_Ry, "Rydberg -> eV", "CODATA 2018"),
    "ry_to_erg": (_Ry * _e * 1e7, "Rydberg -> erg", "CODATA 2018"),
    "seconds_per_year": (3.15576e7, "yr -> s", "Julian year, 365.25 d"),
    "km_to_cm": (1e5, "km -> cm", ""),
    "pc_to_cm": (_au * 648000.0 / math.pi, "pc -> cm", "IAU 2015"),
    "Mpc_to_cm": (_au * 648000.0 / math.pi * 1e6, "Mpc -> cm", "IAU 2015"),
    "stefboltz_erg": (2.0 * math.pi ** 5 * (_kB * 1e7) ** 4
                      / (15.0 * (_h * 1e7) ** 3 * (_c * 1e2) ** 2),
                      "erg s^-1 cm^-2 K^-4", "exact"),
    "N_avogadro": (_NA, "mol^-1", "exact"),
    "Rgas_J": (_NA * _kB, "J K^-1 mol^-1", "exact"),
    "Rgas_kJ": (_NA * _kB * 1e-3, "kJ K^-1 mol^-1", "exact"),
}

globals().update({k: v[0] for k, v in CONSTANTS.items()})


def literal(x):
    """the shortest decimal that reads back as x"""
    return repr(float(x))


def cpp_block():
    """the constants as C++ (inside namespace krome)"""
    out = ["// physical constants (kromex/constants.py: CODATA 2018, IAU),",
           "// e.g. krome::constants::p_mass",
           "namespace constants {"]
    for k, (v, unit, note) in CONSTANTS.items():
        c = "; ".join(x for x in (unit, note) if x)
        out.append("constexpr double %s = %s;%s"
                   % (k, literal(v), ("  // " + c) if c else ""))
    out.append("}  // namespace constants")
    return out


def _note(unit, note):
    return "; ".join(x for x in (unit, note) if x)


def c_block():
    """the KROME_CONST_* macros of runtime/interfaces/C/krome.h"""
    out = ["/* physical constants (kromex/constants.py: CODATA 2018, IAU), the",
           " * numbers the solver uses (C++ krome::constants, Fortran",
           " * krome_constants) */"]
    for k, (v, unit, note) in CONSTANTS.items():
        c = _note(unit, note)
        out.append("#define KROME_CONST_%s %s%s"
                   % (k, literal(v), ("  /* %s */" % c) if c else ""))
    return out


def py_block():
    """the class attributes of runtime/interfaces/python/krome.py _Constants"""
    out = []
    for k, (v, unit, note) in CONSTANTS.items():
        c = _note(unit, note)
        out.append("    %s = %s%s" % (k, literal(v),
                                      ("  # " + c) if c else ""))
    return out


def f90_literal(x):
    s = literal(x)
    return s.replace("e", "d") if "e" in s else s + "d0"


def f90_block():
    """the parameters of runtime/interfaces/fortran/krome.f90 krome_constants"""
    out = []
    for k, (v, unit, note) in CONSTANTS.items():
        c = _note(unit, note)
        out.append("  real(8), parameter :: %s = %s%s"
                   % (k, f90_literal(v), ("  ! " + c) if c else ""))
    return out
