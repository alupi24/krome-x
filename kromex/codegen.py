#!/usr/bin/env python3
"""Code generator: NetworkIR (from knet.py) -> C++ reference backend.

Emits one self-contained header `net_<name>.hpp` with:
  - species table (indices, names, masses, charges)
  - Params struct (runtime parameters with file defaults)
  - rates(Tgas, n, p, k): rate coefficients, piecewise groups clamped per
    network-format.md 3.1 (build option `outside=zero` supported via
    OUTSIDE_ZERO)
  - rhs(n, k, dn): chemical right-hand side
  - jacobian(Tgas, n, p, k, J): analytic dense Jacobian d(dn_i)/d(n_j),
    including d(k)/d(n) terms from density-dependent rate expressions
    (extern/table calls are treated as constant w.r.t. n)

Scope (Step 2 milestone 1): rates of kind 'expr'. photo/xray/surface/auto/
extern-based rates generate a stub k=0 plus a listing in `unsupported()`,
and the generator warns; those become runtime-module calls in later
milestones.

CLI: python3 codegen.py network.knet -o outdir
"""

import argparse
import copy
import math
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import knet  # noqa: E402
import constants  # noqa: E402

# ---------------------------------------------------------------------------
# AST helpers: constant folding constructors + differentiation
# ---------------------------------------------------------------------------

ZERO = {"n": 0.0}
ONE = {"n": 1.0}


def is_num(a, v=None):
    return "n" in a and (v is None or a["n"] == v)


def add(a, b):
    if is_num(a, 0.0):
        return b
    if is_num(b, 0.0):
        return a
    if is_num(a) and is_num(b):
        return {"n": a["n"] + b["n"]}
    return {"b": "+", "l": a, "r": b}


def sub(a, b):
    if is_num(b, 0.0):
        return a
    if is_num(a) and is_num(b):
        return {"n": a["n"] - b["n"]}
    if is_num(a, 0.0):
        return {"u": "-", "x": b}
    return {"b": "-", "l": a, "r": b}


def mul(a, b):
    if is_num(a, 0.0) or is_num(b, 0.0):
        return ZERO
    if is_num(a, 1.0):
        return b
    if is_num(b, 1.0):
        return a
    if is_num(a) and is_num(b):
        return {"n": a["n"] * b["n"]}
    return {"b": "*", "l": a, "r": b}


def div(a, b):
    if is_num(a, 0.0):
        return ZERO
    if is_num(b, 1.0):
        return a
    return {"b": "/", "l": a, "r": b}


def pw(a, b):
    if is_num(b, 1.0):
        return a
    if is_num(b, 0.0):
        return ONE
    return {"b": "**", "l": a, "r": b}


def call(f, *a):
    return {"f": f, "a": list(a), "kind": "builtin"}


import paths  # noqa: E402 (where the runtime and data are)

# Burcat's thermochemical database (data/thermo_burcat.toml) is not
# covered by KROME-X's GPL-3.0: its authors' terms, printed when a project
# uses it and written next to its coefficients in the generated code
BURCAT_TERMS = (
    "Burcat's thermochemical database (E. Goos, A. Burcat & B. Ruscic, "
    "https://respecth.elte.hu/burcat/; cite Burcat & Ruscic 2005), used "
    "here, is not covered by KROME-X's GPL-3.0 license: it is free to use with "
    "citation, but the database or parts of it may not be included in "
    "commercial software without the authors' permission "
    "(Elke.Goos@dlr.de); see data/data_sources.md")
BURCAT_COMMENT = [
    "// Thermochemical data from Burcat's database (E. Goos, A. Burcat &",
    "// B. Ruscic, https://respecth.elte.hu/burcat/; Burcat & Ruscic 2005).",
    "// Not covered by KROME-X's GPL-3.0: free to use with citation; the",
    "// database or parts of it may not be included in commercial software",
    "// without the authors' permission (Elke.Goos@dlr.de). See",
    "// data/data_sources.md in KROME-X."]

DEFAULT_DATA_DIR = paths.ROOT / "data" / "database"


def krome_version():
    """the KROME-X version (v2/VERSION), with the git revision of the KROME-X
    tree as build metadata when it is a checkout: 2.0.0-dev+g04bb53f
    (.dirty with uncommitted changes)"""
    root = paths.ROOT
    try:
        ver = (root / "VERSION").read_text().split()[0]
    except (OSError, IndexError):
        ver = "unknown"
    try:
        rev = subprocess.run(["git", "-C", str(root), "describe", "--always",
                              "--dirty=.dirty", "--abbrev=7",
                              "--exclude=*"], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True,
                             timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        rev = ""
    return ver + ("+g" + rev if rev else "")

# extension functions implemented by runtime modules (runtime/core/dust_rates.hpp, ...).
# Call sites get the implicit state arguments (Tgas, n, p) appended; anything
# not listed here stays a template stub returning 0.
KNOWN_EXTERNS = {
    "H_recombination_on_dust":
        "double H_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "He_recombination_on_dust":
        "double He_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "C_recombination_on_dust":
        "double C_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "Si_recombination_on_dust":
        "double Si_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "O_recombination_on_dust":
        "double O_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "Fe_recombination_on_dust":
        "double Fe_recombination_on_dust(double Tgas, const double* n, "
        "const Params& p);",
    "krate_stickSi":
        "double krate_stickSi(int sp, double Tgas, const double* n, "
        "const Params& p);",
    "krate_evaporation":
        "double krate_evaporation(int sp, double Tgas, const double* n, "
        "const Params& p);",
    "krate_nonthermal_evaporation":
        "double krate_nonthermal_evaporation(int sp, double Gnot, double Av, "
        "double crflux, double yield, double Tgas, const double* n, "
        "const Params& p);",
    "sqimass":
        "double sqimass(int sp, double Tgas, const double* n, "
        "const Params& p);",
    "fselfH2":
        "double fselfH2(double N, double b, double Tgas, const double* n, "
        "const Params& p);",
    "fselfH2_WG11":
        "double fselfH2_WG11(double N, double Tgas, const double* n, "
        "const Params& p);",
    "fselfH2_R14":
        "double fselfH2_R14(double N, double bturb, double Tgas, "
        "const double* n, const Params& p);",
    "H2_solomonLW":
        "double H2_solomonLW(double flux_eV, double Tgas, const double* n, "
        "const Params& p);",
    "dissH2_Martin96":
        "double dissH2_Martin96(double Tgas, const double* n, "
        "const Params& p);",  # runtime/core/special_rates.hpp (MSM96 H2 dissociation)
    # H2 formation on dust (runtime/core/dust_rates.hpp), for H + {H} -> H2
    "krate_H2dust":
        "double krate_H2dust(double Tgas, const double* n, const Params& p);",
    "krate_H2dust_Jura":
        "double krate_H2dust_Jura(double Tgas, const double* n, "
        "const Params& p);",
    "elec_recomb_ST93":
        "double elec_recomb_ST93(double nabund, double nelec, double ntot, "
        "double nucleiH, double Trad, double Tgas, const double* n, "
        "const Params& p);",  # runtime/core/special_rates.hpp (Sasaki & Takahara 93)
    # cosmological recombination as RECFAST 1.5.2 (special_rates.hpp)
    "recfast_H_rec":
        "double recfast_H_rec(double nH0, double Tgas, const double* n, "
        "const Params& p);",
    "recfast_H_ion":
        "double recfast_H_ion(double nH0, double Tgas, const double* n, "
        "const Params& p);",
    "recfast_He_rec":
        "double recfast_He_rec(double nHe0, double nH0, double Tgas, "
        "const double* n, const Params& p);",
    "recfast_He_ion":
        "double recfast_He_ion(double nHe0, double nH0, double Tgas, "
        "const double* n, const Params& p);",
}

# legacy extension functions that codegen DEFINES in the generated code
# (instead of the zero stub): name -> whether calls get the implicit
# (Tgas, n, p) appended. krome_fshield's body depends on --shielding.
GENERATED_EXTERNS = {
    "fHnOj": False,               # CO photodissociation attenuation(Av)
    "n_global": True,             # abundance at the start of the call
    "get_table_Tdust": True,      # alias of the live table_Tdust
    "krome_fshield": True,        # H2 self-shielding (legacy -shielding)
    "get_ratioFluxDraine": True,
    "krate_H2dust_table": True,   # H2 on dust from the dust tables  # photobin field / Draine ISRF
    "stars_screen": True,         # nuclear screening (--screening)
    # kinetic nucleation theory (legacy krome_subs, #KROME_use_cluster_*)
    "cluster_growth_rate": True,
    "general_cluster_growth_rate": True,
    "troe_falloff": False,        # simplified Troe falloff (N = 1, c = 0)
    "shield_dust": True,          # dust shielding exp(-sigma_d N_H)
    "num2col": True,              # legacy num2col(x): printed as num2col_
}

# generated externs whose value does not depend on the current state
# densities beyond their explicit arguments: "const" parameters and
# constants only, "call" constant during one solve() call (n_global: the
# input abundances of the call), "field" the photon field (constant unless
# --evolve-rad). With n-independent arguments their rates need no
# numerical dk/dn, and they may enter rates_fixed.
EXTERN_STATE = {
    "fHnOj": "const",
    "troe_falloff": "const",
    "n_global": "call",
    "get_ratioFluxDraine": "field",
}

# device-callable kernels: every generated function that the rhs /
# Jacobian path calls is KROME_HD (__host__ __device__ under a CUDA / HIP
# compiler, which Kokkos' device back-ends use; empty on a host build). A
# back-end may predefine it (e.g. as KOKKOS_INLINE_FUNCTION's attributes).
KROME_HD_DEF = [
    "#ifndef KROME_HD",
    "#if defined(__CUDACC__) || defined(__HIPCC__)",
    "#define KROME_HD __host__ __device__",
    "#else",
    "#define KROME_HD",
    "#endif",
    "#endif",
    "",
]

# monomer radii [cm] for cluster growth (legacy kromelib get_radius)
NUCLEATION_RADII = {
    "TIO2": 1.62e-8,      # Jeong et al. 2000, doi:10.1088/0953-4075/33/17/319
    "AL2O3": 3.304e-8,    # O-Al-O distance, Archibong et al. 1999
    "MGO": 0.865e-8,      # half a bond length, Farrow et al. 2014
    "SIO": 0.75765e-8,    # bond length, Bromley et al. 2016
}

# atomic numbers (nuclear screening, --screening)
ATOMIC_Z = {e: i + 1 for i, e in enumerate(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn".split())}

# extension functions with NO implicit state arguments (usable inside
# photo cross-section expressions, which only see E_eV)
PLAIN_EXTERNS = {
    "H2_sigmaLW": "double H2_sigmaLW(double E_eV);",
    "H2_sigmaLW_R14": "double H2_sigmaLW_R14(double E_eV);",
}

# universal runtime parameters: a network gets only those its generated
# physics or its own expressions use (plus ALWAYS_PARAMS); a default of
# None means REQUIRED: no default exists that would not misrepresent the
# process (NaN in the generated Params, rejected by the solver)
UNIVERSAL_PARAMS = [
    ("Tdust", 10.0, "dust temperature [K]"),
    ("G0", 0.0, "local radiation field [Habing]"),
    ("d2g", 1.0e-2, "dust-to-gas mass ratio (0.01 = Milky Way)"),
    ("Av", 0.0, "visual extinction [mag]"),
    ("gamma_ad", 5.0 / 3.0, "adiabatic index for the thermal ODE"),
    ("Tcmb", 2.725, "CMB temperature [K]; also the rates' lower temperature "
     "cap: rates are evaluated at min(max(Tgas, Tcmb), 1e9 K), the thermal "
     "terms at min(max(Tgas, 2.73 K), 1e9 K) (the legacy caps)"),
    ("tff", None, "compression timescale [s]: heating n kB T / tff"),
    ("zredshift", 0.0, "redshift"),
    ("hubble0", None, "Hubble constant H0 [s^-1]"),
    ("omega0", None, "matter density parameter"),
    ("visc_nu", None, "kinematic viscosity [cm^2 s^-1]"),
    ("visc_omega", None, "orbital frequency [s^-1] of the sheared flow"),
    ("opratio", 3.0, "H2 ortho:para ratio (metal line cooling colliders)"),
    ("total_Z", 1.0, "metallicity [Z_sun] scaling the ZCIE cooling/heating "
     "tables"),
    ("clump", 1.0, "clumping factor of the Jura H2-on-dust rate"),
    ("Tfloor", 2.725, "cooling-floor temperature [K]: the radiative "
     "coolers are evaluated minus their value at Tfloor"),
    ("vturb", 0.0, "turbulent velocity [cm s^-1] in the line widths of "
     "thermal.escape"),
    ("Lcol", None, "the column-density length [cm] of radiation.num2col = 'length': "
     "N = n Lcol (e.g. the cell size, or a Sobolev or shielding length)"),
]
ALWAYS_PARAMS = ("Tcmb",)  # gamma_ad: scanned (only with --gamma=CONST)

# diatomic constants for --gamma (legacy kromelib get_ve_vib / get_be_rot,
# NIST / Atkins, cm^-1 -> K with hc/k; legacy multiplied by 1.42879, a
# typo for 1.43878, so its temperatures were 0.7 % low)
GAMMA_VIB_K = {s: v * constants.hc_k for s, v in {
    "H2": 4401.213, "HD": 3813.15, "D2": 3115.5, "C2": 1855.0663,
    "C2-": 1781.189, "CH": 2860.7508, "CO": 2169.75589, "CO+": 2214.127,
    "N2": 2358.57, "N2+": 2207.0115, "NH": 3282.72, "NO": 1904.1346,
    "NO+": 2376.72, "O2": 1580.161, "O2+": 1905.892, "OH": 3737.761}.items()}
GAMMA_ROT_K = {s: v * constants.hc_k for s, v in {
    "H2": 60.853, "H2+": 42.9, "HD": 45.644, "D2": 30.443, "N2": 1.9982,
    "O2": 1.4264, "CO": 2.78}.items()}
GAMMA_MODES = ("CONST", "EXACT", "POPOVAS")
# the legacy approximations removed in KROME-X 2.0.0: each is EXACT with
# a degree of freedom frozen or a species left out
GAMMA_REMOVED = {
    "FULL": "H/He/e- and H2 only, H2 always 7/5 (rotation fully excited, "
            "no vibration)",
    "VIB": "classical rotation",
    "ROT": "no vibration",
    "REDUCED": "no vibration, diatomics H2 and CO only"}
# Weingartner & Draine 2001 (ApJS 134, 263) photoelectric heating, eq. (44)
# (Table 2: C0..C6) and grain recombination cooling, eq. (45) (Table 3:
# D0..D4), per (R_V, b_C [1e-5 C/H in the smallest grains], size-
# distribution case, radiation field: B0 = a 3e4 K blackbody, ISRF =
# Mathis+ 1983); --wd01 selects the row
WD01_TABLE = {
    (3.1, 0.0, "A", "B0"): ((5.56, 1.82e-2, 0.00492, 0.03368, 0.557, 0.666,
                             0.532),
                            (0.4800, 1.783, -7.617, 1.655, 0.06326)),
    (3.1, 2.0, "A", "B0"): ((3.41, 4.27, 0.03700, 0.00569, 0.102, 0.494,
                             0.669),
                            (0.5322, 1.325, -6.047, 1.435, 0.05378)),
    (3.1, 4.0, "A", "B0"): ((10.59, 4.05e-3, 0.01234, 0.01391, 0.808, 0.580,
                             0.573),
                            (0.4336, 2.108, -7.359, 1.736, 0.06502)),
    (3.1, 6.0, "A", "B0"): ((3.66, 7.66, 0.00661, 0.01552, 0.094, 0.697,
                             0.482),
                            (0.4270, 2.120, -7.301, 1.786, 0.06731)),
    (4.0, 0.0, "A", "B0"): ((4.30, 1.79e-4, 0.00572, 0.02488, 1.026, 0.701,
                             0.505),
                            (0.5332, 1.379, -7.125, 1.504, 0.05732)),
    (4.0, 2.0, "A", "B0"): ((5.03, 2.27, 0.07441, 0.00345, 0.140, 0.412,
                             0.737),
                            (0.4425, 2.019, -7.755, 1.712, 0.06424)),
    (4.0, 4.0, "A", "B0"): ((3.28, 4.44, 0.00786, 0.01219, 0.102, 0.686,
                             0.500),
                            (0.4847, 1.493, -6.109, 1.528, 0.05751)),
    (5.5, 0.0, "A", "B0"): ((3.67, 2.42e-3, 0.07034, 0.00360, 0.725, 0.365,
                             0.814),
                            (0.4184, 2.129, -9.082, 1.818, 0.06943)),
    (5.5, 1.0, "A", "B0"): ((5.45, 2.44e-2, 0.09594, 0.00225, 0.534, 0.418,
                             0.758),
                            (0.4981, 1.510, -7.237, 1.535, 0.05777)),
    (5.5, 2.0, "A", "B0"): ((4.87, 7.28e-3, 0.00812, 0.01969, 0.656, 0.620,
                             0.534),
                            (0.4851, 1.501, -6.775, 1.533, 0.05762)),
    (5.5, 3.0, "A", "B0"): ((5.93, 9.26e-2, 0.00675, 0.01648, 0.417, 0.683,
                             0.489),
                            (0.4948, 1.442, -6.377, 1.505, 0.05645)),
    (4.0, 0.0, "B", "B0"): ((4.85, 6.59e-2, 0.00437, 0.04234, 0.415, 0.664,
                             0.508),
                            (0.4568, 1.876, -7.810, 1.688, 0.06412)),
    (4.0, 2.0, "B", "B0"): ((5.97, 0.287, 0.00736, 0.02065, 0.302, 0.637,
                             0.528),
                            (0.4597, 1.594, -6.675, 1.599, 0.06094)),
    (4.0, 4.0, "B", "B0"): ((7.97, 0.175, 0.00277, 0.03820, 0.382, 0.771,
                             0.401),
                            (0.5155, 1.340, -5.791, 1.449, 0.05411)),
    (5.5, 0.0, "B", "B0"): ((4.04, 0.653, 0.11310, 0.00259, 0.198, 0.348,
                             0.809),
                            (0.4360, 1.996, -8.372, 1.735, 0.06572)),
    (5.5, 1.0, "B", "B0"): ((3.38, 0.608, 0.00396, 0.02568, 0.192, 0.752,
                             0.439),
                            (0.5140, 1.349, -6.676, 1.472, 0.05536)),
    (5.5, 2.0, "B", "B0"): ((1.10, 3.30, 0.00941, 0.01431, 0.087, 0.643,
                             0.524),
                            (0.4633, 1.690, -7.158, 1.607, 0.06055)),
    (5.5, 3.0, "B", "B0"): ((3.47, 1.77, 0.00494, 0.01953, 0.140, 0.723,
                             0.456),
                            (0.4677, 1.780, -7.173, 1.615, 0.06027)),
    (3.1, 0.0, "A", "ISRF"): ((4.35, 2.63e-5, 0.00242, 0.03003, 1.235,
                               0.827, 0.399),
                              (0.4291, 2.406, -8.357, 1.714, 0.06354)),
    (3.1, 2.0, "A", "ISRF"): ((4.87, 0.948, 0.02576, 0.00766, 0.188, 0.490,
                               0.650),
                              (0.5232, 1.678, -5.942, 1.339, 0.04813)),
    (3.1, 4.0, "A", "ISRF"): ((7.41, 0.772, 0.03895, 0.00606, 0.239, 0.408,
                               0.708),
                              (0.3959, 2.380, -6.554, 1.575, 0.05674)),
    (3.1, 6.0, "A", "ISRF"): ((9.30, 0.248, 0.00697, 0.01848, 0.365, 0.633,
                               0.509),
                              (0.3632, 2.937, -7.601, 1.742, 0.06228)),
}
WD01_DEFAULT = "3.1,6.0,A,ISRF"  # Milky Way dust, the solar-neighbourhood ISRF


def parse_wd01(spec):
    """'RV,bC,case,field' -> the WD01_TABLE key"""
    t = [x.strip() for x in (spec or WD01_DEFAULT).split(",")]
    if len(t) != 4:
        raise ValueError("thermal.wd01 takes RV,bC,case,field (e.g. %s)"
                         % WD01_DEFAULT)
    try:
        key = (float(t[0]), float(t[1]), t[2].upper(), t[3].upper())
    except ValueError:
        raise ValueError("thermal.wd01 %s: RV and bC are numbers" % spec)
    if key not in WD01_TABLE:
        raise ValueError("thermal.wd01 %s: not a row of WD01 Table 2 (%s)" % (
            spec, "; ".join("%g,%g,%s,%s" % k for k in WD01_TABLE)))
    return key


# parameters read by the runtime modules for a configured thermal process
PROCESS_PARAMS = {
    "heat_compress": ("tff",),
    "heat_visc": ("visc_nu", "visc_omega"),
    "cool_expansion": ("hubble0", "omega0", "zredshift"),
    "cool_compton": ("zredshift",),
    "cool_CIE": ("Tcmb",),
    "cool_HD": ("Tcmb",),
    "heat_photoDust": ("G0", "Av"),
    "heat_netPhotoDust": ("G0", "Av"),
    "cool_dust": ("d2g",),
    "cool_dust_bins": ("d2g", "Tcmb"),
    "cooling_floor": ("Tfloor",),
    "dust_evol": ("d2g", "Tcmb"),
}
# ... and by the runtime extern functions a network calls
EXTERN_PARAMS = {
    "elec_recomb_ST93": ("zredshift", "hubble0", "omega0"),
    **{"recfast_%s" % x: ("zredshift", "hubble0", "omega0")
       for x in ("H_rec", "H_ion", "He_rec", "He_ion")},
    "krate_stickSi": ("d2g",),
    "shield_dust": ("d2g",),
    "krate_H2dust": ("d2g",),
    "krate_H2dust_Jura": ("d2g", "clump"),
    **{"%s_recombination_on_dust" % x: ("d2g", "G0")
       for x in ("H", "He", "C", "Si", "O", "Fe")},
}
# consumers of the live dust temperature table_Tdust (the Tdust parameter
# unless --dust-tables supplies it)
TDUST_EXTERNS = ("krate_stickSi", "krate_evaporation", "get_table_Tdust",
                 "krate_H2dust")
# H2 formation on dust rates: CHEM heating adds 4.2 h2heatfac + 0.2 eV per
# H2 formed (legacy HChem nH2dust term)
H2DUST_EXTERNS = ("krate_H2dust", "krate_H2dust_Jura", "krate_H2dust_table")
TDUST_PROCESSES = ("cool_dust",)


def extract_cie_tables(data_root):
    """the Shen et al. 2013 ZCIE tables ship compressed
    (data/coolZ_CIE.tar.xz, as in legacy): extracted on first use next to
    the archive in a checkout, or, for an installed package (or a
    read-only folder), into the user cache ($XDG_CACHE_HOME or ~/.cache,
    /kromex), so that uninstalling leaves nothing behind; returns the
    folder holding them (legacy's tar also had coolZ_CIE2005.dat, which
    nothing read)"""
    import os
    import tarfile
    data_root = Path(data_root)
    arc = data_root / "coolZ_CIE.tar.xz"
    names = ("coolZ_CIE2012.dat", "coolZ_CIE2012NOUV.dat")
    cache = Path(os.environ.get("XDG_CACHE_HOME",
                                Path.home() / ".cache")) / "kromex"
    for d in (data_root, cache):
        if all((d / n).exists() for n in names):
            return d
    installed = "_share" in data_root.resolve().parts
    dest = cache if installed or not os.access(data_root, os.W_OK) \
        else data_root
    dest.mkdir(parents=True, exist_ok=True)
    if arc.exists():
        with tarfile.open(arc, mode="r:xz") as tf:
            for m in tf.getmembers():
                if m.name in names and not (dest / m.name).exists():
                    try:  # Python >= 3.12: the safe extraction filter
                        tf.extract(m, dest, filter="data")
                    except TypeError:
                        tf.extract(m, dest)
    return dest


def resolve_data_path(path, data_dir):
    """a data file given by the user: as given (absolute, or relative to
    the working directory), else relative to the root of v2/, its
    data/cooling, the legacy KROME root or its data directory"""
    p = Path(path)
    if p.is_absolute() or p.exists():
        return p
    root = paths.ROOT
    for base in (root, root / "data" / "cooling", root.parent / "legacy",
                 Path(data_dir).parent):
        if (base / p).exists():
            return base / p
    raise FileNotFoundError("coolant file %s not found (looked in the "
                            "working directory, %s, %s)"
                            % (path, root, root.parent / "legacy"))


def load_photo_db(path):
    """Legacy photoionization database -> list of entries with canonical
    species names, energy limits [eV] and the rate string: sigma_v96(...)
    (Verner et al. 1996, the outer shell), or above the inner-shell edge
    a sum of sigma_vy95(...) (Verner & Yakovlev 1995, one per shell)."""
    def canon(name):
        n = name.strip()
        return "e-" if n.upper() == "E" else n
    entries, cur = [], {}
    if not path.exists():
        return entries
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("@type:"):
            cur = {"type": s.split(":", 1)[1].strip()}
        elif s.startswith("@reacts:"):
            cur["reacts"] = sorted(canon(x)
                                   for x in s.split(":", 1)[1].split(","))
        elif s.startswith("@prods:"):
            cur["prods"] = sorted(canon(x)
                                  for x in s.split(":", 1)[1].split(","))
        elif s.startswith("@limits:"):
            cur["limits"] = [float(d2e(x))
                             for x in s.split(":", 1)[1].split(",")]
        elif s.startswith("@rate:"):
            cur["rate"] = s.split(":", 1)[1].strip()
            entries.append(cur)
    return entries


# photo cross sections from data files (the generator and kromex_doc)

# h c [eV Angstrom]: Angstrom -> eV (SWRI files); the Leiden files give
# the wavelength in nm (their headers; legacy kromelib LEIDEN2KROME)
HC_EV_AA = constants.clight * constants.planck_eV / 1e-8
HC_EV_NM = constants.clight * constants.planck_eV / 1e-7

def leiden_table(folder, rx):
    """Leiden database file <R>__<P1>_<P2>.dat -> (E[eV], sigma[cm^2]),
    ascending in E. Column: photoionisation if e- among the products,
    photodissociation otherwise."""
    rea = rx["reactants"][0]["sp"]
    prods = sorted(p["sp"] for p in rx["products"])
    folder = folder / "leiden_xsecs"
    fname = None
    for f in sorted(folder.glob("%s__*.dat" % rea)):
        toks = f.stem.split("__", 1)[1].split("_")
        fp = sorted("e-" if t.upper() == "E" else t for t in toks)
        if fp == prods:
            fname = f
            break
    if fname is None:
        raise FileNotFoundError(
            "no Leiden file %s__*.dat matching products %s in %s "
            "(download from the Leiden photodissociation database)"
            % (rea, "+".join(prods), folder))
    col = 3 if "e-" in prods else 2
    E, S = [], []
    for line in fname.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        v = s.split()
        E.append(HC_EV_NM / float(v[0]))  # lambda [nm] -> eV
        S.append(float(v[col]))
    pairs = sorted(zip(E, S))
    return ([e for e, _ in pairs], [s for _, s in pairs])

def swri_table(folder, rx, warnings):
    """SWRI (phidrates) file <R>.dat -> (E[eV], sigma[cm^2]) for the
    branch matching the reaction products. Branch header tokens use '/'
    between products and '+' for 'ion + e-'; excited-state labels are
    merged onto the ground state (columns summed), as legacy KROME did.
    Spectral lines (1e-35 sentinel brackets) become boxes, see below."""
    rea = rx["reactants"][0]["sp"]
    prods = sorted(p["sp"] for p in rx["products"])
    fname = folder / "swri_xsecs" / ("%s.dat" % rea)
    if not fname.exists():
        raise FileNotFoundError(
            "no SWRI file %s (download from phidrates.space.swri.edu)"
            % fname)
    header, rows = [], []
    for line in fname.read_text().splitlines():
        toks = line.split()
        if not toks:
            continue
        if toks[0] == "Lambda":
            header, rows = toks, []
            continue
        try:
            rows.append([float(x) for x in toks])
        except ValueError:
            continue
    for state in ("1D", "1S", "3P"):
        header = [h.replace(state, "") for h in header]
    # merge duplicate branches (excited states folded to ground)
    drop = set()
    for i in range(1, len(header)):
        for j in range(i + 1, len(header)):
            if header[i] == header[j] and j not in drop:
                drop.add(j)
                for row in rows:
                    row[i] += row[j]
    branch_prods = []
    for h in header:
        ps = [x for x in h.replace("+", "+/E/").split("/") if x]
        branch_prods.append(sorted("e-" if x.upper() == "E" else x
                                   for x in ps))
    if prods not in branch_prods:
        raise FileNotFoundError(
            "products %s not among SWRI branches %s in %s"
            % ("+".join(prods),
               branch_prods[2:] if len(branch_prods) > 2 else [], fname))
    col = branch_prods.index(prods)
    # 1e-35 rows are SWRI's spectral-line sentinels. A line is a run
    # of values bracketed by sentinel runs on both sides (H2: the
    # Lyman-Werner lines, "1e-35 v 1e-35" triplets): a box of that
    # height over the bracketing sentinel energies, zero outside (as
    # legacy SWRI2KROME; the value mean where a line has several). Any
    # other sentinel (C2: 1436 A; CO, OH, H2O: the last two rows; H2:
    # the outer edges of the line block) marks the edge of the
    # tabulated range and is dropped. Unlike legacy, a continuum in
    # the same branch is kept (legacy replaced it by the lines).
    asc = [r for r in rows[::-1] if r[0] > 0.0]  # ascending energy
    pts = [(HC_EV_AA / r[0], r[col]) for r in asc]
    sent = [v == 1e-35 for _, v in pts]
    pairs, n_lines, q = [], 0, 0
    while q < len(pts):
        if sent[q]:
            q += 1
            continue
        e = q
        while e < len(pts) and not sent[e]:
            e += 1
        vals = pts[q:e]
        if q > 0 and e < len(pts):  # bracketed: a line box
            lo, hi = pts[q - 1][0], pts[e][0]
            vv = sum(v for _, v in vals) / len(vals)
            pairs += [(lo * (1.0 - 1e-9), 0.0), (lo, vv), (hi, vv),
                      (hi * (1.0 + 1e-9), 0.0)]
            n_lines += 1
        else:
            pairs += vals
        q = e
    if n_lines:
        warnings.append("SWRI %s: %d spectral lines embedded as "
                             "finite-width boxes" % (fname.name, n_lines))
    # exact duplicate energies (a continuum point on a box edge) would
    # make a zero-width interpolation interval: keep the first
    pairs.sort()
    dedup = []
    for e, v in pairs:
        if dedup and e == dedup[-1][0]:
            continue
        dedup.append((e, v))
    pairs = dedup
    # the table spans its first to last nonzero point (legacy xmin /
    # xmax): a threshold is a step, not a linear ramp from the last
    # tabulated zero (OH -> OH+: zeros below 19.8 eV); xsec_interp is
    # zero outside the range
    while pairs and pairs[0][1] == 0.0:
        pairs.pop(0)
    while pairs and pairs[-1][1] == 0.0:
        pairs.pop()
    return ([e for e, _ in pairs], [s for _, s in pairs])


# legacy rate-database helper functions (src/krome_subs.f90), emitted
# verbatim into the generated header when a resolved auto rate uses them.
# All are pure functions of Tgas and constants.
AUTO_TFUNCS = {
    "colion_v96": """static inline double colion_v96(double Tgas, double dE, double P,
                                double A, double X, double K) {
  // collisional ionization, Voronov 1997 fit (legacy used kB = 8.621738e-5
  // eV/K, 5e-4 high)
  const double U = dE / (Tgas * constants::boltzmann_eV);
  return A * (1.0 + P * std::sqrt(U)) * std::pow(U, K) * std::exp(-U)
      / (X + U);
}""",
    "recV96": """static inline double recV96(double Tgas, double a, double b) {
  // radiative recombination, Verner & Ferland power law
  return a * std::pow(1e4 / Tgas, b);
}""",
    "recNewV96": """static inline double recNewV96(double Tgas, double r1, double r2,
                               double r3, double r4) {
  // radiative recombination, Verner & Ferland 1996 fit
  const double tt = std::sqrt(Tgas / r3);
  return r1 / (tt * std::pow(tt + 1.0, 1.0 - r2)
               * std::pow(1.0 + std::sqrt(Tgas / r4), 1.0 + r2));
}""",
    "recFeV96": """static inline double recFeV96(double Tgas, double r1, double r2,
                              double r3) {
  // radiative recombination of Fe ions (legacy recFeV96)
  const double tt = std::sqrt(Tgas * 1e-4);
  return r1 / std::pow(tt, r2 + r3 + std::log10(tt));
}""",
}


def load_auto_db(folder):
    """Legacy rate database (data/database/*.dat: hydroHe.dat and
    recombination_badnell.dat first, the rest in sorted order, as legacy):
    entries {reacts, prods (sorted KROME-X names), limits
    (lo, hi) or None, rate (legacy expression string), fname}."""
    import convert_legacy as cl
    entries = []
    if not folder.is_dir():
        return entries
    # hydroHe.dat (the primordial H/He set) and Badnell's recombination
    # (2.0.0: before the astrochemistry databases, whose X+ + e- -> X fits
    # lack the dielectronic part) first, then the rest in sorted order
    first = ("hydroHe.dat", "recombination_badnell.dat")
    files = [folder / n for n in first if (folder / n).is_file()]
    files += [f for f in sorted(folder.iterdir()) if f.name not in first]
    for f in files:
        if not f.is_file() or f.suffix != ".dat" or "~" in f.name:
            continue
        cur = None
        for line in f.read_text(errors="replace").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s == "#BREAK DATABASE":
                break
            low = s.lower()
            if (low.startswith("@photoxsec:") or low.startswith("@cr:")
                    or low.startswith("@photoav:")):
                break
            if not s.startswith("@") or ":" not in s:
                continue
            key, val = s[1:].split(":", 1)
            key = key.strip().lower()
            val = val.strip()
            if key == "type":
                cur = {"limits": None, "fname": f.name}
            elif cur is None:
                continue
            elif key in ("reacts", "prods"):
                names = [cl.canon_species(x) for x in val.split(",")]
                cur[key] = sorted(x for x in names if x)
            elif key == "limits":
                # comparison prefixes (.LE., >, ...) as in legacy readTOpt;
                # KROME-X ranges are half-open, which the strictness does not
                # change for a continuous tiling
                parts = [cl.parse_limit(x) for x in val.split(",")]
                if len(parts) == 2 and any(parts):  # empty = unlimited
                    lo, hi = parts
                    cur["limits"] = (float(lo) if lo else 0.0,
                                     float(hi) if hi else math.inf)
            elif key == "rate":
                cur["rate"] = val
                if "reacts" in cur and "prods" in cur:
                    entries.append(cur)
                cur = None
    return entries


def d2e(s):
    return re.sub(r"(?<=[0-9.])[dD](?=[+\-]?[0-9])", "e", s)


class Codegen:
    # outside-of-range behavior for piecewise T groups (network-format.md 3.1):
    #  'asym' (default): zero below the lower bound, clamp above the upper
    #  'clamp': clamp at both bounds; 'zero': legacy KROME (reaction off)
    # secondary ionization follows Quadri et al. 2025 (arXiv:2505.05556,
    # Appendix B): Dalgarno+99 results refit in the Shull & van Steenberg
    # 1985 form, including the H/H2/He mixture. Fit species and the
    # threshold energies used in the phi_i(E) = (E/Eth - 1) A_i fits
    # (paper eqs. B.10-B.12); metal species are rescaled from H by their
    # X-ray cross-section ratio sigma_A/sigma_H (eq. B.14), given in the
    # config as e.g. C=3.92 (Meijerink & Spaans 2005 value, paper B.15).
    Q25_FIT = {"H": (0, 13.6), "H2": (1, 15.4), "He": (2, 24.6)}
    Q25_KIND_METAL = 3
    SECION_RATIOS_FILE = paths.ROOT / "data" \
        / "secion_xray_ratios.dat"

    def load_secion_ratios(self):
        """Metal secondary-ionization scaling ratios sigma_ei(A)/sigma_ei(H)
        (Meijerink & Spaans 2005 table, see the data file header)."""
        out = {}
        if self.SECION_RATIOS_FILE.exists():
            for line in self.SECION_RATIOS_FILE.read_text().splitlines():
                s = line.split("#")[0].split()
                if len(s) == 2:
                    out[s[0]] = float(s[1])
        return out

    # molecular table coolants (legacy -cooling=CO/OH/H2O/HCN): trilinear
    # tables in (log10 Tgas, log10 (n_H + n_H2), log10 N_coolant) following
    # Neufeld & Kaufman 1993 eqn. 5, data courtesy of K. Omukai (legacy
    # data/cool*.dat, loaded verbatim and embedded in the header)
    MOL_COOLANTS = {"CO": "coolCO.dat", "OH": "coolOH.dat",
                    "H2O": "coolH2O.dat", "HCN": "coolHCN.dat"}

    def __init__(self, net, outside="asym", data_dir=DEFAULT_DATA_DIR,
                 secion=(), mol_cooling=(), num2col="default",
                 metal_cooling=(), zcie="", dust_tables="",
                 gas_total=False, dh=False, cooling=(), heating=(),
                 cooling_floor=False, evolve_rad=False, evolve_Td=False,
                 shielding="", gamma="CONST", h2_opacity="",
                 photo_induced=False, photoelectric="BT94", wd01="",
                 cmb_pumping=False, habing_shield=False,
                 reactant_multiplicity=False, screening=False,
                 binding_energies="KROME", surface_nu="HASEGAWA",
                 surface_competition=True, surface_tunnelling=True,
                 reactive_desorption="NONE",
                 cr_desorption="", surface_bins=0, rtol=None, atol=None,
                 cool_files=(), cool_levels=None, evolve_dust="",
                 custom_atol="", custom_rtol="", cool_escape="",
                 rhs="auto"):
        self.net = net
        # --rhs: the mass-action rhs, fluxes and Jacobian as generated
        # statements (unrolled) or as loops over index tables (loop; legacy
        # -iRHS), auto: loops above RHS_LOOP_MIN reactions
        if rhs not in ("auto", "unrolled", "loop"):
            raise ValueError("rates.rhs must be auto, unrolled or loop")
        self.rhs_mode = rhs
        # --cool-escape: the optical depth of the line coolants (LVG
        # escape probability; "" = optically thin)
        self.cool_escape = (cool_escape or "").upper()
        # project defaults of the solver tolerances (--rtol / --atol,
        # legacy -RTOL / -ATOL): the starting values of SolverOpts::rtol /
        # atol_n, still changeable at runtime
        self.rtol = 1e-8 if rtol is None else float(rtol)
        self.atol = 1e-20 if atol is None else float(atol)
        if not (self.rtol > 0.0 and self.atol > 0.0):
            raise ValueError("solver.rtol and solver.atol must be positive")
        self.surface_bins = surface_bins  # model B: bins (0: model A)
        self.binding_energies = binding_energies  # KROME | SEMENOV
        self.surface_nu = surface_nu              # HASEGAWA | LEGACY
        self.surface_competition = surface_competition
        # barrier reactions: max(thermal, tunnelling) (default), or thermal
        # hopping only (legacy Semenov framework)
        self.surface_tunnelling = surface_tunnelling
        self.reactive_desorption = reactive_desorption  # NONE|FIXED|MINISSALE
        self.cr_desorption = cr_desorption        # CR-rate param name or ""
        self._gen_ext = set()  # per-network generated extern functions
        self.reactant_multiplicity = reactant_multiplicity
        self.screening = screening
        self.outside = outside
        self.secion = [s for s in secion if s]
        # 'auto' = every known coolant present in the network
        if mol_cooling == "auto" or "auto" in mol_cooling:
            mol_cooling = [c for c in self.MOL_COOLANTS if c in net.species]
        self.mol_cooling = [c for c in mol_cooling if c]
        self.num2col = num2col  # 'default' | 'jeans' | 'jeans40' | 'length'
        # metal line cooling: list of coolZ.dat metals, or 'auto' = all
        # whose species exist in the network
        if metal_cooling == "auto" or "auto" in (metal_cooling or ()):
            metal_cooling = "auto"
        else:
            metal_cooling = [m for m in (metal_cooling or ()) if m]
        self.metal_cooling = metal_cooling
        self.zcie = zcie  # '' | '2012' (HM12 UV + heating) | 'nouv'
        self.dust_tables = dust_tables  # '' | 'HM2012' | 'MC'
        self.dh = bool(dh)  # enthalpy heating/cooling (legacy DH)
        # gas+total representation (docs/solver-design.md section 4): each
        # surface slot of a gas<->surface pair stores TOTAL = gas + surface,
        # removing the stiff freezeout/evaporation exchange from it
        self.gas_total = bool(gas_total)
        self.pair_gas, self.total_of = {}, {}
        if self.gas_total:
            for s in sorted(net.species):
                if net.species[s]["phase"] != "gas":
                    base = s.rsplit("_", 1)[0]
                    if base in net.species:
                        self.pair_gas[s] = base
                        self.total_of[base] = s
        self.data_dir = Path(data_dir)
        self.photo_db = load_photo_db(self.data_dir /
                                      "photoionization.dat")
        self.warnings = []
        # per-species defaults (--custom-atol / --custom-rtol, legacy
        # -customATOL / -customRTOL): {species: value}, the starting
        # SolverOpts::atol_species / rtol_species
        self.custom_atol = self.load_species_tols(custom_atol, "atol")
        self.custom_rtol = self.load_species_tols(custom_rtol, "rtol")
        # --cmb-pumping: needed by configure_floor (the pumped metal lines
        # drop out of the cooling floor)
        self.cmb_pumping = bool(cmb_pumping)
        # coolant files (--cool-file, legacy -coolFile) and the level cut
        # (--cool-levels, legacy -coolLevels, every coolant)
        if cool_levels is not None and int(cool_levels) < 1:
            raise ValueError("thermal.levels needs at least 1 (levels 0-1)")
        self.cool_levels = None if cool_levels is None else int(cool_levels)
        self.file_coolants = self.load_cool_files(cool_files or ())
        self.file_active = []
        # thermal processes: fixed at generation time (legacy -cooling= /
        # -heating=); emitted as the constexpr ThermalConfig
        self.thermal = self.configure_thermal(cooling, heating,
                                              cooling_floor)
        # krome_fshield() model (legacy -shielding): WG11 by default since
        # 2.0.0 (legacy's default was none, f = 1, which switched off the
        # shielding a network asks for); NONE keeps f = 1
        shielding = (shielding or "WG11").upper()
        if shielding not in ("NONE", "DB96", "WG11", "R14"):
            raise ValueError("unknown radiation.shielding %s (WG11, DB96, R14, NONE)"
                             % shielding)
        self.shielding = "" if shielding == "NONE" else shielding
        gamma = (gamma or "CONST").upper()
        if gamma in GAMMA_REMOVED:
            raise ValueError("thermal.gamma=%s (%s) was removed in KROME-X "
                             "2.0.0: use EXACT (all the diatomics with "
                             "their exact rotation and vibration) or "
                             "POPOVAS" % (gamma, GAMMA_REMOVED[gamma]))
        if gamma not in GAMMA_MODES:
            raise ValueError("unknown thermal.gamma %s (%s)"
                             % (gamma, ", ".join(GAMMA_MODES)))
        self.gamma = gamma
        h2_opacity = (h2_opacity or "").upper()
        if h2_opacity not in ("", "RIPAMONTI", "OMUKAI"):
            raise ValueError("unknown thermal.h2_opacity %s (RIPAMONTI, OMUKAI)"
                             % h2_opacity)
        self.h2_opacity = h2_opacity
        self.photo_induced = bool(photo_induced)
        # the dust photoelectric heating model of PHOTODUST / PHOTODUSTNET
        photoelectric = (photoelectric or "BT94").upper()
        if photoelectric not in ("BT94", "WD01"):
            raise ValueError("unknown thermal.photoelectric %s (BT94, WD01)"
                             % photoelectric)
        self.photoelectric = photoelectric
        self.wd01 = parse_wd01(wd01)
        if wd01 and photoelectric != "WD01":
            raise ValueError("thermal.wd01 selects a row of the WD01 fit: it "
                             "needs thermal.photoelectric = 'WD01'")
        self.habing_shield = bool(habing_shield)
        self.thermal["photo_induced"] = self.photo_induced
        # what is evolved alongside the chemistry: fixed at generation too
        # (thermal evolution as a whole stays the one runtime switch)
        self.evolve_rad = bool(evolve_rad)
        self.evolve_Td = bool(evolve_Td)
        # --evolve-dust: grain radii of the bins as solver state, changed
        # by accretion of the key element (GROWTH), thermal sputtering
        # (SPUTTER) and sublimation (EVAP, with KELVIN its curvature
        # correction)
        ev = [x.strip().upper() for x in (evolve_dust or "").split(",")
              if x.strip()]
        for x in ev:
            if x not in ("GROWTH", "SPUTTER", "EVAP", "KELVIN"):
                raise ValueError("dust.evolve: unknown process %s "
                                 "(GROWTH, SPUTTER, EVAP, KELVIN)" % x)
        if "KELVIN" in ev and "EVAP" not in ev:
            raise ValueError("dust.evolve: KELVIN corrects the "
                             "evaporation, it needs EVAP")
        self.evolve_dust = ev
        if ev and not ({"C", "Si"} & set(net.species)):
            raise ValueError("dust.evolve needs the gas-phase key "
                             "element of the grains in the network: C "
                             "(carbonaceous) and/or Si (silicates)")
        self.thermal["dust_evol"] = bool(ev)
        if self.evolve_rad and not any(rx["rate"]["kind"] == "photo"
                                       for rx in net.reactions):
            raise ValueError("radiation.evolve needs photo reactions in the "
                             "network (the field is absorbed by them)")
        if self.evolve_Td and not self.thermal["cool_dust_bins"]:
            raise ValueError("dust.evolve_Td needs thermal.cooling = ['DUSTBINS'] (the "
                             "per-bin dust temperatures it evolves)")
        self.spidx = {s: i for i, s in enumerate(sorted(net.species))}
        # does an AST (transitively through vars) depend on n()?
        self._ndep_cache = {}
        # vars hoisted into locals in the current emission scope: name->local
        self._hoisted = {}
        self._jlw_unit = False  # print J_LW as 1 (the unit-J_LW rates)
        self._in_ratefn = False  # printing a rate_<i>() body

    # -- thermal process configuration ------------------------------------------

    # the ThermalConfig members, in emission order
    THERMAL_FLAGS = (
        "cool_atomic", "cool_H2", "cool_HD", "cool_compton", "cool_ff",
        "cool_CIE", "cool_expansion", "cool_cont", "cool_mol", "cool_Z",
        "cool_ZCIE", "cool_Zext", "cool_dust", "cool_dustTab",
        "cool_dust_bins", "cool_dH", "heat_compress", "heat_visc",
        "heat_photo", "heat_CR", "heat_chem", "heat_dH", "heat_photoDust",
        "heat_netPhotoDust", "heat_photoAv", "heat_ZCIE", "cooling_floor",
        "floor_atomic", "floor_H2", "floor_HD", "floor_mol", "floor_Z",
        "floor_ZCIE", "floor_Zext", "cool_chem", "photo_induced",
        "heat_nuclear", "cool_file", "floor_file", "dust_evol")
    # --cooling-floor terms (legacy -useIndividualFloor names): the
    # radiative line coolers, each floored as cool(T) - cool(Tfloor)
    FLOOR_NAMES = {
        "ATOMIC": "floor_atomic", "H2": "floor_H2", "HD": "floor_HD",
        "CO": "floor_mol", "OH": "floor_mol", "H2O": "floor_mol",
        "HCN": "floor_mol", "Z": "floor_Z", "Z_CIE": "floor_ZCIE",
        "Z_CIENOUV": "floor_ZCIE", "Z_EXTENDED": "floor_Zext"}
    COOLING_NAMES = {  # legacy -cooling= names -> ThermalConfig flag
        "ATOMIC": "cool_atomic", "H2": "cool_H2", "HD": "cool_HD",
        "COMPTON": "cool_compton", "FF": "cool_ff", "CIE": "cool_CIE",
        "EXPANSION": "cool_expansion", "CONT": "cool_cont",
        "CHEM": "cool_chem", "DISS": "cool_chem",
        "Z": "cool_Z", "Z_CIE": "cool_ZCIE",
        "Z_CIENOUV": "cool_ZCIE", "Z_EXTENDED": "cool_Zext",
        "CO": "cool_mol", "OH": "cool_mol", "H2O": "cool_mol",
        "HCN": "cool_mol", "DUST": "cool_dust", "DUSTTAB": "cool_dustTab",
        "DUSTBINS": "cool_dust_bins", "DH": "cool_dH"}
    HEATING_NAMES = {  # legacy -heating= names -> ThermalConfig flag
        "COMPRESS": "heat_compress", "VISCOUS": "heat_visc",
        "PHOTO": "heat_photo", "CR": "heat_CR", "CHEM": "heat_chem",
        "DH": "heat_dH", "PHOTOAV": "heat_photoAv",
        "PHOTODUST": "heat_photoDust", "PHOTODUSTNET": "heat_netPhotoDust",
        "NUCLEAR": "heat_nuclear"}
    UNSUPPORTED = {"GH", "XRAY"}
    # processes removed in KROME-X 2.0.0, with what to use instead
    REMOVED = {"H2GP98": "H2GP98 (Galli & Palla 1998) was removed in "
               "KROME-X 2.0.0: use H2 (Glover & Abel 2008, with the same LTE "
               "limit below 2e3 K; legacy -cooling=H2)"}

    def configure_thermal(self, cooling, heating, cooling_floor):
        """validate the thermal.cooling/thermal.heating lists, derive the table
        generation they imply (molecular, metal, ZCIE tables, enthalpies)
        and return the ThermalConfig flags {name: bool}."""
        def names(lst, what):
            out = [x.strip().upper() for x in lst if x and x.strip()]
            for x in out:
                if x in self.UNSUPPORTED:
                    raise ValueError("%s %s is not supported in KROME-X"
                                     % (what, x))
                if x in self.REMOVED:
                    raise ValueError("%s %s" % (what, self.REMOVED[x]))
            return out
        cool = names(cooling, "cooling")
        heat = names(heating, "heating")
        # coolant-file coolants: selected by their flag (legacy names: CII,
        # FeII for atoms and ions, the species for molecules, e.g. O2), by
        # FILE (all of them) or by --metal-cooling (a metal a file defines
        # comes from the file); each is its own cooling term
        byflag = {md["flag"].upper(): md for md in self.file_coolants}
        sel, rest = [], []
        for x in cool:
            if x == "FILE":
                if not self.file_coolants:
                    raise ValueError("cooling FILE needs thermal.coolant_files")
                sel += self.file_coolants
            elif x in byflag:
                sel.append(byflag[x])
            else:
                rest.append(x)
        cool = rest
        if self.metal_cooling and self.metal_cooling != "auto":
            for md in self.file_coolants:
                if md["species"] in self.metal_cooling and md not in sel:
                    self.warnings.append(
                        "metal coolant %s taken from %s (the user file "
                        "overrides the embedded coolZ.dat)"
                        % (md["species"], md["file"]))
                    sel.append(md)
        for md in sel:
            if md in self.file_active:
                continue
            if md["species"] is None:
                self.warnings.append("coolant %s (%s): species %s not in "
                                     "the network, skipped"
                                     % (md["flag"], md["file"], md["name"]))
            else:
                self.file_active.append(md)
        if "ZCIE" in heat:
            raise ValueError(
                "heating ZCIE is not a separate option: the Shen+13 ZCIE "
                "heating always comes with its cooling (thermal.cooling = ['Z_CIE'], "
                "or Z_EXTENDED above 1e4 K)")
        for x in cool:
            if x not in self.COOLING_NAMES:
                raise ValueError("unknown cooling %s (known: %s)" % (
                    x, ", ".join(sorted(self.COOLING_NAMES))))
        for x in heat:
            if x not in self.HEATING_NAMES:
                raise ValueError("unknown heating %s (known: %s)" % (
                    x, ", ".join(sorted(self.HEATING_NAMES))))
        flags = {f: False for f in self.THERMAL_FLAGS}
        for x in cool:
            flags[self.COOLING_NAMES[x]] = True
        for x in heat:
            flags[self.HEATING_NAMES[x]] = True
        flags["cool_file"] = bool(self.file_active)

        # implied generation
        mols = [x for x in cool if x in self.MOL_COOLANTS]
        for m in mols:
            if m not in self.mol_cooling:
                self.mol_cooling.append(m)
        if ("Z" in cool or "Z_EXTENDED" in cool) and not self.metal_cooling:
            self.metal_cooling = "auto"
        if "Z_CIE" in cool and "Z_CIENOUV" in cool:
            raise ValueError("cooling Z_CIE and Z_CIENOUV are exclusive")
        if "Z_CIENOUV" in cool:
            self.zcie = "nouv"
        elif "Z_CIE" in cool or "Z_EXTENDED" in cool:
            self.zcie = self.zcie or "2012"
        # the ZCIE cooling and heating (UV background, 2012 tables) always
        # go together; the NOUV table has no heating
        if "Z_CIE" in cool or "Z_EXTENDED" in cool:
            flags["heat_ZCIE"] = True
        if "DH" in cool or "DH" in heat:
            self.dh = True
        if "DUSTTAB" in cool and not self.dust_tables:
            raise ValueError("cooling DUSTTAB needs dust.tables")
        # Z_EXTENDED is the merged Z + Z_CIE: never count them twice
        # (legacy skipped the separate Z block likewise)
        if flags["cool_Zext"] and (flags["cool_Z"] or flags["cool_ZCIE"]):
            self.warnings.append("Z_EXTENDED already contains Z and Z_CIE;"
                                 " the separate terms are dropped")
            flags["cool_Z"] = flags["cool_ZCIE"] = False
        if (flags["heat_chem"] or flags["cool_chem"]) and \
                (flags["cool_dH"] or flags["heat_dH"]):
            self.warnings.append("CHEM and DH overlap on the H2 channels "
                                 "(legacy alternatives, not additive)")
        self.configure_floor(cooling_floor, cool, flags)
        return flags

    def configure_floor(self, spec, cool, flags):
        """thermal.floor: no terms (False/''), the default set (True or
        'DEFAULT': every configured radiative line cooler), or an explicit
        list of terms (legacy -useIndividualFloor)."""
        self.floor_terms, self.mol_floor = [], []
        if not spec:
            return
        floorable = [x for x in cool if x in self.FLOOR_NAMES
                     and not (x in ("Z", "Z_CIE") and "Z_EXTENDED" in cool)]
        files = bool(self.file_active)
        if spec is True or (isinstance(spec, str)
                            and spec.strip().upper() == "DEFAULT"):
            terms = floorable
            if files:  # the coolant-file coolants are line coolers too
                flags["floor_file"] = True
            if not terms and not files:
                raise ValueError("thermal.floor: none of the configured "
                                 "coolings can be floored (floorable: %s)"
                                 % ", ".join(sorted(self.FLOOR_NAMES)))
        else:
            if isinstance(spec, str):
                spec = spec.split(",")
            terms = [x.strip().upper() for x in spec if x.strip()]
            if "FILE" in terms:
                if not files:
                    raise ValueError("thermal.floor: FILE needs coolant-"
                                     "file coolants (thermal.coolant_files)")
                flags["floor_file"] = True
                terms = [x for x in terms if x != "FILE"]
            for x in terms:
                if x == "CHEM":
                    raise ValueError(
                        "thermal.floor: CHEM has no cooling term to floor "
                        "(the chemical cooling is part of the net "
                        "heating_chem term; legacy coolingChem is "
                        "identically 0, so its CHEM floor was a no-op)")
                if x not in self.FLOOR_NAMES:
                    raise ValueError("thermal.floor: unknown term %s "
                                     "(floorable: %s)" % (x, ", ".join(
                                         sorted(self.FLOOR_NAMES))))
                if x not in floorable:
                    raise ValueError(
                        "thermal.floor: %s is not a configured cooling%s"
                        % (x, " (Z_EXTENDED contains it: floor Z_EXTENDED)"
                           if x in ("Z", "Z_CIE") and "Z_EXTENDED" in cool
                           else " (add it to thermal.cooling)"))
        if self.cmb_pumping:
            # the CMB pumps the metal-line levels: the floor for those lines
            # (Z, the metal part of Z_EXTENDED) would count it twice
            dropped = [x for x in terms if x in ("Z", "Z_EXTENDED")]
            if dropped:
                self.warnings.append("thermal.cmb_pumping: cooling floor dropped "
                                     "for %s (the CMB is in the level "
                                     "balance)" % ", ".join(dropped))
            terms = [x for x in terms if x not in ("Z", "Z_EXTENDED")]
            if flags.get("floor_file"):
                self.warnings.append("thermal.cmb_pumping: cooling floor dropped "
                                     "for the coolant-file coolants (the "
                                     "CMB is in the level balance)")
                flags["floor_file"] = False
            if not terms:
                self.warnings.append("thermal.floor: nothing left to floor"
                                     " with thermal.cmb_pumping (floor off)")
                return
        for x in terms:
            flags[self.FLOOR_NAMES[x]] = True
            if x in self.MOL_COOLANTS and x not in self.mol_floor:
                self.mol_floor.append(x)
        self.floor_terms = terms + (["FILE"] if flags.get("floor_file")
                                    else [])
        if not self.floor_terms:
            return
        flags["cooling_floor"] = True

    def emit_thermal_config(self, w):
        w("// the thermal processes of this build (codegen thermal.cooling /")
        w("// thermal.heating / thermal.floor); only those set true are")
        w("// compiled in")
        w("struct ThermalConfig {")
        for f in self.THERMAL_FLAGS:
            w("  static constexpr bool %s = %s;"
              % (f, "true" if self.thermal[f] else "false"))
        th = self.thermal
        if self.floor_terms:
            w("  // cooling floor terms: %s%s" % (
                ", ".join(self.floor_terms),
                " (molecules: %s)" % ", ".join(self.mol_floor)
                if self.mol_floor else ""))
        needs_cool = any(th[f] for f in (
            "cool_atomic", "cool_H2", "cool_HD", "cool_compton", "cool_CIE",
            "cool_expansion", "cool_cont"))
        w("  // runtime modules the configured processes need")
        w("  static constexpr bool needs_cooling_hpp = %s;"
          % ("true" if needs_cool else "false"))
        w("  static constexpr bool needs_dust_rates_hpp = %s;"
          % ("true" if th["cool_dust"] else "false"))
        w("  static constexpr bool needs_dust_bins_hpp = %s;"
          % ("true" if th["cool_dust_bins"] else "false"))
        w("  static constexpr bool needs_photo_hpp = %s;"
          % ("true" if th["heat_photo"] else "false"))
        w("};")
        w("")
        w("// what is evolved alongside the chemistry, FIXED at generation")
        w("// (codegen radiation.evolve / dust.evolve_Td). Without rad the photon")
        w("// field is still used (photo rates, photoheating, grains), held")
        w("// at the values the caller passes; with rad it is co-evolved as")
        w("// ln Phi and written back. Thermal evolution as a whole is the")
        w("// one runtime switch (SolverOpts::evolve_T); with it off, Td")
        w("// stays fixed too.")
        w("struct EvolveConfig {")
        w("  static constexpr bool rad = %s;"
          % ("true" if self.evolve_rad else "false"))
        w("  static constexpr bool Td = %s;"
          % ("true" if self.evolve_Td else "false"))
        w("  // dust.evolve: the bins' grain radii as state, and the")
        w("  // processes changing them")
        for nm, key in (("dust", None), ("dust_growth", "GROWTH"),
                        ("dust_sputter", "SPUTTER"), ("dust_evap", "EVAP"),
                        ("dust_kelvin", "KELVIN")):
            on = bool(self.evolve_dust) if key is None \
                else key in self.evolve_dust
            w("  static constexpr bool %s = %s;"
              % (nm, "true" if on else "false"))
        w("};")
        w("")

    # -- naming --------------------------------------------------------------

    def load_species_tols(self, path, what):
        """a legacy per-species tolerance file (-customATOL / -customRTOL):
        'SPECIES VALUE' per line (Fortran reals, '#' comments), the species
        by name or by legacy index name (idx_Hj). Legacy skipped species
        not in the network silently; KROME-X warns. atol must be positive,
        rtol non-negative"""
        if not path:
            return {}
        path = resolve_data_path(path, self.data_dir)
        byidx = {"idx_" + self.cid(s).lower(): s for s in self.net.species}
        out = {}
        for k, raw in enumerate(Path(path).read_text().splitlines(), 1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            f = s.split()
            if len(f) < 2:
                raise ValueError("%s:%d: expected 'SPECIES %s', got %r"
                                 % (path, k, what.upper(), s))
            try:
                v = float(_fortran_real(f[1]))
            except ValueError:
                raise ValueError("%s:%d: %r is not a number"
                                 % (path, k, f[1]))
            if not (v > 0.0 if what == "atol" else v >= 0.0):
                raise ValueError("%s:%d: %s %s must be %s" % (
                    path, k, what, f[1],
                    "positive" if what == "atol" else "non-negative"))
            sp = f[0] if f[0] in self.net.species else \
                byidx.get(f[0].lower())
            if sp is None:
                self.warnings.append("solver.%s_file: species %s not in "
                                     "the network, ignored" % (what, f[0]))
                continue
            out[sp] = v  # a repeated species: the last line wins
        return out

    H2_FORMS = ("H2", "H2_ortho", "H2_para")

    def h2_forms(self):
        """the network's H2 species: H2 and / or its nuclear-spin states"""
        return [x for x in self.H2_FORMS if x in self.net.species]

    def gamma_dens(self, s):
        """the density of a species in the adiabatic index: H2 is the total
        of its forms (nH2_) when the network has ortho / para species"""
        if s == "H2" and self.h2_forms() != ["H2"]:
            return "nH2_(n)"
        return "n[IDX_%s]" % self.cid(s)

    def cid(self, species):
        """Species name -> C identifier fragment (legacy-style j/k charges)."""
        s = species.replace("+", "j").replace("-", "k")
        if s == "ek":
            s = "E"
        if s[0].isdigit():
            s = "i" + s
        return re.sub(r"\W", "_", s)

    # -- gas+total representation helpers --------------------------------------

    def nden(self, sp):
        """C expression for a species' physical number density; in
        gas+total mode a surface slot stores gas+surface, so the surface
        density is the difference of the two slots."""
        g = self.pair_gas.get(sp)
        if g is not None:
            return "(n[IDX_%s] - n[IDX_%s])" % (self.cid(sp), self.cid(g))
        return "n[IDX_%s]" % self.cid(sp)

    def _tdust_define(self):
        """the Tdust state symbol always routes through table_Tdust:
        p.Tdust by default (the stub), the live table when dust.tables
        is generated (then Tdust-bearing rates are n-dependent and go
        through the numerical-dk Jacobian path), and the binned-dust
        area-weighted mean whenever a DustBins hook is installed."""
        return "#define Tdust table_Tdust(Tgas, n, p)"

    def den_cols(self, sp):
        """Jacobian columns a species' density depends on: [(slot, sign)]."""
        g = self.pair_gas.get(sp)
        if g is not None:
            return [(sp, 1), (g, -1)]
        return [(sp, 1)]

    def stoich_net(self, rx):
        """net state-slot increments of one reaction as [(species, coeff)];
        in gas+total mode a gas member also increments its pair's total
        slot, and exact cancellations (freezeout/evaporation on the total)
        are dropped here rather than emitted as -=/+= pairs."""
        net_c = {}

        def add(sp, v):
            net_c[sp] = net_c.get(sp, 0) + v
            if sp in self.total_of:
                t = self.total_of[sp]
                net_c[t] = net_c.get(t, 0) + v

        for r in rx["reactants"]:
            add(r["sp"], -1)
        for pr in rx["products"]:
            add(pr["sp"], +1)
        return [(s, v) for s, v in sorted(net_c.items()) if v != 0]

    # -- n-dependence --------------------------------------------------------

    def ndep(self, ast):
        key = id(ast)
        if key in self._ndep_cache:
            return self._ndep_cache[key]
        dep = False
        for nd in knet.walk(ast):
            if nd.get("kind") == "density" or "sp" in nd \
                    or nd.get("s") in ("ntot", "Hnuclei"):
                dep = True
            elif nd.get("s") == "Tdust" and self.dust_tables:
                dep = True  # live Tdust = table_Tdust(Tgas, n, p)
            elif nd.get("kind") == "var":
                dep = dep or self.ndep(self.net.vars[nd["s"]])
        self._ndep_cache[key] = dep
        return dep

    # -- differentiation d/d n(target) ---------------------------------------

    def diff(self, a, target):
        """the derivative of the expression AST a (knet.py's node forms)
        with respect to target, as a new AST (ZERO / ONE for constants):
        target is a species name (d/dn), "Tgas", or "@ntot" / "@Hnuclei"
        for the aggregates (whose own dependence on the densities is the
        Jacobian's low-rank part)"""
        if "n" in a or "str" in a:
            return ZERO
        if "s" in a:
            if a.get("kind") == "var":
                return self.diff(self.net.vars[a["s"]], target)
            # the aggregates ntot / Hnuclei: a derivative with respect to
            # the aggregate itself (target "@ntot" / "@Hnuclei"); their
            # dependence on the densities is the low-rank part U W^T of
            # the Jacobian (jacobian_split), not per-species terms
            if a["s"] in ("ntot", "Hnuclei"):
                return ONE if target == "@" + a["s"] else ZERO
            return ZERO  # other state/param/const are n-independent
        if "sp" in a:
            return ZERO
        if "u" in a:
            return {"u": "-", "x": self.diff(a["x"], target)}
        if "c" in a:
            return ZERO
        if "b" in a:
            l, r = a["l"], a["r"]
            dl, dr = self.diff(l, target), self.diff(r, target)
            if a["b"] == "+":
                return add(dl, dr)
            if a["b"] == "-":
                return sub(dl, dr)
            if a["b"] == "*":
                return add(mul(dl, r), mul(l, dr))
            if a["b"] == "/":
                return sub(div(dl, r), div(mul(l, dr), mul(r, r)))
            if a["b"] == "**":
                if is_num(dl, 0.0) and is_num(dr, 0.0):
                    return ZERO
                # d(l^r) through the underflow-safe runtime helper dpow_:
                # the naive l^r*(dr*log(l)+r*dl/l) is NaN once l underflows
                # to 0 (log(0)*0), e.g. exp(-E/T) factors of density-
                # interpolated rates (Shapiro-Kang H2 dissociation) at low
                # Tgas, which poisons the whole analytic Jacobian
                return {"f": "dpow_", "kind": "builtin",
                        "a": [l, r, dl, dr]}
        if "f" in a:
            f, args = a["f"], a["a"]
            if a.get("kind") == "density":
                # SLOT derivative: in gas+total mode a surface species'
                # density (total - gas) also depends on the gas slot
                sp = args[0].get("sp")
                if target.startswith("@"):
                    return ZERO
                if sp == target:
                    return ONE
                if self.pair_gas.get(sp) == target:
                    return {"n": -1.0}
                return ZERO
            if a.get("kind") == "tfunc":  # functions of Tgas and constants
                return ZERO
            if a.get("kind") in ("extern", "table"):
                if not self.state_free_extern(a):
                    self.warnings.append("d/dn of %s() treated as 0" % f)
                return ZERO
            da = [self.diff(x, target) for x in args]
            x = args[0]
            if f == "exp":
                return mul(a, da[0])
            if f == "log":
                return div(da[0], x)
            if f == "log10":
                return div(da[0], mul(x, {"n": math.log(10.0)}))
            if f == "sqrt":
                return div(da[0], mul({"n": 2.0}, a))
            if f == "tanh":
                return mul(sub(ONE, mul(a, a)), da[0])
            if f == "abs":
                return {"f": "where", "kind": "builtin",
                        "a": [{"c": ">=", "l": x, "r": ZERO},
                              da[0], {"u": "-", "x": da[0]}]}
            if f in ("min", "max"):
                op = "<" if f == "min" else ">"
                return {"f": "where", "kind": "builtin",
                        "a": [{"c": op, "l": args[0], "r": args[1]},
                              da[0], da[1]]}
            if f == "where":
                return {"f": "where", "kind": "builtin",
                        "a": [args[0], da[1], da[2]]}
        raise ValueError("cannot differentiate node %r" % a)

    # -- C++ expression printer ----------------------------------------------

    def cnum(self, v):
        if v == math.inf:
            return "INFINITY"
        s = repr(float(v))
        return s if ("e" in s or "." in s or "inf" in s) else s + ".0"

    # precedence levels for minimal parenthesization:
    #  0 comparisons, 1 add/sub, 2 mul/div/unary-minus, 9 atoms and calls
    def cexpr(self, a, tvar="Tgas"):
        return self._cx(a, tvar)[0]

    def _wrap(self, child, tvar, need):
        s, p = self._cx(child, tvar)
        return "(%s)" % s if p < need else s

    def _cx(self, a, tvar):
        """the C++ text of the expression AST a, and its precedence (1: a
        sum, 2: a product or a unary minus, 9: an atom or a call), from
        which _wrap decides the parentheses; tvar: the C++ name of the
        rate's (clamped) temperature"""
        if "n" in a:
            s = self.cnum(a["n"])
            return s, (2 if s.startswith("-") else 9)
        if "str" in a:
            return '"%s"' % a["str"], 9
        if "s" in a:
            k, s = a.get("kind"), a["s"]
            if k == "state":
                if a.get("unclamped"):  # the actual gas T (1/K_c)
                    return ("Tg_" if self._in_ratefn else "Tgas"), 9
                if s == "Tgas":
                    return tvar, 9
                if s == "Tgas_eV":
                    return "%s*constants::boltzmann_eV" % tvar, 2
                if s == "ntot":
                    return "ntot_(n)", 9
                if s == "Hnuclei":
                    return "Hnuclei_(n)", 9
                if s == "J_LW":  # band-mean J of the photon bins (5.1)
                    return ("1.0" if self._jlw_unit else "photo_J_LW_(p)"), 9
                return s, 9  # Tdust: runtime-provided, Params for now
            if k == "ncol":  # column argument of a column-form shield
                return "Ncol", 9
            if k == "param":
                return "p." + s, 9
            if k == "var":
                if s in self._hoisted:
                    return self._hoisted[s], 9
                return "v_%s(%s, n, p)" % (s, tvar), 9
            if k == "const":
                return self.cnum(knet.CONSTANTS[s]), 9
        if "sp" in a:
            return "IDX_%s" % self.cid(a["sp"]), 9
        if "sps" in a:  # species-list argument (only reaches zero stubs)
            return "0 /* [%s] */" % ", ".join(a["sps"]), 9
        if "u" in a:
            return "-%s" % self._wrap(a["x"], tvar, 3), 2
        if "c" in a:
            return "%s %s %s" % (self._wrap(a["l"], tvar, 1), a["c"],
                                 self._wrap(a["r"], tvar, 1)), 0
        if "b" in a:
            op = a["b"]
            if op == "**":
                r = a["r"]
                if is_num(r) and float(r["n"]).is_integer() \
                        and 2 <= abs(r["n"]) <= 4:
                    base = self._wrap(a["l"], tvar, 9)
                    m = "*".join([base] * int(abs(r["n"])))
                    return (m, 2) if r["n"] > 0 else ("1.0/(%s)" % m, 2)
                # other integer and half-integer literal exponents without
                # std::pow (roundoff-level differences, user-approved):
                # compile-time products, sqrt for the halves
                e = self.lit_num(r)
                if self.cheap_power(e):
                    x = self.cexpr(a["l"], tvar)
                    if e.is_integer():
                        return "ipow_<%d>(%s)" % (int(e), x), 9
                    if (2.0 * e).is_integer() and abs(e) <= 2.5:
                        k = int(abs(e) - 0.5)  # |e| = k + 1/2
                        sq = "std::sqrt(%s)" % x
                        m = sq if k == 0 else "ipow_<%d>(%s) * %s" % (k, x, sq)
                        return ((m, 9 if k == 0 else 2) if e > 0
                                else ("1.0/(%s)" % m, 2))
                return "std::pow(%s, %s)" % (self.cexpr(a["l"], tvar),
                                             self.cexpr(a["r"], tvar)), 9
            if op in ("+", "-"):
                return "%s %s %s" % (self._wrap(a["l"], tvar, 1), op,
                                     self._wrap(a["r"], tvar, 2)), 1
            # * and /: the right operand of / must bind tighter
            return "%s %s %s" % (self._wrap(a["l"], tvar, 2), op,
                                 self._wrap(a["r"], tvar,
                                            3 if op == "/" else 2)), 2
        if "f" in a:
            f, args = a["f"], a["a"]
            if a.get("kind") == "density":
                return self.nden(args[0]["sp"]), 9
            if f == "where":
                return "(%s ? %s : %s)" % tuple(self.cexpr(x, tvar)
                                                for x in args), 9
            if f in ("min", "max"):
                return "std::f%s(%s, %s)" % (f, self.cexpr(args[0], tvar),
                                             self.cexpr(args[1], tvar)), 9
            if f in ("exp", "log", "log10", "sqrt", "tanh", "abs"):
                cf = "std::fabs" if f == "abs" else "std::" + f
                return "%s(%s)" % (cf, self.cexpr(args[0], tvar)), 9
            if a.get("kind") == "tfunc":  # database helper, pure in T
                return "%s(%s)" % (f, ", ".join(self.cexpr(x, tvar)
                                                for x in args)), 9
            if f == "dpow_":  # underflow-safe d(l^r), core helper
                return "dpow_(%s)" % ", ".join(self.cexpr(x, tvar)
                                               for x in args), 9
            if a.get("kind") in ("extern", "table"):
                s = ", ".join(self.cexpr(x, tvar) for x in args)
                if f in KNOWN_EXTERNS or GENERATED_EXTERNS.get(f) \
                        or f in self._gen_ext:
                    # append the implicit state args
                    s = ", ".join(x for x in (s, "%s, n, p" % tvar) if x)
                if f == "num2col":
                    # legacy num2col(x): the network's own num2col_ (the
                    # --num2col method), not a second definition
                    f = "num2col_"
                return "%s(%s)" % (f, s), 9
        raise ValueError("cannot print node %r" % a)

    # -- emission ------------------------------------------------------------

    def _num2col_c(self, x):
        """C statements computing N = num2col(x) into `N` (the legacy
        -columnDensityMethod choice, radiation.num2col)."""
        if self.num2col == "default":
            return ["const double N = 1.87e21 * std::pow(std::fmax(%s, "
                    "1e-40) * 1e-3, 2.0/3.0);" % x]
        if self.num2col == "length":
            return ["const double N = std::fmax(%s, 1e-40) * p.Lcol;" % x]
        T = ("std::fmin(std::fmax(Tgas, p.Tcmb), 4e1)"
             if self.num2col == "jeans40" else "std::fmax(Tgas, p.Tcmb)")
        return ["const double rr = std::fmax(rho_(n), 1e-40);",
                "const double mu = rr / std::fmax(ntot_(n), 1e-40) "
                "/ constants::p_mass;",
                "const double lj = std::sqrt(constants::pi * constants::boltzmann_erg "
                "* %s / rr / constants::p_mass / constants::gravity / mu);" % T,
                "const double N = 0.5 * std::fmax(%s, 1e-40) * lj;" % x]

    def emit_generated_externs(self, w):
        """definitions of the legacy helpers KROME-X generates (GENERATED_EXTERNS)
        that this network uses."""
        used = [nm for nm in sorted(self.net.externs)
                if nm in GENERATED_EXTERNS]
        if not used:
            return
        w("// legacy helper functions used by this network (ported from")
        w("// src/krome_phfuncs.f90 / krome_grfuncs.f90)")
        if "shield_dust" in used:
            hs = [(sp, f) for sp, f in (("H", 1.0), ("H+", 1.0))
                  if sp in self.net.species]
            w("// legacy shield_dust(n, Tgas, gam) (krome_phfuncs): dust")
            w("// shielding exp(-sigma_d N_H), sigma_d = 2e-21 (d2g/0.01) gam")
            w("// cm^2 (Richings et al. 2014; legacy dust2gas_ratio = d2g /")
            w("// 0.01), N_H = N(H) + N(H+) + 2 N(H2) from num2col_ of each")
            w("double inline shield_dust(double gam, double Tgas, "
              "const double* n, const Params& p) {")
            w("  const double sigma_d = 2e-21 * (p.d2g / 1e-2) * gam;")
            w("  double NH = 0.0;")
            for sp, f in hs:
                w("  NH += %s * num2col_(n[%s], Tgas, n, p);"
                  % (self.cnum(f), "IDX_%s" % self.cid(sp)))
            if self.h2_forms():
                w("  NH += 2.0 * num2col_(nH2_(n), Tgas, n, p);")
            w("  return std::exp(-sigma_d * NH);")
            w("}")
        if "stars_screen" in used:
            zs = []
            for sp in sorted(self.net.species):
                z = self.nucleus_Z(sp)
                zs.append(z if z is not None else 0)
            xh = ("n[IDX_H] * species_mass_g[IDX_H] / rho"
                  if "H" in self.net.species else "0.0")
            w("// nuclear charge of each species (0: electrons, molecules)")
            w("inline const int species_Z[NSPEC] = {%s};"
              % ", ".join(map(str, zs)))
            w("// plasma screening factor of a Z1 + Z2 reaction (legacy")
            w("// krome_stars stars_screen): Itoh et al. 1977/79 strong")
            w("// screening, else Graboske et al. 1973 / DeWitt et al. 1973")
            w("// weak, intermediate and strong. Y_i = n_i / (rho N_A)")
            w("// [mol/g], rho from the state. Legacy's power tables stopped")
            w("// at Z = 26 (and read out of bounds at Z = 0); X_H is the H")
            w("// mass fraction (legacy passed the molar abundance)")
            w("static inline double stars_screen(double z1, double z2, "
              "double Tgas, const double* n, const Params& p) {")
            w("  (void)p;")
            w("  const double rho = std::fmax(rho_(n), 1e-300);")
            w("  const double NA = %s;"
              % self.cnum(knet.CONSTANTS["N_avogadro"]))
            w("  double ZY = 0.0, sY = 0.0, sZZt = 0.0, sZ158 = 0.0;")
            w("  double Y[NSPEC];")
            w("  for (int i = 0; i < NSPEC; ++i) {")
            w("    Y[i] = std::fmax(n[i], 0.0) / (rho * NA);")
            w("    ZY += species_Z[i] * Y[i];")
            w("    sY += Y[i];")
            w("  }")
            w("  if (!(ZY > 0.0) || !(sY > 0.0)) return 1.0;")
            w("  // Itoh: strong screening")
            w("  const double Atmp = 1.442250 / std::cbrt(7.567682e24 * rho "
              "* ZY);")
            w("  const double T13 = std::cbrt(Tgas);")
            w("  {")
            w("    const double a = std::cbrt(z1) * Atmp;")
            w("    const double b = std::cbrt(z2) * Atmp;")
            w("    const double c = z1 * z2 * 2.307121e-19 / (0.5 * (a + b) "
              "* constants::boltzmann_erg * Tgas);")
            w("    const double d = z1 * z2 * 5351.69 / (std::cbrt(z1 + z2) "
              "* T13);")
            w("    const double chk = 3.0 * c / d;")
            w("    if (chk >= 0.2)")
            w("      return std::exp(1.25 * c - 0.095 * d * chk * chk);")
            w("  }")
            w("  // Graboske & DeWitt")
            w("  const double XH = %s;" % xh)
            w("  const double fer = 5.4885e7 * rho * (1.0 + XH) / (Tgas * "
              "std::sqrt(Tgas));")
            w("  const double theta = fer <= 3.0")
            w("      ? 1.0 / (1.0 + 0.39716 * fer - 0.00929 * fer * fer)")
            w("      : 1.1447 / (std::pow(fer, 2.0 / 3.0) * (1.0 + 0.28077 "
              "/ fer));")
            w("  for (int i = 0; i < NSPEC; ++i) {")
            w("    sZZt += species_Z[i] * Y[i] * (species_Z[i] + theta);")
            w("    sZ158 += Y[i] * std::pow(species_Z[i], 1.58);")
            w("  }")
            w("  const double pmol = 1.0 / sY;")
            w("  const double Zm = ZY * pmol, Zt = std::sqrt(sZZt * pmol);")
            w("  const double Z158m = sZ158 * pmol;")
            w("  const double lamb0 = 1.88e8 * std::sqrt(rho / (Tgas * Tgas "
              "* Tgas * pmol));")
            w("  const double lamb12 = lamb0 * z1 * z2 * Zt;")
            w("  if (lamb12 <= 0.1) return std::exp(lamb12);  // weak")
            w("  const double l0s = std::pow(lamb0, 2.0 / 3.0);")
            w("  const double Zm13 = std::cbrt(Zm);")
            w("  const double a5 = std::pow(z1 + z2, 5.0 / 3.0) - std::pow(z1,"
              " 5.0 / 3.0) - std::pow(z2, 5.0 / 3.0);")
            w("  const double b4 = std::pow(z1 + z2, 4.0 / 3.0) - std::pow(z1,"
              " 4.0 / 3.0) - std::pow(z2, 4.0 / 3.0);")
            w("  const double c2 = (std::pow(z1 + z2, 2.0 / 3.0) - std::pow("
              "z1, 2.0 / 3.0) - std::pow(z2, 2.0 / 3.0)) / l0s;")
            w("  const double H12s = std::fmin(0.624 * Zm13 * (a5 + 0.316 * b4"
              " * Zm13 + 0.737 * c2 / Zm13) * l0s, 200.0);")
            w("  if (lamb12 <= 5.0) {  // intermediate")
            w("    const double a = std::pow(z1 + z2, 1.86) - std::pow(z1, "
              "1.86) - std::pow(z2, 1.86);")
            w("    const double H12i = 0.38 * a * std::pow(lamb0, 0.86) * "
              "Z158m / (std::pow(Zt, 0.58) * std::pow(Zm, 0.28));")
            w("    if (lamb12 <= 2.0) return std::exp(H12i);")
            w("    return std::exp(std::fmin(H12i, H12s));")
            w("  }")
            w("  return std::exp(H12s);  // strong")
            w("}")
        if "fHnOj" in used:
            w("static inline double fHnOj(double Av) {")
            w("  // CO photodissociation attenuation vs visual extinction")
            w("  return Av <= 15.0 ? std::exp(-2.55 * Av + 0.0165 * Av * Av)")
            w("                    : std::exp(-2.8 * Av);")
            w("}")
        if "krate_H2dust_table" in used:
            w("// H2 formation on dust from the dust tables (defined with")
            w("// the tables below)")
            w("inline double krate_H2dust_table(double Tgas, const double* n,"
              " const Params& p);")
        if "n_global" in used:
            w("KROME_HD static inline double n_global(int sp, double, "
              "const double* n,")
            w("                                       const Params& p) {")
            w("  // abundance at the START of the solver call (legacy")
            w("  // n_global), e.g. a fallback when n(sp) turns negative")
            w("  return p.hk.n_global ? p.hk.n_global[sp] : n[sp];")
            w("}")
        if "get_table_Tdust" in used:
            w("static inline double get_table_Tdust(double Tgas, "
              "const double* n,")
            w("                                     const Params& p) {")
            w("  return table_Tdust(Tgas, n, p);  // live Tdust")
            w("}")
        if "get_ratioFluxDraine" in used:
            w("KROME_HD static inline double get_ratioFluxDraine(double, "
              "const double*,")
            w("                                                  "
              "const Params& p) {")
            w("  // photobin field over the Draine ISRF (legacy")
            w("  // get_integratedFlux()/7.95e-8): sum_b J_b dE_b / E_b =")
            w("  // sum_b Phi_b h / (4 pi), h = planck_eV; 0 without a field")
            w("  double s = 0.0;")
            w("  for (int b = 0; b < p.hk.nbins; ++b) s += p.hk.phi[b];")
            w("  return s * constants::planck_eV / (4.0 * constants::pi) "
              "/ 7.95e-8;")
            w("}")
        if "krome_fshield" in used:
            sh = self.shielding
            w("static inline double krome_fshield(double Tgas, "
              "const double* n,")
            w("                                   const Params& p) {")
            w("  (void)Tgas; (void)n; (void)p;")
            if not sh:
                w("  return 1.0;  // radiation.shielding = 'NONE'")
            else:
                # the H2 column: H2, or its nuclear-spin species (networks
                # with H2_ortho / H2_para: the total H2 shields)
                h2 = [x for x in ("H2", "H2_ortho", "H2_para")
                      if x in self.net.species]
                if not h2:
                    raise ValueError("radiation.shielding needs H2 (or H2_ortho / "
                                     "H2_para) in the network")
                w("  // legacy -shielding=%s (calc_H2shield%s), column from"
                  % (sh, sh))
                w("  // num2col (radiation.num2col=%s) of %s" % (self.num2col,
                                                         " + ".join(h2)))
                w("  const double nH2 = std::fmax(1e-40, nH2_(n));")
                for ln in self._num2col_c("nH2"):
                    w("  " + ln)
                w("  const double N_H2 = 2.0 * N;")
                if sh == "DB96":
                    w("  return std::fmin(1.0, std::pow(N_H2 * 1e-14, "
                      "-0.75));")
                else:
                    w("  const double H_mass = constants::p_mass + "
                      "constants::e_mass;")
                    if sh == "WG11":
                        w("  const double xN = N_H2 * 2e-15;")
                        w("  const double b5 = std::sqrt(constants::boltzmann_erg * Tgas "
                          "/ H_mass) * 1e-5;")
                        w("  return 0.965 / std::pow(1.0 + xN / b5, 1.1)")
                        w("      + 0.035 / std::sqrt(1.0 + xN)")
                        w("      * std::exp(-8.5e-4 * std::sqrt(1.0 + xN));")
                    else:  # R14, with the corrected piecewise (the legacy
                        # .or. made its T >= 4e3 K branch dead code)
                        w("  const double bturb = 7.1e5;  // 7.1 km/s")
                        w("  const double b5 = std::sqrt(constants::boltzmann_erg * Tgas "
                          "/ H_mass + bturb * bturb) * 1e-5;")
                        w("  const double om = 0.013 * std::pow(1.0 + "
                          "std::pow(Tgas / 2.7e3, 1.3), 1.0 / 1.3)")
                        w("      * std::exp(-std::pow(Tgas / 3.9e3, 14.6));")
                        w("  double alpha, Ncrit;")
                        w("  if (Tgas < 3e3) {")
                        w("    alpha = 1.4;")
                        w("    Ncrit = 1.3 * (1.0 + std::pow(Tgas / 6e2, "
                          "0.8));")
                        w("  } else if (Tgas < 4e3) {")
                        w("    alpha = std::pow(Tgas / 4.5e3, -0.8);")
                        w("    Ncrit = std::pow(Tgas / 4.76e3, -3.8);")
                        w("  } else {")
                        w("    alpha = 1.1;")
                        w("    Ncrit = 2.0;")
                        w("  }")
                        w("  const double xN = N_H2 * 1e-14 / Ncrit;")
                        w("  return (1.0 - om) / std::pow(1.0 + xN / b5, "
                          "alpha)")
                        w("      * std::exp(-5e-7 * (1.0 + xN))")
                        w("      + om / std::sqrt(1.0 + xN)")
                        w("      * std::exp(-8.5e-4 * std::sqrt(1.0 + xN));")
            w("}")
        w("")

    def _auto_ast(self, rate):
        """legacy database rate string -> classified KROME-X AST, or (None,
        reason) when it needs machinery KROME-X does not have."""
        import convert_legacy as cl
        s = cl.normalize_expr(rate)
        s = re.sub(r"(?<![A-Za-z0-9_])tgas(?![A-Za-z0-9_])", "Tgas", s,
                   flags=re.I)
        # legacy shortcut symbols, case-insensitive like legacy
        for name in sorted(cl.SHORTCUTS, key=len, reverse=True):
            s = re.sub(r"(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_(])"
                       % re.escape(name), "(%s)" % cl.SHORTCUTS[name], s,
                       flags=re.I)
        try:
            ast = knet.parse_expr(s, 0)
        except knet.ParseError as exc:
            return None, "unparsable rate (%s)" % exc
        used = set()
        for nd in knet.walk(ast):
            if "s" in nd and "kind" not in nd:
                if nd["s"] in ("Tgas", "Tgas_eV"):
                    nd["kind"] = "state"
                elif nd["s"] in knet.CONSTANTS:
                    nd["kind"] = "const"
                else:
                    return None, "symbol %s" % nd["s"]
            elif "f" in nd and "kind" not in nd:
                if nd["f"] in AUTO_TFUNCS:
                    nd["kind"] = "tfunc"
                    used.add(nd["f"])
                elif nd["f"] in knet.BUILTIN_FUNCS:
                    nd["kind"] = "builtin"
                else:
                    return None, "function %s()" % nd["f"]
        self._tfuncs |= used
        return ast, None

    def resolve_auto(self):
        """bare `auto` rates (network-format.md 7): look the reaction up in
        the legacy rate database (data/database/*.dat, in load_auto_db
        order) with legacy's rules: the first matching entry without T
        limits wins; entries with limits are range pieces of one piecewise
        rate; an unlimited entry after limited ones is an overlap error.
        A reaction without an entry, or whose entry needs machinery
        KROME-X lacks (reverse reactions, legacy dust variables), is an
        error naming the cause (surface reactions never get here: knet
        makes their bare auto surface(auto))."""
        if getattr(self, "_auto_done", False):
            return
        self._auto_done = True
        self._tfuncs = set()
        net = self.net
        if not any(rx["rate"]["kind"] == "auto" for rx in net.reactions):
            return
        db = load_auto_db(self.data_dir)
        out = []
        for rx in net.reactions:
            if rx["rate"]["kind"] != "auto":
                out.append(rx)
                continue
            R = sorted(r["sp"] for r in rx["reactants"])
            P = sorted(p0["sp"] for p0 in rx["products"])
            pieces, problem = [], None
            for e in db:
                if e["reacts"] != R or e["prods"] != P:
                    continue
                if e["limits"] is None:
                    if pieces:
                        problem = ("database has limited and unlimited "
                                   "entries (overlap)")
                    else:
                        pieces = [e]
                    break
                pieces.append(e)
            if not pieces and problem is None:
                problem = "no database entry"
            asts = []
            for e in pieces if problem is None else ():
                ast, why = self._auto_ast(e["rate"])
                if ast is None:
                    problem = "%s in %s" % (why, e["fname"])
                    break
                asts.append((e, ast))
            if problem is not None:
                # an error, as in legacy: a rate of 0 would silently drop
                # the reaction
                raise ValueError(
                    "auto rate %s (line %s): %s. Give the reaction an "
                    "explicit rate, remove it, or add it to "
                    "data/database" % (rx["eq"], rx["line"], problem))
            if rx["trange"] and any(e["limits"] for e, _ in asts):
                self.warnings.append("auto rate %s: database T limits "
                                     "replace the network's range"
                                     % rx["eq"])
            for e, ast in asts:
                nrx = dict(rx)
                nrx["rate"] = {"kind": "expr", "ast": ast,
                               "auto_src": e["fname"]}
                nrx["trange"] = e["limits"] or rx["trange"]
                out.append(nrx)
        net.reactions = out
        for i, rx in enumerate(net.reactions):
            rx["id"] = i + 1
            rx.pop("group", None)
        knet.group_reactions(net)

    # -- kinetic nucleation theory ---------------------------------------------

    def resolve_nucleation(self):
        """cluster growth / destruction (legacy krome_subs, Gail &
        Sedlmayr 2013 ch. 13):
          cluster_growth_rate(M, N[, stick]) and
          general_cluster_growth_rate(M, N1, N2[, stick]): generated
            (v_thermal x geometric cross section x stick, monomer radius
            from NUCLEATION_RADII, cluster radius r1 N^(1/3));
          revKc_with_GFE([R...], [P...]): K_c of R -> P from the tabulated
            Gibbs free energies (data/thermochemistry/*.gfe), a generated
            function of the actual gas temperature;
          cluster_destruction_rate(M, N[, stick]) in N-mer -> (N-1)-mer + M:
            detailed balance, cluster_growth_rate(M, N-1[, stick]) times
            revKc_with_GFE of the reaction itself (undefined in legacy:
            react_nucleation never compiled)."""
        if getattr(self, "_nucl_done", False):
            return
        self._nucl_done = True
        net = self.net
        self._gfe = []  # (R list, P list) per revKc_with_GFE site
        used = set()

        def setnode(nd, new):
            nd.clear()
            nd.update(new)

        for rx in net.reactions:
            if rx["rate"]["kind"] != "expr":
                continue
            ast = rx["rate"]["ast"]
            for nd in list(knet.walk(ast)):
                if nd.get("f") == "cluster_destruction_rate":
                    a = nd["a"]
                    if len(a) not in (2, 3):
                        raise ValueError("%s: cluster_destruction_rate(M, N"
                                         "[, stick])" % rx["eq"])
                    grow = {"f": "cluster_growth_rate", "kind": "extern",
                            "a": [a[0], {"b": "-", "l": a[1],
                                         "r": {"n": 1.0}}] + a[2:]}
                    kc = {"f": "revKc_with_GFE", "kind": "extern",
                          "a": [{"sps": [q["sp"] for q in rx["reactants"]
                                         if q["kinetic"]]},
                                {"sps": [q["sp"] for q in rx["products"]]}]}
                    setnode(nd, {"b": "*", "l": grow, "r": kc})
            for nd in list(knet.walk(ast)):
                f = nd.get("f")
                if f == "revKc_with_GFE":
                    a = nd["a"]
                    if len(a) != 2 or "sps" not in a[0] or "sps" not in a[1]:
                        raise ValueError("%s: revKc_with_GFE([reactants], "
                                         "[products])" % rx["eq"])
                    q = len(self._gfe)
                    self._gfe.append((list(a[0]["sps"]), list(a[1]["sps"])))
                    setnode(nd, {"f": "gfeKc_%d" % q, "kind": "tfunc",
                                 "a": [{"s": "Tgas", "kind": "state",
                                        "unclamped": True}]})
                elif f in ("cluster_growth_rate",
                           "general_cluster_growth_rate"):
                    nmin = 2 if f == "cluster_growth_rate" else 3
                    if len(nd["a"]) not in (nmin, nmin + 1) \
                            or "sp" not in nd["a"][0]:
                        raise ValueError("%s: %s(monomer, size%s[, stick])"
                                         % (rx["eq"], f, "" if nmin == 2
                                            else ", size2"))
                    if len(nd["a"]) == nmin:
                        nd["a"].append({"n": 1.0})
                    used.add(nd["a"][0]["sp"])
                    net.externs.add(f)
                elif f == "troe_falloff":
                    net.externs.add(f)
        net.externs.discard("cluster_destruction_rate")
        net.externs.discard("revKc_with_GFE")
        miss = sorted(m for m in used
                      if m.upper() not in NUCLEATION_RADII)
        if miss:
            raise ValueError("cluster growth: no monomer radius for %s "
                             "(known: %s)" % (", ".join(miss),
                                              ", ".join(NUCLEATION_RADII)))
        need = sorted({x for R, P in self._gfe for x in R + P})
        self._gfe_data = {}
        miss = []
        for sp in need:
            path = self.data_dir.parent / "thermochemistry" / (
                sp.upper() + ".gfe")
            if not path.exists():
                miss.append(sp)
                continue
            rows = [ln.split() for ln in path.read_text().splitlines()
                    if ln.strip() and not ln.lstrip().startswith("#")]
            T = [float(r[0]) for r in rows]
            G = [float(r[1]) for r in rows]
            last = max(i for i, g in enumerate(G) if not math.isnan(g))
            if any(math.isnan(g) for g in G[:last]):
                raise ValueError("%s.gfe: NaN inside the table"
                                 % sp.upper())
            if last < len(G) - 1:  # CO2, OH: no data above 3000 K
                self.warnings.append(
                    "%s.gfe: no data above %g K, G held at that value"
                    % (sp.upper(), T[last]))
                G[last + 1:] = [G[last]] * (len(G) - last - 1)
            self._gfe_data[sp] = (T, G)
        if miss:
            raise ValueError("revKc_with_GFE: no data/thermochemistry/"
                             "<NAME>.gfe for %s" % ", ".join(miss))

    def emit_nucleation(self, w):
        net = self.net
        if any(x in net.externs for x in ("cluster_growth_rate",
                                         "general_cluster_growth_rate")):
            w("// monomer radii for cluster growth [cm] (legacy get_radius)")
            w("inline const double nucleation_radius[NSPEC] = {%s};"
              % ", ".join(self.cnum(NUCLEATION_RADII.get(s.upper(), 0.0))
                          for s in sorted(net.species)))
            w("// kinetic nucleation theory (Gail & Sedlmayr 2013 ch. 13):")
            w("// k = v_thermal(reduced mass) x pi (r_a + r_b)^2 x stick, a")
            w("// cluster of N monomers having radius r1 N^(1/3) (legacy")
            w("// cluster_growth_rate / general_cluster_growth_rate)")
            w("static inline double general_cluster_growth_rate(int m, "
              "double N1, double N2, double stick, double Tgas, "
              "const double* n, const Params& p) {")
            w("  (void)n; (void)p;")
            w("  const double im = 1.0 / species_mass_g[m];")
            w("  const double ir = im / N1 + im / N2;")
            w("  const double v = std::sqrt(8.0 * constants::boltzmann_erg * Tgas * ir "
              "/ constants::pi);")
            w("  const double r1 = nucleation_radius[m];")
            w("  const double ra = r1 * std::cbrt(N1), rb = r1 * "
              "std::cbrt(N2);")
            w("  return v * constants::pi * (ra + rb) * (ra + rb) * stick;")
            w("}")
            w("// monomer addition to an N-cluster")
            w("static inline double cluster_growth_rate(int m, double N, "
              "double stick, double Tgas, const double* n, const Params& p)"
              " {")
            w("  return general_cluster_growth_rate(m, 1.0, N, stick, Tgas, "
              "n, p);")
            w("}")
        if "troe_falloff" in net.externs:
            w("// simplified Troe falloff (legacy troe_falloff: N = 1, c = 0)")
            w("static inline double troe_falloff(double k0, double kinf, "
              "double Fc, double m) {")
            w("  // k0 or kinf underflowed (low T): the rate is 0 (legacy")
            w("  // computed 0/0 = NaN)")
            w("  if (!(k0 * m > 0.0) || !(kinf > 0.0)) return 0.0;")
            w("  const double rm = k0 * m / kinf;")
            w("  const double lr = std::log10(rm);")
            w("  return k0 * m / (1.0 + rm) * std::pow(Fc, 1.0 / (1.0 + lr "
              "* lr));")
            w("}")
        gfe = getattr(self, "_gfe", [])
        if not gfe:
            return
        sps = sorted(self._gfe_data)
        T0 = self._gfe_data[sps[0]][0]
        for sp in sps:
            if self._gfe_data[sp][0] != T0:
                raise ValueError("revKc_with_GFE: %s.gfe has a different "
                                 "temperature grid" % sp.upper())
        nT = len(T0)
        vals = [v for sp in sps for v in self._gfe_data[sp][1]]
        w("// Gibbs free energies [kJ/mol] on a uniform T grid (data/")
        w("// thermochemistry/*.gfe; legacy read every other row, and")
        w("// read zeros above ~3040 K): gfe_tab[s * GFE_NT + i]")
        w("constexpr int GFE_NT = %d;" % nT)
        w("constexpr double GFE_T0 = %s, GFE_DT = %s;"
          % (self.cnum(T0[0]), self.cnum(T0[1] - T0[0])))
        self.emit_data_array(w, "gfe_tab", vals,
                             "Gibbs free energies, %d species" % len(sps))
        w("// G [kJ/mol] of table slot s at T: linear, clamped to the grid")
        w("static inline double gfe_(int s, double T) {")
        w("  double x = (T - GFE_T0) / GFE_DT;")
        w("  x = std::fmin(std::fmax(x, 0.0), GFE_NT - 1.0);")
        w("  int i = static_cast<int>(x);")
        w("  if (i > GFE_NT - 2) i = GFE_NT - 2;")
        w("  const double f = x - i;")
        w("  const double* g = gfe_tab + static_cast<size_t>(s) * GFE_NT;")
        w("  return g[i] + f * (g[i + 1] - g[i]);")
        w("}")
        slot = {sp: j for j, sp in enumerate(sps)}
        for q, (R, P) in enumerate(gfe):
            dn = len(P) - len(R)
            dg = " + ".join("gfe_(%d, T)" % slot[x] for x in P) + " - (" \
                + " + ".join("gfe_(%d, T)" % slot[x] for x in R) + ")"
            w("// K_c of %s -> %s from the Gibbs free energies"
              % (" + ".join(R), " + ".join(P)))
            w("static inline double gfeKc_%d(double T) {" % q)
            w("  T = std::fmax(T, 1e-10);")
            w("  // exponent capped at 700 (see revKc_*)")
            w("  return std::exp(std::fmin(-(%s) / (constants::Rgas_kJ * T), "
              "700.0))%s;"
              % (dg, (" * std::pow(constants::boltzmann_erg * T * 1e-6, %d.0)" % -dn)
                      if dn else ""))
            w("}")
        w("")

    # -- nuclear networks ------------------------------------------------------

    def resolve_multiplicity(self):
        """rates.reactant_multiplicity (legacy -nuclearMult, opt-in): a
        reaction with a species c times among its kinetic reactants has
        flux k prod(n) / c! (identical particles counted once per set).
        Folded into the rate expression, so rhs, Jacobian, fluxes and
        heating agree; done before the reverses are built, whose rate
        k_fwd/c! / K_c then keeps detailed balance. Legacy's flag was a
        no-op (build_RHS was called again without it) and its Jacobian
        kept the curly reactants the rhs dropped."""
        if not self.reactant_multiplicity or getattr(self, "_mult_done",
                                                     False):
            return
        self._mult_done = True
        for rx in self.net.reactions:
            if rx["rate"]["kind"] != "expr":
                continue
            fac = 1
            for c in self.kin_counts(rx).values():
                fac *= math.factorial(c)
            if fac > 1:
                rx["mult"] = fac
                rx["rate"]["ast"] = {"b": "/", "l": rx["rate"]["ast"],
                                     "r": {"n": float(fac)}}

    def nucleus_Z(self, sp):
        """atomic number of a single-nucleus species (0 for e- and other
        species without charge-carrying nuclei), or None for molecules"""
        comp = self.net.species[sp]["composition"] or {}
        if not comp:
            return 0
        if sum(comp.values()) != 1:
            return None
        return ATOMIC_Z.get(next(iter(comp)))

    def resolve_screening(self):
        """rates.screening (legacy krome_stars stars_screen, never wired into
        the ODE there): the rate of every reaction whose kinetic reactants
        are two or three charged nuclei is multiplied by the plasma
        screening factor, f(Z1, Z2) or f(Z1, Z2) f(Z1 + Z2, Z3) (legacy
        3-body form). The factor depends on the composition (an extern:
        numerical dk/dn in the Jacobian)."""
        if not self.screening or getattr(self, "_screen_done", False):
            return
        self._screen_done = True
        nsc = 0
        for rx in self.net.reactions:
            if rx["rate"]["kind"] != "expr":
                continue
            kin = [r["sp"] for r in rx["reactants"] if r["kinetic"]]
            if len(kin) not in (2, 3):
                continue
            zs = [self.nucleus_Z(x) for x in kin]
            if None in zs:
                self.warnings.append("rates.screening: %s has a molecular "
                                     "reactant, not screened" % rx["eq"])
                continue
            if 0 in zs:
                continue

            def call(a, b):
                return {"f": "stars_screen", "kind": "extern",
                        "a": [{"n": float(a)}, {"n": float(b)}]}
            f = call(zs[0], zs[1])
            if len(zs) == 3:
                f = {"b": "*", "l": f, "r": call(zs[0] + zs[1], zs[2])}
            rx["rate"]["ast"] = {"b": "*", "l": rx["rate"]["ast"], "r": f}
            nsc += 1
        if nsc:
            self.net.externs.add("stars_screen")
        else:
            self.warnings.append("rates.screening: no reaction between charged "
                                 "nuclei")

    # -- reverse reactions (network-format.md 3.2) ---------------------------

    def load_thermo_nist(self):
        """data/thermoNIST.dat (NIST-JANAF Shomate): NAME -> list of ranges
        (Tmin, Tmax, [A..G]); the trailing H (= H298) is not needed, as
        A..F give the absolute enthalpy (legacy load_thermochemistry)."""
        for base in (self.data_dir.parent, self.data_dir,
                     self.data_dir.parent.parent):
            path = base / "thermoNIST.dat"
            if path.exists():
                break
        else:
            return {}
        out = {}
        for row in path.read_text().splitlines():
            srow = row.strip()
            if not srow or srow[0] == "#":
                continue
            a = srow.split(":")
            v = [float(x) for x in a[1:]]
            out.setdefault(a[0].upper(), []).append((v[0], v[1], v[2:9]))
        for k in out:
            out[k].sort()
        return out

    def thermo_all(self):
        if not hasattr(self, "_thermo_all"):
            self._thermo_all = (self.load_thermo(), self.load_thermo_nist())
        return self._thermo_all

    def emit_gRT(self, w, species):
        """G/RT = H/RT - S/R of a species at T (standard state 1 bar):
        NASA 7-coefficient polynomials (Burcat, thermo_burcat.toml) when
        present, else NIST Shomate (thermoNIST.dat). Outside the data
        the temperature is clamped to the nearest limit (legacy used
        H = S = 0 there)."""
        nasa, nist = self.thermo_all()
        key = species.upper()
        name = "gRT_%s" % self.cid(species)
        w("// G/RT of %s (%s)" % (species, "NASA, Burcat" if species in nasa
                                   else "NIST-JANAF Shomate"))
        w("static inline double %s(double T) {" % name)
        if species in nasa:
            c = nasa[species]
            lo, mid, hi = c[0], c[1], c[2]
            w("  T = std::fmin(std::fmax(T, %s), %s);"
              % (self.cnum(lo), self.cnum(hi)))
            for tag, p in (("if (T <= %s) {" % self.cnum(mid), c[10:17]),
                           ("} else {", c[3:10])):
                w("  " + tag)
                w("    const double a[7] = {%s};"
                  % ", ".join(self.cnum(x) for x in p))
                w("    const double H = a[0] + a[1]*T/2.0 + a[2]*T*T/3.0 "
                  "+ a[3]*T*T*T/4.0 + a[4]*T*T*T*T/5.0 + a[5]/T;")
                w("    const double S = a[0]*std::log(T) + a[1]*T "
                  "+ a[2]*T*T/2.0 + a[3]*T*T*T/3.0 + a[4]*T*T*T*T/4.0 "
                  "+ a[6];")
                w("    return H - S;")
            w("  }")
        else:
            rng = nist[key]
            w("  T = std::fmin(std::fmax(T, %s), %s);"
              % (self.cnum(rng[0][0]), self.cnum(rng[-1][1])))
            for j, (tmin, tmax, p) in enumerate(rng):
                cond = ("if (T <= %s) {" % self.cnum(tmax)
                        if j < len(rng) - 1 else "{")
                w("  " + cond)
                w("    const double a[7] = {%s};"
                  % ", ".join(self.cnum(x) for x in p))
                w("    const double t = T * 1e-3;")
                w("    // H [kJ/mol], S [J/mol/K] (legacy had t^5 in the "
                  "D term)")
                w("    const double H = a[0]*t + a[1]*t*t/2.0 "
                  "+ a[2]*t*t*t/3.0 + a[3]*t*t*t*t/4.0 - a[4]/t + a[5];")
                w("    const double S = a[0]*std::log(t) + a[1]*t "
                  "+ a[2]*t*t/2.0 + a[3]*t*t*t/3.0 - a[4]/(2.0*t*t) + a[6];")
                w("    return H / (constants::Rgas_kJ * T) - S / constants::Rgas_J;")
                w("  }")
        w("}")

    def resolve_reverse(self):
        """`A <-> B` (network-format.md 3.2): the reverse reaction B -> A is
        added right after each forward piece, with the forward's
        temperature range and conditions and the rate from detailed balance,
          k_rev = k_fwd / K_c,  K_c = exp(-dG/RT) (P0 / kB T)^dn,
        P0 = 1 bar, dG and dn = n_products - n_reactants of the forward
        reaction. Legacy -reverse applied the (kB T / P0)^dn factor twice
        (in the generated rate string and again in revKc)."""
        net = self.net
        if not any(rx.get("reverse") for rx in net.reactions):
            return
        nasa, nist = self.thermo_all()
        out, self._reverse = [], []
        missing = set()
        for rx in net.reactions:
            out.append(rx)
            if not rx.get("reverse"):
                continue
            sps = [r["sp"] for r in rx["reactants"]] + \
                  [p0["sp"] for p0 in rx["products"]]
            for sp in sps:
                if sp not in nasa and sp.upper() not in nist:
                    missing.add(sp)
            q = len(self._reverse)
            self._reverse.append(rx)
            rev = {"reactants": [{"sp": p0["sp"], "kinetic": True}
                                 for p0 in rx["products"]],
                   "products": [{"sp": r["sp"], "kinetic": True}
                                for r in rx["reactants"]],
                   "trange": rx["trange"], "line": rx["line"],
                   "qeff": None, "reverse_of": rx["eq"],
                   "nuclear": rx.get("nuclear", False), "lw": False}
            rev["eq"] = "%s -> %s" % (
                " + ".join(r["sp"] for r in rev["reactants"]),
                " + ".join(p0["sp"] for p0 in rev["products"]))
            if rx["rate"]["kind"] == "expr":
                rev["rate"] = {"kind": "expr", "ast": {
                    "b": "*", "l": copy.deepcopy(rx["rate"]["ast"]),
                    "r": {"f": "revKc_%d" % q, "kind": "tfunc",
                          "a": [{"s": "Tgas", "kind": "state",
                                 "unclamped": True}]}}}
            else:
                self.warnings.append("reverse of %s: the forward rate is "
                                     "stubbed, so is the reverse" % rx["eq"])
                rev["rate"] = {"kind": "stub"}
            out.append(rev)
        if missing:
            raise ValueError("<-> reactions need thermochemical data "
                             "(thermo_burcat.toml / thermoNIST.dat) for: %s"
                             % ", ".join(sorted(missing)))
        net.reactions = out
        for i, rx in enumerate(net.reactions):
            rx["id"] = i + 1
            rx.pop("group", None)
        knet.group_reactions(net)

    def emit_reverse(self, w):
        """G/RT of the species of the <-> reactions and the factor
        1/K_c(T) of each reverse rate (a function of T alone)"""
        revs = getattr(self, "_reverse", [])
        if not revs:
            return
        w("// ---- reverse reactions: 1/K_c from the Gibbs energies (spec "
          "3.2)")
        need = sorted({r["sp"] for rx in revs
                       for r in rx["reactants"] + rx["products"]})
        if any(sp in self.thermo_all()[0] for sp in need):
            self.burcat_used = True
            for ln in BURCAT_COMMENT:
                w(ln)
        for sp in need:
            self.emit_gRT(w, sp)
        for q, rx in enumerate(revs):
            P = [p0["sp"] for p0 in rx["products"]]
            R = [r["sp"] for r in rx["reactants"]]
            dn = len(P) - len(R)
            dg = " + ".join("gRT_%s(T)" % self.cid(x) for x in P) + " - (" \
                + " + ".join("gRT_%s(T)" % self.cid(x) for x in R) + ")"
            w("// 1/K_c of %s: exp(dG/RT) (kB T / 1 bar)^%d" % (rx["eq"], dn))
            w("static inline double revKc_%d(double T) {" % q)
            w("  T = std::fmax(T, 1e-10);")
            w("  // exponent capped at 700: far below the data range 1/K_c")
            w("  // would overflow, and k_fwd (often 0 there) times inf is NaN")
            if dn:
                w("  return std::exp(std::fmin(%s, 700.0)) * std::pow("
                  "constants::boltzmann_erg * T * 1e-6, %d.0);" % (dg, dn))
            else:
                w("  return std::exp(std::fmin(%s, 700.0));" % dg)
            w("}")
        w("")

    def reverse_report(self):
        """rates.check_reverse (legacy -checkReverse): for every reaction, its
        reverse in the network (exact swap), and the forward reaction
        enthalpy at 300 K [K] from the thermo data. '*' marks a missing
        reverse that is likely relevant: no electrons, >= 2 reactants and
        >= 2 products, |dH| <= 3000 K (legacy's criterion, with a proper
        electron test; legacy's used a small built-in enthalpy table and
        stopped on unknown species)."""
        nasa, nist = self.thermo_all()

        def h_RT(sp, T=300.0):
            c = nasa.get(sp)
            if c is not None:
                a = c[10:17] if T <= c[1] else c[3:10]
                return (a[0] + a[1] * T / 2 + a[2] * T * T / 3
                        + a[3] * T ** 3 / 4 + a[4] * T ** 4 / 5 + a[5] / T)
            r = nist.get(sp.upper())
            if r is not None:
                t = T * 1e-3
                a = next((q[2] for q in r if T <= q[1]), r[-1][2])
                H = (a[0] * t + a[1] * t * t / 2 + a[2] * t ** 3 / 3
                     + a[3] * t ** 4 / 4 - a[4] / t + a[5])
                return H / (constants.Rgas_kJ * T)
            return None

        rxs = self.net.reactions
        keys = {}
        for rx in rxs:
            keys.setdefault((tuple(sorted(r["sp"] for r in rx["reactants"])),
                             tuple(sorted(p0["sp"] for p0 in
                                          rx["products"]))), rx)
        lines = ["# reverse-reaction check (codegen rates.check_reverse)",
                 "# columns: flag, reaction, dH(300 K) [K], reverse in the "
                 "network"]
        for rx in rxs:
            R = tuple(sorted(r["sp"] for r in rx["reactants"]))
            P = tuple(sorted(p0["sp"] for p0 in rx["products"]))
            rev = keys.get((P, R))
            hs = [h_RT(x) for x in P + R]
            if None in hs:
                dH = None
            else:
                dH = 300.0 * (sum(hs[:len(P)]) - sum(hs[len(P):]))
            flag = " "
            if rev is None and dH is not None and "e-" not in R + P \
                    and len(R) >= 2 and len(P) >= 2 and abs(dH) <= 3000.0:
                flag = "*"
            lines.append("%s %-40s %12s  %s" % (
                flag, rx["eq"], "n/a" if dH is None else "%.1f" % dH,
                rev["eq"] if rev is not None else "NO REVERSE FOUND"))
        return "\n".join(lines) + "\n"

    # -- grain-surface chemistry (network-format.md 5.3) ---------------------

    SURF_APP2 = 9e-16     # area of a surface site, (3 A)^2 [cm^2] (Hocuk+15)
    SURF_TUNNEL_A = 1e-8  # barrier width for tunnelling [cm] (legacy ar)
    HBAR = constants.hbar_erg  # erg s
    FDIFF = {"KROME": 2.0 / 3.0, "SEMENOV": 0.77}  # E_diff / E_bind

    def surf_base(self, sp):
        """gas-phase name of a surface species (X_dust -> X, H_c_dust -> H;
        a per-bin clone maps to its original, surface.bins)"""
        if sp in getattr(self, "_bin_of", {}):
            sp = self._bin_of[sp][0]
            return self._surf_base_orig[sp]
        s = self.net.species[sp]
        if s["phase"] == "gas":
            return sp
        base = sp[:sp.lower().rfind("_" + s["phase"])]
        if s["label"] == "c" and base.endswith("_c"):
            base = base[:-2]
        return base

    def load_ebind_set(self):
        """binding energies [K] (bare, ice) of the surface.binding_energies set,
        by gas-phase name"""
        fn = {"KROME": "Ebare_ice.dat",
              "SEMENOV": "Ebare_ice_Semenov.dat"}[self.binding_energies]
        out = {}
        f = self.data_dir.parent / fn
        for line in f.read_text().splitlines():
            s = line.split("#")[0].split()
            if len(s) == 3:
                out[s[0]] = (float(s[1]), float(s[2]))
        # the data files spell some names in capitals (Semenov's HE, SI,
        # MGH2): an exact match first, then a case-insensitive one, as
        # legacy matched them. The gas names of the surface species count
        # too: an ice may be in the network without its gas phase (a
        # reduced network)
        up = {k.upper(): v for k, v in out.items()}
        names = set(self.net.species) | {self.surf_base(sp)
                                         for sp in self.net.species}
        for sp in sorted(names):
            if sp not in out and sp.upper() in up:
                out[sp] = up[sp.upper()]
        out.update(self.net.ebind)  # network ebind directives win
        return out

    def load_surface_db(self):
        """data/database/surface_2body.dat: {(reactant bases, products as
        (base, on surface)): (p, pesc_ice, pesc_bare)} (legacy
        dust_2body_rate arguments 1, 8, 9)"""
        path = self.data_dir / "surface_2body.dat"
        out = {}
        if not path.exists():
            return out
        cur = {}
        for line in path.read_text().splitlines():
            s = line.strip()
            if s.startswith("@reacts:"):
                cur = {"R": [x.strip() for x in s[8:].split(",")]}
            elif s.startswith("@prods:"):
                cur["P"] = [x.strip() for x in s[7:].split(",")]
            elif s.startswith("@rate:") and "R" in cur:
                body = s[6:].strip()
                args = body[body.index("(") + 1:body.rindex(")")]
                parts, depth, tok = [], 0, ""
                for ch in args:
                    depth += ch == "("
                    depth -= ch == ")"
                    if ch == "," and depth == 0:
                        parts.append(tok)
                        tok = ""
                    else:
                        tok += ch
                parts.append(tok)
                num = [eval(x.lower().replace("d", "e"), {"__builtins__": {}})
                       if "(" not in x else None for x in parts]

                def base(x):
                    return x[:-5] if x.lower().endswith("_dust") else x
                key = (tuple(sorted(base(x) for x in cur["R"])),
                       tuple(sorted((base(x), x.lower().endswith("_dust"))
                                    for x in cur["P"])))
                out[key] = (num[0], num[7], num[8])
        return out

    def bin_name(self, sp, j):
        """per-bin clone of a surface species: X_dust -> X_b<j>_dust"""
        ph = self.net.species[sp]["phase"]
        i = sp.lower().rfind("_" + ph)
        return "%s_b%d%s" % (sp[:i], j, sp[i:])

    def expand_surface_bins(self):
        """surface.bins = NB (legacy per-bin surface species): every
        surface species becomes NB clones, one per dust bin of the runtime
        DustBins population (which must have NB bins: flag -107), and every
        reaction involving surface species one reaction per bin. Rate
        expressions that read a surface density see the sum over bins."""
        net = self.net
        self._bin_of, self._surf_base_orig, self._orig_label = {}, {}, {}
        nb = self.surface_bins
        if not nb:
            return
        if self.gas_total:
            raise ValueError("surface.bins is not supported with "
                             "dust.gas_total")
        if net.surface_grains is not None:
            raise ValueError("surface.bins needs the binned population; "
                             "the network declares surface_grains single("
                             "...) (line %d)" % net.surface_grains["line"])
        surf_sp = [s for s in sorted(net.species)
                   if net.species[s]["phase"] != "gas"]
        if not surf_sp:
            self.warnings.append("surface.bins: no surface species")
            return
        for s in surf_sp:
            self._surf_base_orig[s] = self.surf_base(s)
            self._orig_label[s] = net.species[s]["label"]
        clones = {s: [self.bin_name(s, j) for j in range(nb)]
                  for s in surf_sp}
        for s in surf_sp:
            rec = net.species.pop(s)
            for j, nm in enumerate(clones[s]):
                r2 = dict(rec)
                r2["name"] = nm
                net.species[nm] = r2
                self._bin_of[nm] = (s, j)

        def sum_bins(ast):
            for nd in knet.walk(ast):
                if nd.get("kind") == "density" and \
                        nd["a"][0].get("sp") in clones:
                    terms = [{"f": "n", "kind": "density",
                              "a": [{"sp": c}]} for c in clones[
                                  nd["a"][0]["sp"]]]
                    tot = terms[0]
                    for t in terms[1:]:
                        tot = {"b": "+", "l": tot, "r": t}
                    nd.clear()
                    nd.update(tot)
        for v in net.vars.values():
            sum_bins(v)
        out = []
        for rx in net.reactions:
            sps = [q["sp"] for q in rx["reactants"] + rx["products"]]
            if not any(x in clones for x in sps):
                if rx["rate"].get("ast") is not None:
                    sum_bins(rx["rate"]["ast"])
                out.append(rx)
                continue
            for j in range(nb):
                r2 = copy.deepcopy(rx)
                for q in r2["reactants"] + r2["products"]:
                    if q["sp"] in clones:
                        q["sp"] = clones[q["sp"]][j]
                r2["bin"] = j
                r2["eq"] = "%s -> %s" % (
                    " + ".join(q["sp"] + ("" if q["kinetic"] else "{}")
                               for q in r2["reactants"]),
                    " + ".join(q["sp"] for q in r2["products"]))
                if r2["rate"].get("ast") is not None:
                    sum_bins(r2["rate"]["ast"])
                out.append(r2)
        net.reactions = out
        for fr in net.freeze:
            fr["species"] = [c for x in fr["species"]
                             for c in clones.get(x, [x])]
        for i, rx in enumerate(net.reactions):
            rx["id"] = i + 1
            rx.pop("group", None)
        knet.group_reactions(net)

    def resolve_surface(self):
        """grain-surface chemistry on the dust population (network-format.md
        5.3): every surface(...) reaction becomes a generated rate function.
          gas X -> X_dust             adsorption (sticking, HM79)
          X_dust -> X                 thermal (+ cosmic-ray) desorption
          A_dust + B_dust -> products Langmuir-Hinshelwood, with reaction-
                                      diffusion competition
          H_dust <-> H_c_dust, H_dust + H_c_dust, H_c_dust + H_c_dust
                                      chemisorption (legacy tables)
        On the effective MRN population, or area-weighted over the bins of
        a DustBins population at runtime."""
        net = self.net
        self._surf = []
        todo = [rx for rx in net.reactions if rx["rate"]["kind"] == "surface"]
        if not todo:
            return
        if self.gas_total and any(
                sum(1 for q in rx["reactants"]
                    if net.species[q["sp"]]["phase"] != "gas") == 2
                for rx in todo):
            raise ValueError("two-body surface chemistry is not supported "
                             "with dust.gas_total")
        eb = self.load_ebind_set()
        db = self.load_surface_db()
        need_eb = set()
        extra = []   # (after rx, new rx) reactive-desorption branches

        def call(name, args=()):
            self._gen_ext.add(name)
            net.externs.add(name)
            return {"f": name, "kind": "extern", "a": list(args)}

        def surf(sp):
            return net.species[sp]["phase"] != "gas"

        def chem(sp):
            orig = self._bin_of[sp][0] if sp in self._bin_of else sp
            return surf(sp) and self._orig_label.get(orig, net.species[
                sp]["label"]) == "c"

        for rx in todo:
            R = [q["sp"] for q in rx["reactants"]]
            P = [q["sp"] for q in rx["products"]]
            r = rx["rate"]
            kind = None
            bn = {"n": float(rx.get("bin", -1))}  # -1: the whole population
            if len(R) == 1 and len(P) == 1:
                a, b = R[0], P[0]
                if not surf(a) and surf(b) and not chem(b) \
                        and self.surf_base(b) == a:
                    kind, spec = "ads", {"sp": a}
                elif surf(a) and not chem(a) and not surf(b) and (
                        self.surf_base(a) == b
                        or net.species[a]["composition"]
                        == net.species[b]["composition"]):
                    # into its gas counterpart, or an isomer of it (Semenov:
                    # HC2NC_DUST -> HC3N): the adsorbed species' binding
                    # energy and mass
                    kind, spec = "des", {"sp": a}
                    need_eb.add(self.surf_base(a))
                elif self.surf_base(a) == "H" and self.surf_base(b) == "H" \
                        and surf(a) and surf(b) and chem(a) != chem(b):
                    kind = "chem"
                    spec = {"mode": "CP" if chem(a) else "PC"}
            elif len(R) == 2 and all(surf(x) for x in R):
                if any(chem(x) for x in R):
                    if [self.surf_base(x) for x in R] == ["H", "H"] and \
                            [self.surf_base(x) for x in P] == ["H2"]:
                        kind = "chem2"
                        spec = {"mode": "CC" if all(chem(x) for x in R)
                                else "PCCP"}
                else:
                    kind = "2body"
                    need_eb.update(self.surf_base(x) for x in R)
                    key = (tuple(sorted(self.surf_base(x) for x in R)),
                           tuple(sorted((self.surf_base(x), surf(x))
                                        for x in P)))
                    ent = db.get(key)
                    spec = {"A": self.surf_base(R[0]),
                            "B": self.surf_base(R[1]), "branch":
                            r.get("branch", 1.0), "rx": rx, "key": key}
                    if "Ea" in r:
                        spec["Ea"] = r["Ea"]
                    elif "prob" in r:
                        spec["prob"] = r["prob"]
                    elif ent is not None:
                        spec["prob"] = ent[0]
                    else:
                        raise ValueError(
                            "surface(auto) %s (line %s): not in "
                            "data/database/surface_2body.dat; give "
                            "surface(Ea=...) or surface(prob=...)"
                            % (rx["eq"], rx["line"]))
                    spec["gas_prod"] = not any(surf(x) for x in P)
            if kind is None:
                raise ValueError(
                    "surface(...) %s (line %s): not a recognized surface "
                    "process (network-format.md 5.3)" % (rx["eq"], rx["line"]))
            q = len(self._surf)
            spec["kind"] = kind
            self._surf.append(spec)
            if kind != "ads" and r.get("stick") is not None:
                raise ValueError("%s: surface(stick=...) is only for "
                                 "adsorption (gas X -> X_dust)" % rx["eq"])
            if kind == "ads" and r.get("stick") is not None:
                spec["kind"] = "ads_s"
                r2 = {"b": "*", "l": call("surf_adsorption_ns",
                                          [{"sp": spec["sp"]}, bn]),
                      "r": r["stick"]}
            elif kind == "ads":
                r2 = call("surf_adsorption", [{"sp": spec["sp"]}, bn])
            elif kind == "des":
                r2 = call("surf_desorption", [{"sp": spec["sp"]}, bn])
            elif kind in ("chem", "chem2"):
                r2 = call("surf_chem_%s" % spec["mode"], [bn])
            else:
                r2 = call("surf_2body_%d" % q, [{"n": 0.0}, bn])
            rx["rate"] = {"kind": "expr", "ast": r2}
        # reactive desorption: fractions per two-body reaction
        two = [s for s in self._surf if s["kind"] == "2body"]
        keys = {s["key"] for s in two}
        for s in two:
            s["cd"] = (0.0, 0.0)  # desorbing fraction (bare, ice)
            rx = s["rx"]
            if s["gas_prod"] and self.reactive_desorption != "NONE":
                raise ValueError("%s: with surface.reactive_desorption the gas "
                                 "branch is generated; write only the "
                                 "surface-product reaction" % rx["eq"])
            ent = db.get(s["key"])
            if self.reactive_desorption == "NONE":
                # legacy database split, only when both branches are in the
                # network (legacy applied 1 - pesc to the surface branch
                # even without the gas branch, losing that flux)
                if s["gas_prod"]:
                    if ent is None:
                        self.warnings.append("%s: gas-product branch without "
                                             "a database fraction: 1"
                                             % rx["eq"])
                        s["mode"] = 0
                    else:
                        s["cd"] = (ent[2], ent[1])
                        s["mode"] = 2
                else:
                    gkey = (s["key"][0], tuple(sorted(
                        (b, False) for b, _ in s["key"][1])))
                    gent = db.get(gkey)
                    if gkey in keys and gent is not None:
                        s["cd"] = (gent[2], gent[1])
                        s["mode"] = 1
                    else:
                        s["mode"] = 0
                continue
            # automatic gas branch
            if self.reactive_desorption == "FIXED":
                gkey = (s["key"][0], tuple(sorted(
                    (b, False) for b, _ in s["key"][1])))
                gent = db.get(gkey)
                if gent is None:
                    s["mode"] = 0
                    continue
                s["cd"] = (gent[2], gent[1])
            else:  # MINISSALE
                cd = self.minissale(s)
                if cd is None:
                    s["mode"] = 0
                    continue
                s["cd"] = cd
            s["mode"] = 1
            prods = [{"sp": self.surf_base(q["sp"]), "kinetic": True}
                     for q in rx["products"]]
            missing = [p0["sp"] for p0 in prods if p0["sp"] not in net.species]
            if missing:
                self.warnings.append("%s: no gas species %s for the "
                                     "desorbing branch" % (rx["eq"],
                                                           ", ".join(missing)))
                s["mode"] = 0
                continue
            q = self._surf.index(s)
            g = {"reactants": [dict(x) for x in rx["reactants"]],
                 "products": prods, "trange": rx["trange"],
                 "line": rx["line"], "qeff": None, "lw": False,
                 "rate": {"kind": "expr",
                          "ast": call("surf_2body_%d" % q, [
                              {"n": 2.0}, {"n": float(rx.get("bin", -1))}])}}
            g["eq"] = "%s -> %s" % (" + ".join(x["sp"] for x in
                                              g["reactants"]),
                                    " + ".join(x["sp"] for x in prods))
            if "bin" in rx:
                g["bin"] = rx["bin"]
            extra.append((rx, g))
        for s in two:
            if s.get("mode"):
                s["rx"]["rate"]["ast"]["a"][0] = {"n": float(s["mode"])}
        if extra:
            out = []
            after = {id(a): [] for a, _ in extra}
            for a, g in extra:
                after[id(a)].append(g)
            for rx in net.reactions:
                out.append(rx)
                out.extend(after.get(id(rx), []))
            net.reactions = out
            for i, rx in enumerate(net.reactions):
                rx["id"] = i + 1
                rx.pop("group", None)
            knet.group_reactions(net)
        miss = sorted(x for x in need_eb if x not in eb)
        if miss:
            other = ("SEMENOV" if self.binding_energies == "KROME"
                     else "KROME")
            raise ValueError("surface chemistry: no binding energy (%s set) "
                             "for %s: the network may need the other set "
                             "(surface.binding_energies = \"%s\"; its "
                             "header says), or `ebind X = E_bare, E_ice` "
                             "lines" % (self.binding_energies,
                                        ", ".join(miss), other))
        self._surf_eb = eb
        if any(s["kind"] == "des" for s in self._surf) and self.cr_desorption \
                and self.cr_desorption not in net.params and \
                self.cr_desorption not in knet.IMPLICIT_PARAMS:
            raise ValueError("surface.cr_desorption=%s: the network declares no "
                             "such parameter (param %s = ...)"
                             % (self.cr_desorption, self.cr_desorption))

    def minissale(self, s):
        """chemical (reactive) desorption fraction of a one-product surface
        reaction, Minissale et al. 2016 (A&A 585, A24): exp(-E_b N /
        (eps dH_R)), eps = ((M - m)/(M + m))^2, m = 120 amu, N = 3 atoms;
        on water ice their suggested reduction, /10. dH_R [K]: the
        exothermicity from Burcat's enthalpies at 300 K."""
        rx = s["rx"]
        prods = [self.surf_base(q["sp"]) for q in rx["products"]]
        if len(prods) != 1:
            self.warnings.append("%s: Minissale reactive desorption is for "
                                 "one product, none" % rx["eq"])
            return None
        nasa, _ = self.thermo_all()

        def h(sp):
            c = nasa.get(sp)
            if c is None:
                return None
            # c: [Tlow, Tcommon, Thigh, 7 upper-range, 7 lower-range
            # coefficients] (load_thermo); the range that holds 300 K
            a = c[10:17] if 300.0 <= c[1] else c[3:10]
            T = 300.0
            return T * (a[0] + a[1] * T / 2 + a[2] * T * T / 3
                        + a[3] * T ** 3 / 4 + a[4] * T ** 4 / 5 + a[5] / T)
        hs = [h(x) for x in [self.surf_base(q["sp"]) for q in
                             rx["reactants"]] + prods]
        if None in hs:
            self.warnings.append("%s: no thermochemical data, no reactive "
                                 "desorption" % rx["eq"])
            return None
        exo = sum(hs[:-1]) - hs[-1]
        if exo <= 0.0:
            return (0.0, 0.0)
        P = prods[0]
        sp = self.net.species[P]
        natoms = sum((sp["composition"] or {}).values())
        M, m = sp["mass_amu"], 120.0
        eps = ((M - m) / (M + m)) ** 2
        Eb = self._eb_of(P)[0]
        cd = math.exp(-Eb * 3 * natoms / (eps * exo))
        return (cd, cd / 10.0)

    def _eb_of(self, base):
        return self.load_ebind_set().get(base, (0.0, 0.0))

    def emit_surface(self, w):
        surfs = getattr(self, "_surf", [])
        nb = self.surface_bins if getattr(self, "_bin_of", None) else 0
        w("// surface.bins: per-bin surface species (0: one population);")
        w("// the runtime DustBins must then have SURF_NB bins (flag -107)")
        w("constexpr int SURF_NB = %d;" % nb)
        if not surfs:
            return
        net = self.net
        sp_sorted = sorted(net.species)
        eb = self._surf_eb
        base = {s: self.surf_base(s) for s in sp_sorted}
        w("// ---- grain-surface chemistry (network-format.md 5.3; binding "
          "energies %s," % self.binding_energies)
        w("// nu %s, E_diff = %.4g E_bind, reaction-diffusion competition %s)"
          % (self.surface_nu, self.FDIFF[self.binding_energies],
             "on" if self.surface_competition else "off"))
        w("constexpr double SURF_APP2 = %s;  // site area [cm^2]"
          % self.cnum(self.SURF_APP2))
        w("constexpr double SURF_FDIFF = %s;"
          % self.cnum(self.FDIFF[self.binding_energies]))
        w("inline const double surf_Eb_bare[NSPEC] = {%s};"
          % ", ".join(self.cnum(eb.get(base[s], (0.0, 0.0))[0])
                      for s in sp_sorted))
        w("inline const double surf_Eb_ice[NSPEC] = {%s};"
          % ", ".join(self.cnum(eb.get(base[s], (0.0, 0.0))[1])
                      for s in sp_sorted))
        ice = [s for s in sp_sorted if net.species[s]["phase"] != "gas"
               and base[s] == "H2O" and net.species[s]["label"] != "c"]
        self._surf_ice = ice  # the direct density dependence (surf_fice_)
        w("// characteristic vibration frequency [1/s]: %s"
          % ("legacy 1e12" if self.surface_nu == "LEGACY" else
             "Hasegawa et al. 1992, sqrt(2 n_s E_b / (pi^2 m))"))
        w("static inline double surf_nu_(int sp, double Eb) {")
        if self.surface_nu == "LEGACY":
            w("  (void)sp; (void)Eb;")
            w("  return 1e12;")
        else:
            w("  return std::sqrt(2.0 / SURF_APP2 * Eb * constants::boltzmann_erg / "
              "(constants::pi * constants::pi * species_mass_g[sp]));")
        w("}")
        w("// surface sites per volume [cm^-3]: 4 pi sum n_d a^2 / app2, the")
        w("// bins' area moment or the MRN population holding rho d2g")
        w("// (bin >= 0: that bin alone, surface.bins)")
        w("static inline double surf_sites_(const double* n, const Params& p,"
          " int bin) {")
        sg = self.net.surface_grains
        if sg is not None:
            # one grain size (network surface_grains single(...)), legacy
            # Semenov framework: S = 4 pi n_d a^2 / app2
            w("  // surface_grains single(...): n_d grains cm^-3 of radius^2 "
              "a2")
            w("  (void)n; (void)bin;")
            w("  return 4.0 * constants::pi * (%s) * (%s) / SURF_APP2;"
              % (self.cexpr(sg["n"]), self.cexpr(sg["a2"])))
            w("}")
        if sg is None:
            w("  if (bin >= 0)")
            w("    return 4.0 * constants::pi * p.hk.dust.aw[bin] * rho_(n) "
              "* p.d2g / SURF_APP2;")
            w("  if (p.hk.dust.nb > 0) {")
            w("    double s = 0.0;")
            w("    for (int b = 0; b < p.hk.dust.nb; ++b) "
              "s += p.hk.dust.aw[b];")
            w("    return 4.0 * constants::pi * s * rho_(n) * p.d2g / "
              "SURF_APP2;")
            w("  }")
            w("  const double amin = 5e-7, amax = 2.5e-5, pexp = -3.5, "
              "rho0 = 3.0;")
            w("  const double p3 = pexp + 3.0, p4 = pexp + 4.0;")
            w("  const double na2 = rho_(n) * p.d2g * 3.0 / (4.0 * "
              "constants::pi * rho0) * (p4 / p3) * (std::pow(amax, p3) - "
              "std::pow(amin, p3)) / (std::pow(amax, p4) - std::pow(amin, "
              "p4));")
            w("  return 4.0 * constants::pi * na2 / SURF_APP2;")
            w("}")
        w("// water-ice coverage (monolayers, capped at 1; legacy did not")
        w("// clip)")
        w("static inline double surf_fice_(const double* n, double S, "
          "int bin) {")
        if ice and nb:
            w("  constexpr int iw[%d] = {%s};" % (nb, ", ".join(
                "IDX_%s" % self.cid(x) for x in sorted(
                    ice, key=lambda x: self._bin_of[x][1]))))
            w("  if (!(S > 0.0) || bin < 0) return 0.0;")
            w("  return std::fmin(std::fmax(n[iw[bin]] / S, 0.0), 1.0);")
        elif ice:
            # the ice density in either representation (--gas-total: the
            # slot holds gas + ice)
            w("  (void)bin;")
            w("  if (!(S > 0.0)) return 0.0;")
            w("  return std::fmin(std::fmax((%s) / S, 0.0), 1.0);"
              % self.nden(ice[0]))
        else:
            w("  (void)n; (void)S; (void)bin;")
            w("  return 0.0;  // no water ice species")
        w("}")
        w("// area-weighted average of f(Tdust) over the bins, or f at the")
        w("// live Tdust of the effective population")
        w("template <class F>")
        w("static inline double surf_avg_(double Tgas, const double* n, "
          "const Params& p, int bin, F f) {")
        if sg is not None:  # one grain size: the cell's Tdust
            w("  (void)bin;")
        else:
            w("  if (bin >= 0) return f(p.hk.dust.Td[bin]);")
            w("  if (p.hk.dust.nb > 0) {")
            w("    double s = 0.0, w = 0.0;")
            w("    for (int b = 0; b < p.hk.dust.nb; ++b) {")
            w("      s += p.hk.dust.aw[b] * f(p.hk.dust.Td[b]);")
            w("      w += p.hk.dust.aw[b];")
            w("    }")
            w("    return s / w;")
            w("  }")
        w("  return f(table_Tdust(Tgas, n, p));")
        w("}")
        kinds = {s["kind"] for s in surfs}
        if "ads" in kinds:
            w("// adsorption [1/s]: sum_b n_d,b pi a_b^2 v_th S(Tgas, Td,b)")
            w("// (HM79 sticking; the real d2g in both populations)")
            w("inline double surf_adsorption(int sp, double binf, "
              "double Tgas, const double* n, const Params& p) {")
            w("  const int bin = static_cast<int>(binf);")
            w("  const double Tg100 = Tgas * 1e-2;")
            w("  auto stick = [&](double Td) {")
            w("    return 1.0 / (1.0 + 0.4 * std::sqrt(Tg100 + Td * 1e-2) + "
              "0.2 * Tg100 + 0.08 * Tg100 * Tg100);")
            w("  };")
            w("  const double v = %s * std::sqrt(Tgas / "
              "species_mass_g[sp]);  // sqrt(8 kB / pi)" % self.cnum(
                  math.sqrt(8.0 * constants.boltzmann_erg / constants.pi)))
            w("  // pi sum n_d a^2 = S app2 / 4")
            w("  return v * surf_sites_(n, p, bin) * SURF_APP2 / 4.0 * "
              "surf_avg_(Tgas, n, p, bin, stick);")
            w("}")
        if "ads_s" in kinds:
            w("// adsorption without the sticking factor, v_th sum n_d pi "
              "a^2 [1/s]:")
            w("// times the reaction's surface(stick=...)")
            w("inline double surf_adsorption_ns(int sp, double binf, "
              "double Tgas, const double* n, const Params& p) {")
            w("  const int bin = static_cast<int>(binf);")
            w("  const double v = %s * std::sqrt(Tgas / "
              "species_mass_g[sp]);  // sqrt(8 kB / pi)" % self.cnum(
                  math.sqrt(8.0 * constants.boltzmann_erg / constants.pi)))
            w("  return v * surf_sites_(n, p, bin) * SURF_APP2 / 4.0;")
            w("}")
        if "des" in kinds:
            cr = self.cr_desorption
            w("// thermal desorption nu exp(-E/Td), bare and ice surfaces by")
            w("// the ice coverage%s" % (
                ", + cosmic-ray spot heating (Hasegawa & Herbst 1993):"
                if cr else ""))
            if cr:
                w("// f70 nu exp(-E / 70 K), f70 = 3.16e-19 %s / 1.3e-17"
                  % cr)
            w("inline double surf_desorption(int sp, double binf, "
              "double Tgas, const double* n, const Params& p) {")
            w("  const int bin = static_cast<int>(binf);")
            w("  const double fi = surf_fice_(n, surf_sites_(n, p, bin), "
              "bin);")
            w("  const double Eb = surf_Eb_bare[sp], Ei = surf_Eb_ice[sp];")
            w("  const double nb = surf_nu_(sp, Eb), ni = surf_nu_(sp, Ei);")
            w("  auto kT = [&](double Td) {")
            w("    return (1.0 - fi) * nb * std::exp(-Eb / Td) + fi * ni * "
              "std::exp(-Ei / Td);")
            w("  };")
            if cr:
                w("  return surf_avg_(Tgas, n, p, bin, kT) + 3.16e-19 * p.%s"
                  " / 1.3e-17 * kT(70.0);" % cr)
            else:
                w("  return surf_avg_(Tgas, n, p, bin, kT);")
            w("}")
        for q, s in enumerate(surfs):
            if s["kind"] != "2body":
                continue
            A, B = s["A"], s["B"]
            iA = "IDX_%s" % self.cid(A if A in net.species else next(
                x for x in sp_sorted if base[x] == A))
            iB = "IDX_%s" % self.cid(B if B in net.species else next(
                x for x in sp_sorted if base[x] == B))
            rx = s["rx"]
            ident = [x["sp"] for x in rx["reactants"]]
            same = ident[0] == ident[1]
            if "Ea" in s:
                def mass(b):
                    sp = b if b in net.species else next(
                        x for x in sp_sorted if base[x] == b)
                    return net.species[sp]["mass_amu"] * constants.amu
                mu = mass(A) * mass(B) / (mass(A) + mass(B))
                tun = math.exp(-2.0 * self.SURF_TUNNEL_A / self.HBAR
                               * math.sqrt(2.0 * mu * s["Ea"]
                                         * constants.boltzmann_erg))
                if s["Ea"] <= 0:
                    prob = "1.0"
                elif self.surface_tunnelling:
                    prob = ("std::fmax(std::exp(-%s / Td), %s)"
                            % (self.cnum(s["Ea"]), self.cnum(tun)))
                else:  # --no-surface-tunnelling: thermal only
                    prob = "std::exp(-%s / Td)" % self.cnum(s["Ea"])
            else:
                prob = self.cnum(s["prob"])
            cdb, cdi = s["cd"]
            w("// %s: Langmuir-Hinshelwood%s%s" % (
                rx["eq"], ", identical reactants (1/2)" if same else "",
                "" if not s.get("mode") else
                ", reactive desorption (bare %.3g, ice %.3g)" % (cdb, cdi)))
            w("// mode 0: whole rate; 1: the surface branch (1 - f); 2: the")
            w("// desorbing branch (f), f = f_bare CD_bare + f_ice CD_ice")
            w("inline double surf_2body_%d(double mode, double binf, "
              "double Tgas, const double* n, const Params& p) {" % q)
            w("  const int bin = static_cast<int>(binf);")
            w("  const double S = surf_sites_(n, p, bin);")
            w("  if (!(S > 0.0)) return 0.0;")
            w("  const double fi = surf_fice_(n, S, bin), fb = 1.0 - fi;")
            w("  auto k1 = [&](double Td, double EA, double EB) {")
            w("    const double nA = surf_nu_(%s, EA), nB = surf_nu_(%s, EB);"
              % (iA, iB))
            w("    const double hA = nA * std::exp(-SURF_FDIFF * EA / Td);")
            w("    const double hB = nB * std::exp(-SURF_FDIFF * EB / Td);")
            w("    const double prob = %s;" % prob)
            if self.surface_competition:
                w("    // reaction-diffusion competition (Chang+07, Garrod &")
                w("    // Pauly 11)")
                w("    const double Pr = std::fmax(nA, nB) * prob;")
                w("    const double dA = nA * std::exp(-EA / Td), dB = nB * "
                  "std::exp(-EB / Td);")
                w("    const double den = Pr + hA + hB + dA + dB;")
                w("    const double kap = den > 0.0 ? Pr / den : 0.0;")
            else:
                w("    const double kap = prob;  // legacy: no competition")
            w("    return kap * (hA + hB);")
            w("  };")
            w("  auto kT = [&](double Td) {")
            w("    return fb * k1(Td, surf_Eb_bare[%s], surf_Eb_bare[%s])"
              % (iA, iB))
            w("        + fi * k1(Td, surf_Eb_ice[%s], surf_Eb_ice[%s]);"
              % (iA, iB))
            w("  };")
            w("  double k = %s * surf_avg_(Tgas, n, p, bin, kT) / S;"
              % self.cnum(s["branch"] * (0.5 if same else 1.0)))
            w("  const double f = fb * %s + fi * %s;"
              % (self.cnum(cdb), self.cnum(cdi)))
            w("  if (mode == 1.0) k *= 1.0 - f;")
            w("  else if (mode == 2.0) k *= f;")
            w("  return k;")
            w("}")
        chem = {s["mode"] for s in surfs if s["kind"] in ("chem", "chem2")}
        if chem:
            tabs = self.load_chemisorption()
            T0, dT, _ = tabs["meta"]
            w("// H chemisorption rates [1/s] vs Tdust (legacy")
            w("// surface_chemisorption_rates.dat; Cazaux & Tielens 2004,")
            w("// Iqbal+ 2010, Hocuk+ 2014): linear on a uniform grid, clamped")
            w("// (legacy indexed with (N-1)/(N dx) and had no bounds check)")
            w("constexpr double SURF_CHEM_T0 = %s, SURF_CHEM_DT = %s;"
              % (self.cnum(T0), self.cnum(dT)))
            for m in ("PC", "CP", "CC"):
                self.emit_data_array(w, "surf_chem_tab_%s" % m, tabs[m])
            N = len(tabs["PC"])
            w("static inline double surf_chem_rate_(const double* t, double "
              "Td) {")
            w("  double x = (Td - SURF_CHEM_T0) / SURF_CHEM_DT;")
            w("  x = std::fmin(std::fmax(x, 0.0), %d.0);" % (N - 1))
            w("  int i = static_cast<int>(x);")
            w("  if (i > %d) i = %d;" % (N - 2, N - 2))
            w("  return t[i] + (x - i) * (t[i + 1] - t[i]);")
            w("}")
            for m in sorted(chem):
                if m in ("PC", "CP"):
                    w("// H %s: %s [1/s]" % (
                        "physisorbed -> chemisorbed" if m == "PC" else
                        "chemisorbed -> physisorbed", m))
                    w("inline double surf_chem_%s(double binf, double Tgas, "
                      "const double* n, const Params& p) {" % m)
                    w("  return surf_avg_(Tgas, n, p, static_cast<int>(binf),"
                      " [](double Td) { return surf_chem_rate_("
                      "surf_chem_tab_%s, Td); });" % m)
                    w("}")
                elif m == "PCCP":
                    w("// H + H_c -> H2: (PC + CP) / S (legacy listed the two")
                    w("// channels and used only the first, CP)")
                    w("inline double surf_chem_PCCP(double binf, double Tgas, "
                      "const double* n, const Params& p) {")
                    w("  const int bin = static_cast<int>(binf);")
                    w("  const double S = surf_sites_(n, p, bin);")
                    w("  if (!(S > 0.0)) return 0.0;")
                    w("  return surf_avg_(Tgas, n, p, bin, [](double Td) { "
                      "return "
                      "surf_chem_rate_(surf_chem_tab_PC, Td) + "
                      "surf_chem_rate_(surf_chem_tab_CP, Td); }) / S;")
                    w("}")
                else:
                    w("// H_c + H_c -> H2: CC / S")
                    w("inline double surf_chem_CC(double binf, double Tgas, "
                      "const double* n, const Params& p) {")
                    w("  const int bin = static_cast<int>(binf);")
                    w("  const double S = surf_sites_(n, p, bin);")
                    w("  if (!(S > 0.0)) return 0.0;")
                    w("  return surf_avg_(Tgas, n, p, bin, [](double Td) { "
                      "return "
                      "surf_chem_rate_(surf_chem_tab_CC, Td); }) / S;")
                    w("}")
        w("")

    def load_chemisorption(self):
        """data/surface_chemisorption_rates.dat: blocks MODE, N, Tmin, dT,
        then N rates [1/s]"""
        path = self.data_dir.parent / "surface_chemisorption_rates.dat"
        rows = [ln.strip() for ln in path.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
        out, i = {}, 0
        while i < len(rows):
            mode = rows[i]
            N = int(float(rows[i + 1].replace("d", "e")))
            T0 = float(rows[i + 2].replace("d", "e"))
            dT = float(rows[i + 3].replace("d", "e"))
            out[mode] = [float(x.replace("d", "e"))
                         for x in rows[i + 4:i + 4 + N]]
            out["meta"] = (T0, dT, N)
            i += 4 + N
        return out

    PARAMS_SLOT = "//@@PARAMS@@"
    H2OP_SLOT = "//@@H2OP@@"

    def h2op_lines(self):
        """nH2o_ / nH2p_ / h2_ofrac_: ortho and para H2 for the processes
        that need them (the H2or/H2pa colliders, the adiabatic index): the
        network's own H2_ortho / H2_para species, else H2 split by the
        parameter opratio (both when the network has both). Emitted only
        when used (self._need_h2op), so opratio exists only then."""
        h2 = self.h2_forms()
        out = ["// ortho / para H2: the network's own species, else H2 split "
               "by opratio"]
        for fn, spin, frac in (("nH2o_", "H2_ortho",
                                "p.opratio / (p.opratio + 1.0)"),
                               ("nH2p_", "H2_para", "1.0 / (p.opratio + 1.0)")):
            parts = []
            if spin in h2:
                parts.append("n[IDX_%s]" % self.cid(spin))
            if "H2" in h2:
                parts.append("n[IDX_H2] * p.opratio / (p.opratio + 1.0)"
                             if spin == "H2_ortho" else
                             "n[IDX_H2] / (p.opratio + 1.0)")
            out.append("inline double %s(const double* n, const Params& p) {"
                       % fn)
            out.append("  (void)n; (void)p;")
            out.append("  return %s;" % (" + ".join(parts) or "0.0"))
            out.append("}")
        out.append("// the ortho fraction of H2 (a = opr / (opr + 1))")
        out.append("inline double h2_ofrac_(const double* n, const Params& p) "
                   "{")
        out.append("  (void)n; (void)p;")
        if h2 == ["H2"]:
            out.append("  return p.opratio / (p.opratio + 1.0);")
        else:
            out.append("  return nH2o_(n, p) / std::fmax(nH2o_(n, p) "
                       "+ nH2p_(n, p), 1e-300);")
        out.append("}")
        return out

    def active_params(self, sec):
        """the runtime parameters this network has: ALWAYS_PARAMS, the
        network's own params, every p.<name> the generated code reads, and
        those the runtime modules read for the configured processes and
        the externs used; Tdust only if something consumes table_Tdust
        and no dust table supplies it."""
        net = self.net
        known = {nm for nm, _, _ in UNIVERSAL_PARAMS} | set(net.params)
        used = set(ALWAYS_PARAMS) | set(net.params)
        pat = re.compile(r"\bp\.([A-Za-z_]\w*)\b")
        for lines in sec.values():
            for ln in lines:
                used.update(m for m in pat.findall(ln) if m in known)
        used.discard("Tdust")  # decided below (the stub reads it lazily)
        for f, names in PROCESS_PARAMS.items():
            if self.thermal.get(f):
                used.update(names)
        for f in net.externs:
            used.update(EXTERN_PARAMS.get(f, ()))
        wants_td = (any(f in net.externs for f in TDUST_EXTERNS)
                    or bool(getattr(self, "_surf", None))
                    or any(self.thermal.get(f) for f in TDUST_PROCESSES)
                    or any(nd.get("kind") == "state" and nd.get("s")
                           == "Tdust" for ast in self._all_asts()
                           for nd in self.deep_walk(ast)))
        if "Tdust" in net.params or (wants_td and not self.dust_tables):
            used.add("Tdust")
        return used

    def _all_asts(self):
        net = self.net
        out = [rx["rate"]["ast"] for rx in net.reactions
               if rx["rate"].get("ast") is not None]
        out += [rx["qeff"] for rx in net.reactions if rx.get("qeff")]
        out += [c["expr"] for c in net.cooling + net.heating]
        out += list(net.vars.values())
        return out

    def emit_params(self, sec):
        """the Params struct and its name tables, in place of PARAMS_SLOT"""
        net = self.net
        used = self.active_params(sec)
        univ = [(nm, dv, cm) for nm, dv, cm in UNIVERSAL_PARAMS
                if nm not in net.params and nm in used]
        self.pnames = sorted(net.params) + [nm for nm, _, _ in univ]
        self.param_required = {nm for nm, dv, _ in univ if dv is None}
        self.param_required_by = {}
        req_by = {}
        for f, names in PROCESS_PARAMS.items():
            if self.thermal.get(f):
                for nm in names:
                    req_by.setdefault(nm, []).append(
                        self.thermal_option_name(f))
        for f in sorted(net.externs):
            for nm in EXTERN_PARAMS.get(f, ()):
                req_by.setdefault(nm, []).append(f + "()")
        if self.num2col == "length":
            req_by.setdefault("Lcol", []).append("radiation.num2col = 'length'")
        for nm in self.param_required:
            self.param_required_by[nm] = ", ".join(req_by.get(nm, ["?"]))
        out = ["// binned-dust hook (dust-chemistry coupling): when a DustBins",
               "// population drives a run, the solver installs its per-bin",
               "// temperatures and area weights (xdust_w * a^2) here, and",
               "// every dust-temperature consumer (table_Tdust, sticking,",
               "// evaporation) evaluates area-weighted bin averages instead",
               "// of the table/parameter Tdust. nb = 0: no bins installed.",
               "struct DustBinsHook {",
               "  int nb = 0;",
               "  const double* Td = nullptr;  // per-bin dust temperature",
               "  const double* aw = nullptr;  // per-bin area weight",
               "  const int* isC = nullptr;    // per-bin grain type (1 "
               "carbon, 0 silicate)",
               "};",
               "// the per-cell state the solver sets for helpers used inside",
               "// the rates: n_global = the abundances at the start of the",
               "// call (legacy n_global), the photon field (phi, nbins bins",
               "// with edges elo/ehi [eV]), the binned-dust hook. Not a",
               "// parameter: set by solve() / the device backend / the",
               "// queries in their own copy of Params",
               "struct StateHooks {",
               "  const double* n_global = nullptr;",
               "  const double* phi = nullptr;",
               "  int nbins = 0;",
               "  const double* elo = nullptr;",
               "  const double* ehi = nullptr;",
               "  DustBinsHook dust;",
               "};",
               "// runtime parameters of this network (per cell); a NaN",
               "// default marks a REQUIRED parameter (no default exists),",
               "// rejected by the solver until it is set",
               "struct Params {"]
        for name, dflt in sorted(net.params.items()):
            out.append("  double %s = %s;" % (name, self.cnum(dflt)))
        for nm, dv, cm in univ:
            val = ("std::numeric_limits<double>::quiet_NaN()" if dv is None
                   else self.cnum(dv))
            extra = ("; REQUIRED by %s" % ", ".join(req_by.get(nm, ["?"]))
                     if dv is None else "")
            out.append("  double %s = %s;  // %s%s" % (nm, val, cm, extra))
        out.append("  StateHooks hk;  // solver-set per-cell state (not a "
                   "parameter)")
        out.append("};")
        out.append("constexpr int NPAR = %d;" % len(self.pnames))
        out.append("inline const char* param_names[NPAR] = {%s};"
                   % ", ".join('"%s"' % s for s in self.pnames))
        out.append("// what a REQUIRED parameter is needed for (\"\" = has "
                   "a default)")
        out.append("inline const char* param_required_by[NPAR] = {%s};"
                   % ", ".join('"%s"' % (", ".join(req_by.get(s, ["?"]))
                                         if s in self.param_required
                                         else "") for s in self.pnames))
        out.append("inline double* param_ptr(Params& p, int i) {")
        out.append("  switch (i) {")
        for i, s in enumerate(self.pnames):
            out.append("    case %d: return &p.%s;" % (i, s))
        out.append("  }")
        out.append("  return nullptr;")
        out.append("}")
        out.append("// a parameter row (param_names order: the C interface's "
                   "par[] rows) into a")
        out.append("// Params, on the host or the device (device-resident "
                   "calls build each")
        out.append("// cell's parameters from its row on the device)")
        out.append("KROME_HD inline void params_from_row(const double* row, "
                   "Params& p) {")
        out.append("  (void)row; (void)p;")
        for i, s_ in enumerate(self.pnames):
            out.append("  p.%s = row[%d];" % (s_, i))
        out.append("}")
        out.append("// index of the first REQUIRED parameter still unset "
                   "(NaN), or -1")
        out.append("KROME_HD inline int params_missing(const Params& p) {")
        out.append("  (void)p;")
        for i, s in enumerate(self.pnames):
            if s in self.param_required:
                out.append("  if (p.%s != p.%s) return %d;" % (s, s, i))
        out.append("  return -1;")
        out.append("}")
        for s in self.pnames:
            out.append("#define KROME_HAS_PAR_%s 1" % s)
        core = sec["core"]
        k = core.index(self.PARAMS_SLOT)
        core[k:k + 1] = out

    def thermal_option_name(self, flag):
        """ThermalConfig flag -> the option that enables it (section.key)"""
        if flag == "cooling_floor":
            return "thermal.floor"
        for opt, table in (("thermal.cooling", self.COOLING_NAMES),
                           ("thermal.heating", self.HEATING_NAMES)):
            for nm, fl in table.items():
                if fl == flag:
                    return "%s = ['%s']" % (opt, nm)
        return flag

    def resolve_rate_scale(self):
        """network `rate_scale = expr` (network-format.md 6.2): fold the factor
        into every expression rate, after auto/surface resolution, so the
        Jacobian, the piecewise groups and the GPU path see ordinary
        expressions. The unscaled AST is kept for selfcheck."""
        net = self.net
        if net.rate_scale is None or getattr(self, "_scale_done", False):
            return
        self._scale_done = True
        for rx in net.reactions:
            r = rx["rate"]
            if r["kind"] != "expr":
                continue
            r["ast_unscaled"] = r["ast"]
            r["ast"] = {"b": "*", "l": copy.deepcopy(net.rate_scale),
                        "r": r["ast"]}

    def emit_composition(self, w, sp_sorted):
        """element composition matrix and reaction labels (query API and
        conservation): elements in sorted order, opaque species all zero"""
        net = self.net
        elems = sorted({e for s in sp_sorted
                        for e in (net.species[s]["composition"] or {})})
        self.elements = elems
        w("// element composition: species_elements[i*NELEM + e] = atoms of")
        w("// element e in species i (opaque species: 0, species_opaque)")
        w("constexpr int NELEM = %d;" % len(elems))
        w("inline const char* element_names[%d] = {%s};"
          % (max(len(elems), 1), ", ".join('"%s"' % e for e in elems)
             if elems else '""'))
        rows = []
        for s in sp_sorted:
            comp = net.species[s]["composition"] or {}
            rows.append(", ".join(str(comp.get(e, 0)) for e in elems))
        w("inline const int species_elements[%d] = {%s};"
          % (max(len(elems) * len(sp_sorted), 1),
             ", ".join(rows) if elems else "0"))
        w("inline const bool species_opaque[NSPEC] = {%s};"
          % ", ".join("true" if net.species[s]["opaque"] else "false"
                      for s in sp_sorted))

        def label(rx):
            tr = rx.get("trange")
            s = rx["eq"]
            if tr:
                s += " [%s, %s]" % tuple(self.cnum(v) for v in tr)
            return s.replace("\\", "\\\\").replace('"', '\\"')
        w("// reaction labels (equation, plus the T range of a piece)")
        w("inline const char* reaction_names[%d] = {%s};"
          % (max(len(net.reactions), 1),
             ", ".join('"%s"' % label(rx) for rx in net.reactions)
             if net.reactions else '""'))

    def emit_num2col(self, w):
        """num2col_(): number -> column density (legacy num2col, the
        radiation.num2col method), used by the molecular coolants, the optical
        depth query and the Habing dust shielding"""
        w("// number -> column density (legacy num2col, method '%s')"
          % self.num2col)
        w("double inline num2col_(double x, double Tgas, "
          "const double* n, const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if self.num2col == "default":
            w("  return 1.87e21 * "
              "std::pow(std::fmax(x, 1e-40) * 1e-3, 2.0/3.0);")
        elif self.num2col == "length":
            w("  // the caller's length per cell (the parameter Lcol)")
            w("  return std::fmax(x, 1e-40) * p.Lcol;")
        else:
            if self.num2col == "jeans40":
                w("  const double T = std::fmin(std::fmax(Tgas, p.Tcmb),"
                  " 4e1);")
            else:
                w("  const double T = std::fmax(Tgas, p.Tcmb);")
            w("  // 0.5 * n * Jeans length (legacy get_jeans_length;")
            w("  // the constants pi, kB, p_mass, G)")
            w("  const double rhogas = std::fmax(rho_(n), 1e-40);")
            w("  const double mu = rhogas / std::fmax(ntot_(n), 1e-40)"
              " / constants::p_mass;")
            w("  const double lj = std::sqrt(constants::pi "
              "* constants::boltzmann_erg * T / rhogas")
            w("      / constants::p_mass / constants::gravity / mu);")
            w("  return 0.5 * std::fmax(x, 1e-40) * lj;")
        w("}")
        terms = []
        for sp_, f in (("H", ""), ("H+", ""), ("H2", "2.0 * ")):
            if sp_ in self.net.species:
                terms.append("%snum2col_(n[IDX_%s], Tgas, n, p)"
                             % (f, self.cid(sp_)))
        w("// H nuclei column N_H = N(H) + N(H+) + 2 N(H2) from num2col (the")
        w("// Habing dust shielding and the dust column of get_opacity)")
        w("inline double NH_col_(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        w("  return %s;" % (" + ".join(terms) if terms else "0.0"))
        w("}")
        w("// dust shielding of the photon-field Habing flux (legacy")
        w("// -shieldHabingDust, radiation.habing_dust_shielding): exp(-sigma_d N_H),")
        w("// sigma_d = 2e-21 (d2g/0.01) 0.665 cm^2 (Richings+14), N_H the")
        w("// num2col columns of H, H+ and 2 H2")
        w("inline double habing_shield_(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if self.habing_shield:
            if not any(x in self.net.species for x in ("H", "H+", "H2")):
                raise ValueError("radiation.habing_dust_shielding needs H, H+ or H2")
            w("  return std::exp(-2e-21 * (p.d2g / 1e-2) * 0.665 "
              "* NH_col_(Tgas, n, p));")
        else:
            w("  return 1.0;")
        w("}")
        w("")

    def emit_h2_extra(self, w):
        """the H2 cooling opacity h2_opacity_() (legacy -H2opacity) applied
        to the H2 cooling."""
        mode = self.h2_opacity
        w("// optically thick suppression of the H2 line cooling (legacy")
        w("// -H2opacity; thermal.h2_opacity=%s)" % (mode or "none"))
        if mode == "OMUKAI":
            xs, ys, zs = self.load_u2d("escape_H2.dat")
            w("// escape probability table, K. Omukai 2014 (data/escape_H2.dat)")
            w("// in log10 N (num2col of ntot, radiation.num2col) and log10 Tgas, the")
            w("// legacy uniform-grid bilinear fit with linear extrapolation")
            self.emit_data_array(w, "h2esc_x", xs, "log10 N [cm^-2]")
            self.emit_data_array(w, "h2esc_t", ys, "log10 Tgas [K]")
            self.emit_data_array(w, "h2esc_z", zs, "log10 escape fraction")
        w("inline double h2_opacity_(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if mode == "RIPAMONTI":
            w("  // Ripamonti & Abel 2004")
            w("  return std::fmin(1.0, std::pow(std::fmax(1.25e-10 * ntot_(n),"
              " 1e-40), -0.45));")
        elif mode == "OMUKAI":
            nx, ny = len(xs), len(ys)
            for ln in self._num2col_c("ntot_(n)"):
                w("  " + ln)
            w("  const double xx = std::log10(N), yy = std::log10(Tgas);")
            w("  const double xmul = %s, ymul = %s;"
              % (self.cnum(1.0 / (xs[1] - xs[0])),
                 self.cnum(1.0 / (ys[1] - ys[0]))))
            w("  int i = (int)((xx - h2esc_x[0]) * xmul);")
            w("  i = i < 0 ? 0 : (i > %d ? %d : i);" % (nx - 2, nx - 2))
            w("  int j = (int)((yy - h2esc_t[0]) * ymul);")
            w("  j = j < 0 ? 0 : (j > %d ? %d : j);" % (ny - 2, ny - 2))
            w("  const double pf = (xx - h2esc_x[i]) * xmul;")
            w("  const double zl = pf * (h2esc_z[(i+1)*%d + j] - "
              "h2esc_z[i*%d + j]) + h2esc_z[i*%d + j];" % (ny, ny, ny))
            w("  const double zr = pf * (h2esc_z[(i+1)*%d + j+1] - "
              "h2esc_z[i*%d + j+1]) + h2esc_z[i*%d + j+1];" % (ny, ny, ny))
            w("  return std::pow(10.0, (yy - h2esc_t[j]) * ymul * (zr - zl) "
              "+ zl);")
        else:
            w("  return 1.0;")
        w("}")
        w("")

    def _natoms(self, s):
        """atoms of a species as legacy counts them (the electron is one;
        H2 counts as a diatomic also when only its spin forms exist)"""
        if s == "H2" and s not in self.net.species:
            return 2
        sp = self.net.species[s]
        if sp["opaque"]:
            return 0
        if s == "e-":
            return 1
        return sum((sp["composition"] or {}).values())

    def emit_gamma(self, w):
        """gamma_index(): the adiabatic index of the thermal ODE, legacy
        -gamma (krome_gadiab.f90 + kromeobj makeGammaAdiabatic):
          CONST    the per-cell parameter gamma_ad (legacy DEFAULT: 5/3)
          EXACT    1 + ysum/gsum (Omukai & Nishi 98, Boley+07): monoatomic
                   species 3/2, each diatomic 1/2 (3 + rot + vib), exact
                   rotation (partition functions; H2 ortho/para) and
                   vibration; polyatomics and diatomics without constants
                   are skipped, as legacy
          POPOVAS  H2 and CO only, from the Popovas & Jorgensen (2016) ro-vibrational
                   partition function tables (data/partition), legacy
                   gamma_pop
        (legacy FULL, VIB, ROT and REDUCED were removed in 2.0.0:
        GAMMA_REMOVED)"""
        mode = self.gamma
        if mode == "POPOVAS":
            return self.emit_gamma_popovas(w)
        w("// adiabatic index of the thermal ODE (thermal.gamma=%s)" % mode)
        sp_sorted = sorted(self.net.species)
        if mode == "CONST":
            w("inline double gamma_index(double Tgas, const double* n, "
              "const Params& p) {")
            w("  (void)Tgas; (void)n;")
            w("  return p.gamma_ad;")
            w("}")
            w("")
            return
        mono = [s for s in sp_sorted if self._natoms(s) == 1]
        dis = []  # EXACT: (species, theta_vib, theta_rot), None if unknown
        # H2 in any form is one diatomic (density nH2_, its ortho fraction)
        spin = [x for x in ("H2_ortho", "H2_para") if x in self.net.species]
        if spin and "H2" not in sp_sorted:
            sp_sorted = sorted(sp_sorted + ["H2"])
        for s in sp_sorted:
            if s in spin:
                continue
            if self._natoms(s) != 2:
                continue
            ve, be = GAMMA_VIB_K.get(s), GAMMA_ROT_K.get(s)
            if ve is None and be is None:
                continue
            if ve is None:
                self.warnings.append("thermal.gamma=%s: no vibrational constant for "
                                     "%s" % (mode, s))
            if be is None:
                self.warnings.append("thermal.gamma=%s: no rotational constant for "
                                     "%s" % (mode, s))
            dis.append((s, ve, be))
        ves = [ve for _, ve, _ in dis if ve is not None]
        if any(be is not None for s, _, be in dis if s != "H2"):
            w("// rotational heat capacity C_rot/k = d/dT (T^2 dlnZ/dT) of the")
            w("// diatomic partition function (11 levels, T >= 10 K as legacy),")
            w("// EXACTLY as the variance of E_j/T: legacy took a numerical")
            w("// second difference (dT = 1e-5 T), whose ~1e-6 roundoff noise")
            w("// would corrupt the difference-quotient Jacobian")
            w("inline double gamma_rot_(double Tin, double brot) {")
            w("  const double ib = brot / std::fmax(Tin, 10.0);")
            w("  double z = 0.0, e1 = 0.0, e2 = 0.0;")
            w("  for (int j = 0; j <= 10; ++j) {")
            w("    const double E = j * (j + 1.0) * ib;")
            w("    const double g = (2.0 * j + 1.0) * std::exp(-E);")
            w("    z += g; e1 += g * E; e2 += g * E * E;")
            w("  }")
            w("  const double m = e1 / z;")
            w("  return e2 / z - m * m;")
            w("}")
        if any(s == "H2" and be is not None for s, _, be in dis):
            w("// H2 with an ortho:para ratio opr (legacy zfop: rotational")
            w("// constant 85.4 K, Z = Zp^b (Zo e^(2 theta/T))^a, a = opr/(opr+1)):")
            w("// C_rot/k = b Var_para(E/T) + a Var_ortho(E/T), exact")
            w("inline double gamma_rotop_(double Tin, double a) {")
            w("  const double ib = 85.4 / std::fmax(Tin, 10.0);")
            w("  double zp = 0.0, p1 = 0.0, p2 = 0.0, zo = 0.0, o1 = 0.0, "
              "o2 = 0.0;")
            w("  for (int j = 0; j <= 10; j += 2) {")
            w("    const int j1 = j + 1;")
            w("    const double Ep = j * (j + 1.0) * ib;")
            w("    const double Eo = j1 * (j1 + 1.0) * ib;")
            w("    const double gp = (2.0 * j + 1.0) * std::exp(-Ep);")
            w("    const double go = 3.0 * (2.0 * j1 + 1.0) * std::exp(-Eo);")
            w("    zp += gp; p1 += gp * Ep; p2 += gp * Ep * Ep;")
            w("    zo += go; o1 += go * Eo; o2 += go * Eo * Eo;")
            w("  }")
            w("  const double b = 1.0 - a;  // a: the ortho fraction")
            w("  const double mp = p1 / zp, mo = o1 / zo;")
            w("  return b * (p2 / zp - mp * mp) + a * (o2 / zo - mo * mo);")
            w("}")
        w("inline double gamma_index(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if ves:
            w("  // avoid overflowing exp(x = theta_vib / T) (legacy Tgas_vib)")
            w("  const double iTv = 1.0 / std::fmax(Tgas, %s);"
              % self.cnum(min(ves) * 1e-2))
        terms = []
        for s, ve, be in dis:
            cs = self.cid(s)
            vpart = "0.0"
            if ve is not None:
                w("  const double x_%s = %s * iTv, e_%s = std::exp(x_%s);"
                  % (cs, self.cnum(ve), cs, cs))
                vpart = ("2.0 * x_%s * x_%s * e_%s / ((e_%s - 1.0) * "
                         "(e_%s - 1.0))" % (cs, cs, cs, cs, cs))
            rpart = "2.0"
            if be is not None:
                rpart = ("2.0 * gamma_rotop_(Tgas, h2_ofrac_(n, p))"
                         if s == "H2"
                         else "2.0 * gamma_rot_(Tgas, %s)" % self.cnum(be))
                if s == "H2":
                    self._need_h2op = True
            w("  const double gi_%s = 0.5 * (3.0 + %s + %s);"
              % (cs, rpart, vpart))
            terms.append(cs)
        w("  const double mosum = %s;"
          % (" + ".join("n[IDX_%s]" % self.cid(s) for s in mono) or "0.0"))
        dens = {self.cid(s): self.gamma_dens(s) for s, _, _ in dis}
        w("  const double ysum = mosum%s;"
          % "".join(" + %s" % dens[c] for c in terms))
        w("  const double gsum = 1.5 * mosum%s;"
          % "".join(" + %s * gi_%s" % (dens[c], c) for c in terms))
        w("  return gsum > 0.0 ? 1.0 + ysum / gsum : 5.0 / 3.0;")
        w("}")
        w("")

    def load_partition(self, name):
        """data/partition/part<name>.dat: (Tmin, dT, z[]) on a uniform grid"""
        f = self.data_dir.parent / "partition" / ("part%s.dat" % name)
        T, z = [], []
        for line in f.read_text().splitlines():
            s = line.split("#")[0].split()
            if len(s) == 2:
                T.append(float(s[0]))
                z.append(float(s[1]))
        dT = T[1] - T[0]
        if any(abs((T[i + 1] - T[i]) - dT) > 1e-9 * dT
               for i in range(len(T) - 1)):
            raise ValueError("%s: non-uniform temperature grid" % f)
        return T[0], dT, z

    def emit_gamma_popovas(self, w):
        """thermal.gamma = 'POPOVAS' (legacy gamma_pop / gamma_pop_H2): monoatomic
        species 3/2, H2 and CO 1/2 (3 + 2 C) with C/k = d/dT (T^2 dlnZ/dT)
        of the tabulated partition functions (Popovas & Jorgensen 2016; T =
        10..20000 K every 10 K). The tables are labelled rotational but
        hold the ro-vibrational levels: CO has C = 1 at 100-300 K and ~2 at
        3000 K, and C falls again above ~5000 K (finite level set). No
        separate vibrational term, as legacy. As legacy:
        C_rot by forward differences at the two grid nodes around T,
        linearly interpolated; H2 as Z = Z_even^b (3 Z_odd e^(2 85.4/T))^a,
        a = opr/(opr+1). Fixed: legacy passed (min, dT) in swapped order
        (harmless, both 10 K), sized the CO table 641 of its 2000 values
        and read past the table end above ~6400 K (CO) / 19970 K; here T is
        clamped to the last node that has three neighbours."""
        sp_sorted = sorted(self.net.species)
        mono = [s for s in sp_sorted if self._natoms(s) == 1]
        dis = [s for s in ("H2", "CO") if s in self.net.species
               or (s == "H2" and self.h2_forms())]
        w("// Popovas & Jorgensen (2016) rotational partition functions: C_rot/k =")
        w("// d/dT (T^2 dlnZ/dT) by forward differences at the grid nodes")
        w("// around T, linearly interpolated (legacy gamma_pop)")
        tabs = {}
        for s in dis:
            for part in (("H2even", "H2odd") if s == "H2" else ("CO",)):
                t0, dT, z = self.load_partition(part)
                tabs[part] = (t0, dT, len(z))
                self.emit_data_array(w, "gpop_%s" % part, z,
                                     "rotational Z(T), T = %g + %g i K"
                                     % (t0, dT))
        for s in dis:
            part = "H2even" if s == "H2" else "CO"
            t0, dT, nz = tabs[part]
            args = "double Tin, double a" if s == "H2" else "double Tin"
            w("KROME_HD inline double gamma_pop_%s_(%s) {" % (s, args))
            if s == "H2":
                w("  const double b = 1.0 - a;  // a: the ortho fraction")
                w("  auto lnz = [&](int i) {")
                w("    const double Ti = %s + %s * i;" % (self.cnum(t0),
                                                          self.cnum(dT)))
                w("    return b * std::log(gpop_H2even[i]) + a * std::log(3.0 "
                  "* gpop_H2odd[i] * std::exp(2.0 * 85.4 / Ti));")
                w("  };")
            else:
                w("  (void)0;")
                w("  auto lnz = [&](int i) { return std::log(gpop_CO[i]); };")
            w("  const double T0 = %s, dT = %s;" % (self.cnum(t0),
                                                     self.cnum(dT)))
            w("  const double T = std::fmin(std::fmax(Tin, T0), T0 + dT * %d);"
              % (nz - 4))
            w("  const int i = static_cast<int>((T - T0) / dT);")
            w("  auto cv = [&](int k) {  // C_rot at node k")
            w("    const double Tk = T0 + dT * k;")
            w("    const double e1 = Tk * Tk * (lnz(k + 1) - lnz(k)) / dT;")
            w("    const double e2 = (Tk + dT) * (Tk + dT) * (lnz(k + 2) - "
              "lnz(k + 1)) / dT;")
            w("    return (e2 - e1) / dT;")
            w("  };")
            w("  const double c1 = cv(i), c2 = cv(i + 1);")
            w("  return c1 + (c2 - c1) * (T - (T0 + dT * i)) / dT;")
            w("}")
        w("inline double gamma_index(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        for s in dis:
            cs = self.cid(s)
            rp = ("gamma_pop_H2_(Tgas, h2_ofrac_(n, p))" if s == "H2"
                  else "gamma_pop_CO_(Tgas)")
            if s == "H2":
                self._need_h2op = True
            w("  const double gi_%s = 0.5 * (3.0 + 2.0 * %s);" % (cs, rp))
        w("  const double mosum = %s;"
          % (" + ".join("n[IDX_%s]" % self.cid(s) for s in mono) or "0.0"))
        w("  const double ysum = mosum%s;"
          % "".join(" + %s" % self.gamma_dens(s) for s in dis))
        w("  const double gsum = 1.5 * mosum%s;"
          % "".join(" + %s * gi_%s" % (self.gamma_dens(s), self.cid(s))
                    for s in dis))
        w("  return gsum > 0.0 ? 1.0 + ysum / gsum : 5.0 / 3.0;")
        w("}")
        w("")

    def emit_freeze(self, w):
        """network `freeze` directives (network-format.md 6.2) ->
        frozen_species(): the solver zeroes those rows of the full RHS
        (chemistry, secondary ionization, radiation coupling) and of the
        analytic Jacobian, as
        legacy @ODEmodifier dn(idx) = 0 did after the ODE evaluation."""
        net = self.net
        nmax = sum(len(f["species"]) for f in net.freeze)
        w("// freeze directives: species whose dn/dt is held at zero while")
        w("// the condition holds (the condition is constant during a call)")
        w("constexpr int NFREEZE_MAX = %d;" % nmax)
        w("KROME_HD inline int frozen_species(const Params& p, int* idx) {")
        w("  (void)p; (void)idx;")
        w("  int m = 0;")
        for f in net.freeze:
            rows = " ".join("idx[m++] = IDX_%s;" % self.cid(s)
                            for s in f["species"])
            if f["cond"] is None:
                w("  %s  // line %d" % (rows, f["line"]))
            else:
                w("  if (%s) {  // line %d" % (self.cexpr(f["cond"]),
                                               f["line"]))
                w("    " + rows)
                w("  }")
        w("  return m;")
        w("}")
        w("")

    # generated output sections, in dependency order; each can be written
    # to its own file (generate_files) or concatenated (generate).
    # Big embedded tables go into a separate 'data' bucket: in split mode
    # a net_<name>_data.cpp translation unit (compiled ONCE — including
    # the headers stays cheap), in single-file mode appended at the end
    # (single-TU use only).
    SECTIONS = ("core", "rates", "thermal", "jacobian")
    DATA_THRESHOLD = 512  # values; larger arrays go to the data bucket

    def emit_data_array(self, w, name, values, comment=""):
        """emit a double array: small ones inline where they are used,
        large ones as extern declarations backed by the data bucket."""
        lines, cur = [], []
        for v in values:
            cur.append(self.cnum(v))
            if len(cur) == 8:
                lines.append(", ".join(cur))
                cur = []
        if cur:
            lines.append(", ".join(cur))
        if len(values) <= self.DATA_THRESHOLD:
            w("inline const double %s[%d] = {%s" % (name, len(values),
                                                    ""))
            for i, ln in enumerate(lines):
                w("  %s%s" % (ln, "," if i < len(lines) - 1 else ""))
            w("};")
            return
        w("#ifndef KROME_TABLES_INLINE")
        w("extern const double %s[%d];%s" % (name, len(values),
                                             ("  // " + comment)
                                             if comment else ""))
        w("#endif")
        # the same table in accessor form, for builds that inline the big
        # tables (device builds: KROME_TABLES_INLINE, net_<name>_tables.hpp)
        t = self.__dict__.setdefault("_devtab", [])
        t.append("// %s%s" % (name, (": " + comment) if comment else ""))
        t.append("inline const double %s[%d] = {" % (name, len(values)))
        for i, ln in enumerate(lines):
            t.append("  %s%s" % (ln, "," if i < len(lines) - 1 else ""))
        t.append("};")
        d = self._data_w
        d("// %s%s" % (name, (": " + comment) if comment else ""))
        d("const double %s[%d] = {" % (name, len(values)))
        for i, ln in enumerate(lines):
            d("  %s%s" % (ln, "," if i < len(lines) - 1 else ""))
        d("};")

    ARRAY_RE = re.compile(
        r"^inline const (double|int|bool) (\w+)\[([^\]]+)\] = \{(.*?)\};"
        r"([^\n]*)$", re.M | re.S)

    @classmethod
    def devicify(cls, text):
        """constant arrays as device-safe accessors: every
        `inline const T name[N] = {...};` (T = double, int, bool) becomes a
        KROME_HD function holding a function-local static constexpr array
        (one instance per host / device compilation, which nvcc, clang
        CUDA and HIP all accept), plus a macro so that every `name[i]`
        access stays as written. Names (const char*) and the function
        pointer tables are host-only setup data and stay as they are."""
        def sub(m):
            t, name, n, init, tail = m.groups()
            return ("KROME_HD inline const %s* %s_arr_() {%s\n"
                    "  static constexpr %s t_[%s] = {%s};\n"
                    "  return t_;\n"
                    "}\n"
                    "#define %s (::krome::%s_arr_())"
                    % (t, name, tail, t, n, init, name, name))
        text = cls.ARRAY_RE.sub(sub, text)

        # device-callable functions: every function head at column 0
        # (definitions and the runtime-extern prototypes) gets KROME_HD,
        # except the host-only ones (HOST_ONLY_FUNCS)
        def fsub(m):
            if m.group(4) in cls.HOST_ONLY_FUNCS:
                return m.group(0)
            return "KROME_HD " + m.group(0)
        return cls.FUNC_RE.sub(fsub, text)

    FUNC_RE = re.compile(
        r"^((?:static |inline |constexpr )*)(double|void|int|bool|"
        r"const double\*|const int\*|const bool\*)((?: inline)?) "
        r"(\w+)\(", re.M)
    # host-only generated functions: the dense jacobian() wrapper (a
    # std::vector scratch; the host solver's entry point)
    HOST_ONLY_FUNCS = {"jacobian"}

    def tables_header(self):
        """net_<name>_tables.hpp: the big tables in accessor form, for the
        builds that inline them (KROME_TABLES_INLINE: device builds)"""
        body = getattr(self, "_devtab", [])
        return "\n".join(
            ["// generated by codegen.py from network '%s' (big tables, "
             "inline" % self.net.name,
             "// accessor form for KROME_TABLES_INLINE builds) -- do not "
             "edit", "#pragma once"] + KROME_HD_DEF
            + ["namespace krome {", ""]
            + [self.devicify("\n".join(body))] + ["", "} // namespace krome"]
        ) + "\n"

    def generate(self):
        """single self-contained header (net_<name>.hpp); the big-table
        definitions are appended, so include it from ONE translation unit
        (the split output has no such restriction)."""
        sec = self._generate_sections()
        out = ["// generated by codegen.py from network '%s' -- do not edit"
               % self.net.name,
               "#pragma once", "#include <cmath>", "#include <limits>",
               "#include <vector>", ""] + KROME_HD_DEF + [
               "namespace krome {",
               ""]
        if getattr(self, "_devtab", None):
            out += ["#ifdef KROME_TABLES_INLINE",
                    self.devicify("\n".join(self._devtab)), "#endif", ""]
        for s in self.SECTIONS:
            out.append(self.devicify("\n".join(sec[s])))
        if sec["data"]:
            out += ["", "// big embedded tables (single-file mode: defined",
                    "// here, so include this header from ONE translation",
                    "// unit; split mode puts these in net_<name>_data.cpp)",
                    "#ifndef KROME_TABLES_INLINE"]
            out += sec["data"]
            out.append("#endif")
        out.append("} // namespace krome")
        return "\n".join(out) + "\n"

    def generate_files(self, fortran=False):
        """split output: {filename: content}. net_<name>.hpp is the
        umbrella users include; each section file chain-includes its
        predecessor, so the big embedded tables (thermal) and the analytic
        Jacobian live in files of their own."""
        self.fortran = fortran
        sec = self._generate_sections()
        base = "net_%s" % self.net.name
        files = {}
        prev = None
        for s in self.SECTIONS:
            out = ["// generated by codegen.py from network '%s' "
                   "(section: %s) -- do not edit" % (self.net.name, s),
                   "#pragma once", "#include <cmath>", "#include <limits>",
                   "#include <vector>"] + KROME_HD_DEF
            if prev:
                out.append('#include "%s_%s.hpp"' % (base, prev))
            elif getattr(self, "_devtab", None):
                # device builds inline the big tables (accessor form)
                out += ["#ifdef KROME_TABLES_INLINE",
                        '#include "%s_tables.hpp"' % base, "#endif"]
            out += ["", "namespace krome {", "",
                    self.devicify("\n".join(sec[s])),
                    "", "} // namespace krome"]
            files["%s_%s.hpp" % (base, s)] = "\n".join(out) + "\n"
            prev = s
        # big embedded tables: their own translation unit, compiled once
        # (the headers carry only extern declarations); always emitted so
        # build scripts can link it unconditionally
        files[base + "_data.cpp"] = "\n".join(
            ["// generated by codegen.py from network '%s' (big embedded"
             % self.net.name,
             "// tables) -- do not edit; compile and link this once",
             "// (krome.hpp for the tables' extern declarations and for the",
             "// runtime functions the network headers declare: a unit that",
             "// sees their declarations without the definitions is an",
             "// error with GCC, 'inline function used but never defined')",
             '#include "krome.hpp"', "",
             "// device builds (KROME_TABLES_INLINE) take the tables from",
             "// %s_tables.hpp: this unit is then empty, so a build can" % base,
             "// always compile and link it",
             "#ifndef KROME_TABLES_INLINE", "namespace krome {", ""]
            + (sec["data"] or ["// (no tables above the size threshold)"])
            + ["", "} // namespace krome", "#endif"]) + "\n"
        if getattr(self, "_devtab", None):
            files[base + "_tables.hpp"] = self.tables_header()
        files.update(self.interface_files(fortran))
        files[base + ".hpp"] = "\n".join([
            "// generated by codegen.py from network '%s' -- do not edit"
            % self.net.name,
            "// umbrella header; the sections chain-include each other:",
            "//   %s_core.hpp      indices, species data, Params, helpers"
            % base,
            "//   %s_rates.hpp     tables, vars, photo, rates(), rhs()"
            % base,
            "//   %s_thermal.hpp   generated cooling/heating (+ tables)"
            % base,
            "//   %s_jacobian.hpp  rhs_scale(), analytic jacobian()" % base,
            "//   %s_data.cpp      big embedded tables (compile + link "
            "once)" % base,
            "#pragma once",
            '#include "%s_%s.hpp"' % (base, self.SECTIONS[-1])]) + "\n"
        return files

    def network_md(self, options):
        """NETWORK.md: what this project was generated with and the tables
        a caller needs (species, parameters, reactions). options: the
        options that differ from the defaults, "section.key = value"."""
        net = self.net
        sp = sorted(net.species)
        univ = {nm: (dv, cm) for nm, dv, cm in UNIVERSAL_PARAMS}

        def cell(x):
            return str(x).replace("|", "\\|")
        out = ["# Network %s" % net.name, "",
               "Generated by KROME-X %s codegen (do not edit: "
               "regeneration overwrites this file)." % krome_version(), "",
               "## Generation options", "",
               "The physics below is fixed in the generated code; at "
               "runtime only the solver numerics and `evolve_T` change. "
               "Regenerate with `python3 -m kromex krome.toml` (every option is "
               "recorded in krome.toml, in this directory).", "", "```"]
        out += options or ["(defaults)"]
        out += ["```", "",
                "## Species (NSPEC = %d)" % len(sp), "",
                "Alphabetical (code-point) order. Use the constants, "
                "never the numbers: C/C++ `IDX_<c>` (0-based), C "
                "`KROME_IDX_<c>` (0-based), Fortran `krome_idx_<c>` "
                "(1-based), Python `net.idx(name)`.", "",
                "| index | Fortran | species | constant `<c>` | "
                "mass [amu] | charge |",
                "|---:|---:|---|---|---:|---:|"]
        for i, s in enumerate(sp):
            d = net.species[s]
            out.append("| %d | %d | %s | %s | %.6g | %d |" % (
                i, i + 1, cell(s), self.cid(s), d["mass_amu"], d["charge"]))
        req_by = dict(zip(self.pnames, [""] * len(self.pnames)))
        req_by.update(self.param_required_by)
        out += ["", "## Parameters (NPAR = %d)" % len(self.pnames), "",
                "Per-cell runtime parameters: C++ `Params` fields "
                "(`p.<name>`), C/Fortran `par[]` rows in this order "
                "(`KROME_PAR_<name>` / `krome_par_<name>`), Python "
                "`par={\"<name>\": ...}`. REQUIRED parameters have no "
                "default: the solver rejects the cell (flag -105) until "
                "they are set.", "",
                "| index | name | default | meaning |",
                "|---:|---|---|---|"]
        for i, nm in enumerate(self.pnames):
            if nm in net.params:
                dflt, what = self.cnum(net.params[nm]), "network parameter"
            else:
                dv, what = univ[nm]
                dflt = ("**REQUIRED** by %s" % req_by[nm] if dv is None
                        else self.cnum(dv))
            out.append("| %d | %s | %s | %s |" % (i, nm, cell(dflt),
                                                   cell(what)))
        out += ["", "## Reactions (NREA = %d)" % len(net.reactions), "",
                "The order of `get_rates` / `get_fluxes` (0-based).", "",
                "| index | reaction |", "|---:|---|"]
        for i, rx in enumerate(net.reactions):
            lab = rx["eq"]
            if rx.get("trange"):
                lab += " [%s, %s]" % tuple(self.cnum(v)
                                           for v in rx["trange"])
            out.append("| %d | %s |" % (i, cell(lab)))
        return "\n".join(out) + "\n"

    def interface_files(self, fortran=False):
        """network-specific pieces of the library interface (runtime/interfaces/C/krome.h,
        runtime/interfaces/C/krome.cpp; optional bindings in runtime/interfaces/fortran,
        runtime/interfaces/python):
          krome_network.hpp  selects this network for krome.cpp
          krome_net.h        C constants: KROME_NSPEC, KROME_IDX_<sp>
                             (0-based), KROME_PAR_<param> (0-based)
          krome_net.f90      only with fortran (codegen project.bindings = ['fortran']): module
                             krome_net: krome_nspec, krome_idx_<sp>,
                             krome_par_<param> (1-based)
        The Python binding needs no generated file: it reads the network
        description from the library."""
        net = self.net
        sp_sorted = sorted(net.species)
        nelem = len({e for x in sp_sorted
                     for e in (net.species[x]["composition"] or {})})
        pnames = self.pnames
        head = ("generated by codegen.py from network '%s' -- do not edit"
                % net.name)
        out = {}
        out["krome_network.hpp"] = "\n".join([
            "// " + head,
            "// this project's network, included by krome.hpp (the C++ API) and",
            "// compiled into the C interface (krome.cpp)",
            "#pragma once",
            '#include "net_%s.hpp"' % net.name]) + "\n"

        h = ["/* " + head + " */",
             "#ifndef KROME_NET_H", "#define KROME_NET_H",
             "#define KROME_NETWORK \"%s\"" % net.name,
             "#define KROME_NSPEC %d" % len(sp_sorted),
             "#define KROME_NREA %d" % len(net.reactions),
             "#define KROME_NPAR %d" % len(pnames),
             "#define KROME_NELEM %d  /* elements (krome_element_name)"
             " */" % nelem,
             "/* species indices into n[cell][KROME_NSPEC] (0-based) */"]
        for i, s in enumerate(sp_sorted):
            h.append("#define KROME_IDX_%s %d  /* %s */"
                     % (self.cid(s), i, s))
        h.append("/* parameter indices into par[cell][KROME_NPAR] "
                 "(0-based) */")
        for i, s in enumerate(pnames):
            h.append("#define KROME_PAR_%s %d" % (s, i))
        h.append("#endif")
        out["krome_net.h"] = "\n".join(h) + "\n"
        if not fortran:
            return out

        # Fortran identifiers are case-insensitive: skip colliding names
        # (e.g. CO vs Co) instead of emitting an uncompilable module
        def fortran_safe(names):
            seen = {}
            for nm in names:
                seen.setdefault(nm.lower(), []).append(nm)
            return {nm for group in seen.values() if len(group) == 1
                    for nm in group}
        sp_ids = ["krome_idx_" + self.cid(s) for s in sp_sorted]
        par_ids = ["krome_par_" + s for s in pnames]
        ok = fortran_safe(sp_ids + par_ids)
        f = ["! " + head,
             "! network constants for the krome module "
             "(runtime/interfaces/fortran/krome.f90);",
             "! all indices are 1-BASED: n(krome_idx_H, cell),",
             "! par(krome_par_crate, cell)",
             "module krome_net",
             "  implicit none",
             "  character(len=*), parameter :: krome_network = \"%s\""
             % net.name,
             "  integer, parameter :: krome_nspec = %d" % len(sp_sorted),
             "  integer, parameter :: krome_nrea = %d" % len(net.reactions),
             "  integer, parameter :: krome_npar = %d" % len(pnames),
             "  integer, parameter :: krome_nelem = %d" % nelem]
        for i, (s, ident) in enumerate(zip(sp_sorted, sp_ids)):
            if ident in ok:
                f.append("  integer, parameter :: %s = %d  ! %s"
                         % (ident, i + 1, s))
            else:
                f.append("  ! %s: name clashes case-insensitively, use "
                         "krome_species_idx(\"%s\")" % (s, s))
                self.warnings.append("Fortran constant for species %s "
                                     "skipped (case-insensitive clash)" % s)
        for i, (s, ident) in enumerate(zip(pnames, par_ids)):
            if ident in ok:
                f.append("  integer, parameter :: %s = %d" % (ident, i + 1))
            else:
                f.append("  ! param %s: name clashes case-insensitively, "
                         "use krome_param_idx(\"%s\")" % (s, s))
        f.append("end module krome_net")
        out["krome_net.f90"] = "\n".join(f) + "\n"
        return out

    def _emit_unrolled_fluxes_rhs(self, w):
        """fluxes() (the reaction rates k * prod(n) [cm^-3 s^-1]) and rhs(),
        one statement per reaction"""
        net = self.net
        w("void inline fluxes(const double* n, const double* k, "
          "double* fl) {")
        w("  (void)n; (void)k; (void)fl;")
        for rx in net.reactions:
            i = rx["id"] - 1
            fac = ["k[%d]" % i] + [self.nden(r["sp"])
                                   for r in rx["reactants"] if r["kinetic"]]
            w("  fl[%d] = %s;" % (i, "*".join(fac)))
        w("}")
        w("")

        # rhs()
        w("void inline rhs(const double* n, const double* k, double* dn) {")
        w("  for (int i = 0; i < NSPEC; ++i) dn[i] = 0.0;")
        w("  double flux;")
        w("  (void)n; (void)k; (void)flux;  // an empty network")
        for rx in net.reactions:
            i = rx["id"] - 1
            fac = ["k[%d]" % i] + [self.nden(r["sp"])
                                   for r in rx["reactants"] if r["kinetic"]]
            w("  flux = %s;  // %s" % ("*".join(fac), rx["eq"]))
            if self.gas_total:
                # net increments: a gas member also updates its total
                # slot; exact freezeout/evaporation cancellations on the
                # total are dropped instead of emitted as -=/+= noise
                for sp, v in self.stoich_net(rx):
                    op = "+=" if v > 0 else "-="
                    mag = "" if abs(v) == 1 else "%d.0*" % abs(v)
                    w("  dn[IDX_%s] %s %sflux;" % (self.cid(sp), op, mag))
            else:
                for r in rx["reactants"]:
                    w("  dn[IDX_%s] -= flux;" % self.cid(r["sp"]))
                for pr in rx["products"]:
                    w("  dn[IDX_%s] += flux;" % self.cid(pr["sp"]))
        w("}")

    def _generate_sections(self):
        """the generated C++, as lines per section (SECTIONS, each a file of
        the split output, net_<name>_<section>.hpp, chained in this order)
        plus the big tables (data): the network is resolved first (auto,
        nucleation, reverse, surface, ...), then the sections are written in
        turn, `w` appending to the current one:
          core      species and parameter tables, the summing helpers, the
                    photo / LW / dust declarations and the externs
          rates     the rate functions, rates_fixed / rates_var, fluxes()
                    and rhs()
          thermal   the network's own heating and cooling terms (custom,
                    nuclear, cosmic rays, chemical, free-free, ...)
          jacobian  the analytic Jacobian (jacobian_split), its scaling and
                    the dk/dn difference quotients"""
        net = self.net
        self.resolve_auto()
        self.resolve_nucleation()
        self.resolve_multiplicity()
        self.resolve_reverse()
        self.resolve_screening()
        self.expand_surface_bins()
        self.resolve_surface()
        self.resolve_rate_scale()
        sp_sorted = sorted(net.species)
        nspec, nrea = len(sp_sorted), len(net.reactions)
        sec = {s: [] for s in self.SECTIONS}
        sec["data"] = []
        self._data_w = sec["data"].append

        # ==== core section ===================================================
        w = sec["core"].append

        w('inline const char* NETWORK_NAME = "%s";' % net.name)
        w("// the KROME-X version that generated this code (codegen")
        w("// --version; + the git revision of the tree)")
        w('inline const char* KROME_VERSION = "%s";' % krome_version())
        for ln in constants.cpp_block():
            w(ln)
        w("// the project's default solver tolerances (codegen solver.rtol /")
        w("// solver.atol): SolverOpts::rtol and atol_n start from them")
        w("constexpr double DEFAULT_RTOL = %s;" % self.cnum(self.rtol))
        w("constexpr double DEFAULT_ATOL_N = %s;" % self.cnum(self.atol))
        w("// the Fortran binding was generated (project.bindings = ['fortran']): the error")
        w("// report then also writes the failing call in Fortran")
        w("constexpr bool FORTRAN_BINDING = %s;"
          % ("true" if getattr(self, "fortran", False) else "false"))
        w("constexpr int NSPEC = %d;" % nspec)
        w("constexpr int NREA  = %d;" % nrea)
        for i, s in enumerate(sp_sorted):
            w("constexpr int IDX_%s = %d;" % (self.cid(s), i))
        w("// the per-species default tolerances (codegen solver.atol_file /")
        w("// solver.rtol_file): SolverOpts::atol_species / rtol_species start")
        w("// from them; -1 takes the scalar, nullptr when not given")
        for what, tab in (("ATOL", self.custom_atol),
                          ("RTOL", self.custom_rtol)):
            if not tab:
                w("constexpr const double* DEFAULT_%s_SPECIES_PTR = nullptr;"
                  % what)
                continue
            w("inline constexpr double DEFAULT_%s_SPECIES[NSPEC] = {" % what)
            for s in sp_sorted:
                w("    %s,  // %s" % (self.cnum(tab.get(s, -1.0)), s))
            w("};")
            w("constexpr const double* DEFAULT_%s_SPECIES_PTR ="
              " DEFAULT_%s_SPECIES;" % (what, what))
        w("")
        w("// a network is included (krome.hpp then does not include")
        w("// krome_network.hpp)")
        w("#define KROME_NETWORK_INCLUDED 1")
        w("// species-presence flags for runtime modules (the analog of the")
        w("// legacy #IFKROME_hasX pragmas)")
        for s in sp_sorted:
            w("#define KROME_HAS_%s 1" % self.cid(s))
        h2 = self.h2_forms()
        if h2:
            w("// total H2 density (H2, or its nuclear-spin species): the")
            w("// processes that need n(H2) read it through nH2_")
            w("#define KROME_HAS_H2_ANY 1")
            w("inline double nH2_(const double* n) { return %s; }"
              % " + ".join("n[IDX_%s]" % self.cid(x) for x in h2))
        w("")
        w("inline const char* species_names[NSPEC] = {%s};"
          % ", ".join('"%s"' % s for s in sp_sorted))
        w("inline const double species_mass_amu[NSPEC] = {%s};"
          % ", ".join(self.cnum(net.species[s]["mass_amu"])
                      for s in sp_sorted))
        w("inline const int species_charge[NSPEC] = {%s};"
          % ", ".join(str(net.species[s]["charge"]) for s in sp_sorted))
        w("inline const double species_mass_g[NSPEC] = {%s};"
          % ", ".join(self.cnum(net.species[s]["mass_amu"] * constants.amu)
                      for s in sp_sorted))
        self.emit_composition(w, sp_sorted)
        ebind = self.load_ebind()
        w("// binding energies [K] on bare grain / water ice "
          "(data/Ebare_ice.dat); 1e99 = unknown (evaporation -> 0)")
        w("inline const double species_Ebind_bare[NSPEC] = {%s};"
          % ", ".join(self.cnum(ebind.get(s, (1e99, 1e99))[0])
                      for s in sp_sorted))
        w("inline const double species_Ebind_ice[NSPEC] = {%s};"
          % ", ".join(self.cnum(ebind.get(s, (1e99, 1e99))[1])
                      for s in sp_sorted))
        w("")
        w(self.PARAMS_SLOT)  # filled at the end: needs the generated code
        w(self.H2OP_SLOT)    # ortho/para H2 helpers, when used
        w("")
        w("inline const double species_natoms_H[NSPEC] = {%s};"
          % ", ".join(self.cnum((net.species[s]["composition"] or {})
                                .get("H", 0)) for s in sp_sorted))
        # gas+total abundance representation (dust.gas_total;
        # docs/solver-design.md sect. 4): the surface slot of each
        # gas/surface pair stores the TOTAL (gas + surface) abundance, so
        # the summing helpers below weight the gas member 0 (state_weight)
        # and nothing is counted twice. Without dust.gas_total, or without
        # pairs, every weight is 1 and no state_weight table is emitted.
        pairs = []
        for s in sp_sorted:
            sp = net.species[s]
            if sp["phase"] != "gas":
                base = s.rsplit("_", 1)[0]
                if base in net.species:
                    pairs.append((base, s))
        if pairs:
            w("// gas<->surface pairs (for the gas+total representation): %s"
              % ", ".join("%s/%s" % p for p in pairs))
        weighted = self.gas_total and bool(self.pair_gas)
        if weighted:
            w("// GAS+TOTAL representation: the surface slot of each pair")
            w("// stores TOTAL = gas + surface, so summing over all slots")
            w("// must skip the gas member (weight 0) to avoid double")
            w("// counting; surface density = total slot - gas slot.")
            w("#define KROME_GAS_TOTAL 1")
            gaset = set(self.total_of)
            w("inline const double state_weight[NSPEC] = {%s};"
              % ", ".join("0.0" if s in gaset else "1.0"
                          for s in sp_sorted))
            w("static inline double ntot_(const double* n) {")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i] * "
              "state_weight[i];")
            w("  return s;")
            w("}")
            w("static inline double Hnuclei_(const double* n) {")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i] * "
              "state_weight[i] * species_natoms_H[i];")
            w("  return s;")
            w("}")
            w("static inline double rho_(const double* n) {  // g cm^-3")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i] * "
              "state_weight[i] * species_mass_g[i];")
            w("  return s;")
            w("}")
        else:
            w("// summing helpers; per-species weights become non-trivial "
              "in")
            w("// the gas+total abundance representation")
            w("static inline double ntot_(const double* n) {")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i];")
            w("  return s;")
            w("}")
            w("static inline double Hnuclei_(const double* n) {")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i] * "
              "species_natoms_H[i];")
            w("  return s;")
            w("}")
            w("static inline double rho_(const double* n) {  // g cm^-3")
            w("  double s = 0.0;")
            w("#pragma omp simd reduction(+:s)")
            w("  for (int i = 0; i < NSPEC; ++i) s += n[i] * "
              "species_mass_g[i];")
            w("  return s;")
            w("}")
        w("// the aggregates of the state the rates may depend on (0 ntot, 1")
        w("// Hnuclei, 2 rho): the Jacobian keeps their coupling as the")
        w("// low-rank part U W^T (jacobian_split), W_aj = agg_w_(a, j) =")
        w("// d aggregate_a / d y_j, the weights of the summing helpers")
        w("constexpr int NAGG = 3;")
        w("static inline double agg_w_(int a, int j) {")
        sw = " * state_weight[j]" if weighted else ""
        w("  return (a == 0 ? 1.0 : a == 1 ? species_natoms_H[j] "
          ": species_mass_g[j])%s;" % sw)
        w("}")
        w("// the gas-phase key elements of the grain types (dust.evolve:")
        w("// accretion and return; -1 if not in the network)")
        w("constexpr int DUST_KEY_C = %s;" % ("IDX_C" if "C" in sp_sorted
                                             else "-1"))
        w("constexpr int DUST_KEY_SI = %s;" % ("IDX_Si" if "Si" in sp_sorted
                                              else "-1"))
        # legacy krome_scale_Z: the neutral atomic metals, Asplund+09
        try:
            solar = self.load_solar()
        except FileNotFoundError:
            solar = {}
        metals = [s for s in sp_sorted if s in solar
                  and s not in ("H", "He")]
        w("// legacy krome_scale_Z: every neutral atomic metal of the network")
        w("// (a species named as an element, not H or He) set to")
        w("// n_H 10^(Z + log10(n/n_H)_sun), the Asplund+09 photospheric")
        w("// abundances (data/asplund.dat), Z = log10(Z/Z_sun), floored at")
        w("// 1e-40; nothing without H nuclei. Metals here: %s"
          % (", ".join(metals) or "none"))
        w("inline void scale_Z(double* n, double Z) {")
        if metals and "H" in sp_sorted:
            w("  const double Htot = Hnuclei_(n);")
            for m in metals:
                w("  n[IDX_%s] = std::fmax(Htot * std::pow(10.0, Z + (%s)), "
                  "1e-40);" % (self.cid(m), self.cnum(solar[m])))
        else:
            w("  (void)n; (void)Z;")
        w("}")
        w("")
        w("// underflow-safe d(l^r) = l^r (dr log l + r dl/l) for the")
        w("// analytic Jacobian: returns the correct 0 limit when the base")
        w("// underflows to 0 (e.g. exp(-E/T) at low Tgas) or is")
        w("// transiently non-positive, where the naive form is NaN")
        w("static inline double dpow_(double l, double r, double dl, "
          "double dr) {")
        w("  const double plr = std::pow(l, r);")
        w("  if (l <= 0.0 || plr == 0.0 || !std::isfinite(plr)) "
          "return 0.0;")
        w("  return plr * (dr * std::log(l) + r * dl / l);")
        w("}")
        w("// x^N for an integer literal exponent: products unrolled at")
        w("// compile time (binary powering), instead of std::pow")
        w("template <int N>")
        w("static inline double ipow_(double x) {")
        w("  if constexpr (N < 0) {")
        w("    return 1.0 / ipow_<-N>(x);")
        w("  } else if constexpr (N == 0) {")
        w("    (void)x;")
        w("    return 1.0;")
        w("  } else if constexpr (N == 1) {")
        w("    return x;")
        w("  } else {")
        w("    const double h = ipow_<N / 2>(x);")
        w("    if constexpr (N % 2 == 1) return h * h * x;")
        w("    else return h * h;")
        w("  }")
        w("}")
        w("")
        self.emit_thermal_config(w)
        w("// per-cell state the solver sets for the helpers used inside the")
        w("// rates (legacy n_global, the photobin field, the binned-dust")
        w("// hook) lives in Params::hk (StateHooks, defined with Params),")
        w("// set in the solver's own copy of the parameters: per cell and")
        w("// device-safe")
        w("// mean intensity J_nu [eV cm^-2 s^-1 Hz^-1 sr^-1] of photon bin b,")
        w("// constant across the bin as legacy photoBinJ: the photon flux")
        w("// Phi_b = 4 pi J_b ln(Ehi/Elo) / h")
        w("KROME_HD inline double photo_Jbin_(const Params& p, int b) {")
        w("  return p.hk.phi[b] * constants::planck_eV")
        w("      / (4.0 * constants::pi * std::log(p.hk.ehi[b] "
          "/ p.hk.elo[b]));")
        w("}")
        w("// ... at E [eV]: the bin containing E (legacy get_photoIntensity);")
        w("// 0 outside the bins or without a field")
        w("KROME_HD inline double photo_J_(const Params& p, double E) {")
        w("  for (int b = 0; b < p.hk.nbins; ++b)")
        w("    if (p.hk.elo[b] <= E && E <= p.hk.ehi[b])")
        w("      return photo_Jbin_(p, b);")
        w("  return 0.0;")
        w("}")
        w("// Lyman-Werner band (11.2-13.6 eV) mean intensity of the bins,")
        w("// J_LW = sum_b J_b ln(overlap_b) / ln(13.6/11.2): each bin")
        w("// constant in J (legacy photoBinJ); a flat J gives J_LW = J")
        w("constexpr double LW_EMIN = %s, LW_EMAX = %s;"
          % tuple(self.cnum(x) for x in knet.LW_BAND))
        w("KROME_HD inline double photo_J_LW_(const Params& p) {")
        w("  double s = 0.0;")
        w("  for (int b = 0; b < p.hk.nbins; ++b) {")
        w("    const double a = std::fmax(p.hk.elo[b], LW_EMIN);")
        w("    const double c = std::fmin(p.hk.ehi[b], LW_EMAX);")
        w("    if (c > a) s += photo_Jbin_(p, b) * std::log(c / a);")
        w("  }")
        w("  return s / std::log(LW_EMAX / LW_EMIN);")
        w("}")
        w("")
        w("// area-weighted mean Tdust of the binned-dust hook")
        w("KROME_HD static inline double dustbins_Td_mean_(const Params& p) {")
        w("  const DustBinsHook& h = p.hk.dust;")
        w("  double s = 0.0, w = 0.0;")
        w("  for (int b = 0; b < h.nb; ++b) {")
        w("    s += h.aw[b] * h.Td[b];")
        w("    w += h.aw[b];")
        w("  }")
        w("  return s / w;")
        w("}")
        w("")
        # extension functions: known ones get prototypes (implemented by the
        # runtime modules, e.g. runtime/core/dust_rates.hpp, included after this
        # header); unknown ones are stubbed to 0 and their rates are wrong
        # by construction (excluded from validation)
        for name in sorted(net.externs):
            if name in KNOWN_EXTERNS:
                w("#define KROME_EXTERN_%s 1" % name)
                w("inline " + KNOWN_EXTERNS[name])
            elif name in PLAIN_EXTERNS:
                w("#define KROME_EXTERN_%s 1" % name)
                w("inline " + PLAIN_EXTERNS[name])
            elif name in GENERATED_EXTERNS or name in self._gen_ext:
                pass  # defined in the rates section (emit_generated_externs)
            elif name in ("sigma_v96", "sigma_vy95") and any(
                    rx["rate"]["kind"] == "photo" for rx in net.reactions):
                pass  # the Verner fits are generated with the photo tables
            else:
                self.warnings.append("extern %s() stubbed to 0" % name)
                w("template <class... A> static inline double %s(A...) "
                  "{ return 0.0; }  // STUB" % name)

        # the Tdust state symbol routes through table_Tdust (defined in
        # the thermal section): parameter, table, or binned-dust hook
        # a template: the Params member it may read (Tdust) exists only
        # when something consumes the dust temperature
        w("template <class P>")
        w("double table_Tdust(double Tgas, const double* n, const P& p);")
        w("")

        # ==== rates section ==================================================
        w = sec["rates"].append
        # Tdust as state symbol: the p.Tdust parameter, or the live
        # table_Tdust(Tgas, n) when --dust-tables is generated (undef'd
        # per section so the name never leaks into user code)
        w(self._tdust_define())
        w("")
        self.emit_tables(w)
        # legacy rate-database helper functions used by resolved auto rates
        for name in sorted(getattr(self, "_tfuncs", ())):
            w(AUTO_TFUNCS[name])
        # num2col_ first: the generated shield_dust / num2col use it
        self.emit_num2col(w)
        self.emit_generated_externs(w)
        self.emit_reverse(w)
        self.emit_nucleation(w)
        self.emit_surface(w)

        # vars as inline functions (evaluated at the caller's T: this is what
        # makes per-group clamped temperatures work)
        for name in net.var_order:
            w("static inline double v_%s(double Tgas, const double* n, "
              "const Params& p) {" % name)
            w("  (void)Tgas; (void)n; (void)p;")
            w("  return %s;" % self.cexpr(net.vars[name]))
            w("}")
        w("")

        # rate expressions are emitted directly inside rates(); a callable
        # rate_<i>() is generated ONLY for extern/table-bearing rates, whose
        # dk/dn Jacobian terms need to re-evaluate the rate with perturbed
        # densities (jac_dk_lines)
        stubs = []
        photo = []  # (index-in-photo-list order) dicts
        for rx in net.reactions:
            i = rx["id"] - 1
            kind = rx["rate"]["kind"]
            if kind == "expr":
                if self.has_extern(rx["rate"]["ast"]):
                    # Tgas: the (group-clamped) rate temperature; Tg_: the
                    # actual gas temperature (1/K_c of reverse rates)
                    w("// %s  (line %d)" % (rx["eq"], rx["line"]))
                    w("static inline double rate_%d(double Tgas, "
                      "const double* n, const Params& p, double Tg_) {" % i)
                    w("  (void)Tgas; (void)n; (void)p; (void)Tg_;")
                    self._in_ratefn = True
                    w("  return %s;" % self.cexpr(rx["rate"]["ast"]))
                    self._in_ratefn = False
                    w("}")
            elif kind == "photo" and rx["rate"]["xsec"]["kind"] == "h2pd":
                # one photo channel per ground vibrational level j: its
                # lines as bin cross sections, weighted by the level
                # population z_j(T) as the channel's shielding factor
                for j, (lo, hi) in enumerate(self.h2pd_levels()):
                    if lo is None:
                        continue
                    photo.append({
                        "i": i, "ch": len(photo), "rx": rx, "code": "0.0",
                        "reactant": "H2", "Emin": lo, "Emax": hi,
                        "binned": j, "shield": {
                            "f": "h2pd_pop", "kind": "extern",
                            "a": [{"n": float(j)}]}})
                self._gen_ext.add("h2pd_pop")
                net.externs.add("h2pd_pop")
            elif kind == "photo":
                ph = self.resolve_photo(rx)
                ph["ch"] = len(photo)
                ph["shield"] = rx["rate"].get("shield")
                photo.append(ph)
                if ph["code"] is None and "table" not in ph:
                    stubs.append((i, rx))
            else:
                stubs.append((i, rx))
        if stubs:
            self.warnings.append("%d unsupported rates stubbed to 0 (%s)"
                                 % (len(stubs), ", ".join(sorted(
                                     {r["rate"]["kind"] for _, r in stubs}))))
        w("")
        self.emit_photo(w, photo)
        self.emit_rad_channels(w, photo)
        self.emit_secion(w)
        w("constexpr int N_UNSUPPORTED = %d;" % len(stubs))
        w("inline const int unsupported_rates[%d] = {%s};"
          % (max(len(stubs), 1),
             ", ".join(str(i) for i, _ in stubs) if stubs else "-1"))
        w("")

        # rates(): piecewise groups with clamp (or zero) semantics, split
        # into the rates that depend on Tgas and the parameters only
        # (rates_fixed: constant during a fixed-Tgas solve, computed once
        # per solve there) and the rest (rates_var: densities, externs,
        # tables, the live Tdust). Every k is written by exactly one of
        # them, so rates() is the same code as a single function.
        def rate_body(rx, tv):
            i = rx["id"] - 1
            if self.has_extern(rx["rate"]["ast"]):
                expr = "rate_%d(%s, n, p, Tgas)" % (i, tv)
            else:
                expr = self.cexpr(rx["rate"]["ast"], tv)
            return ["k[%d] = %s;  // %s" % (i, expr, rx["eq"])]

        for part, sel in (("fixed", self.rate_is_fixed),
                          ("var", lambda rx: not self.rate_is_fixed(rx))):
            w("void inline rates_%s(double Tgas, const double* n, "
              "const Params& p, double* k) {" % part)
            w("  (void)Tgas; (void)n; (void)p; (void)k;")
            self.emit_grouped(w, rate_body, hoist=True, select=sel)
            w("}")
        nfix = sum(1 for rx in net.reactions
                   if rx["rate"]["kind"] == "expr" and self.rate_is_fixed(rx))
        w("// rates of Tgas and the parameters only (rates_fixed)")
        w("constexpr int NREA_FIXED = %d;" % nfix)
        w("void inline rates(double Tgas, const double* n, const Params& p, "
          "double* k) {")
        w("  for (int i = 0; i < NREA; ++i) k[i] = 0.0;")
        w("  rates_fixed(Tgas, n, p, k);")
        w("  rates_var(Tgas, n, p, k);")
        w("}")
        w("")
        self.emit_freeze(w)

        loop = self.rhs_loop()
        if loop:  # --rhs=loop: fluxes() and rhs() over the index tables
            self.emit_loop_tables(w)
            for fn, body in (("fluxes", "fl"), ("rhs", "dn")):
                w("void inline %s(const double* n, const double* k, "
                  "double* %s) {" % (fn, body))
                if fn == "rhs":
                    w("  for (int i = 0; i < NSPEC; ++i) dn[i] = 0.0;")
                w("  for (int r = 0; r < NREA; ++r) {")
                w("    double flux = k[r];")
                w("    for (int q = lt_kf_off[r]; q < lt_kf_off[r + 1]; ++q)")
                w("      flux *= KROME_LT_DEN(lt_kf_sp[q], lt_kf_pr[q]);")
                if fn == "fluxes":
                    w("    fl[r] = flux;")
                else:
                    w("    for (int q = lt_sc_off[r]; q < lt_sc_off[r + 1]; "
                      "++q)")
                    w("      KROME_LT_ADD(dn[lt_sc_row[q]], lt_sc_c[q], "
                      "flux);")
                w("  }")
                w("}")
                w("")
        # fluxes() and rhs(), one statement per reaction (--rhs=loop: the
        # loops over the index tables above instead)
        if not loop:
            self._emit_unrolled_fluxes_rhs(sec["rates"].append)
        w = sec["rates"].append
        w("#undef Tdust")

        # ==== thermal section ================================================
        w = sec["thermal"].append
        w(self._tdust_define())
        w("")
        # network-defined cooling/heating terms (network-format.md 6.1):
        # volumetric [erg cm^-3 s^-1]; standard thermal processes (atomic,
        # H2, ...) come from runtime modules and are added by the thermal
        # solver
        w("double inline cooling_custom(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        w("  double L = 0.0;")
        for c in net.cooling:
            w("  L += %s;  // %s" % (self.cexpr(c["expr"]), c["name"]))
        w("  return L;")
        w("}")
        w("double inline heating_custom(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        w("  double H = 0.0;")
        for c in net.heating:
            w("  H += %s;  // %s" % (self.cexpr(c["expr"]), c["name"]))
        w("  return H;")
        w("}")
        # --heating=NUCLEAR: the energy released per reaction event (the
        # ; qeff= annotation, MeV) times the reaction flux
        qrx = [rx for rx in net.reactions if rx.get("qeff") is not None]
        if self.thermal["heat_nuclear"] and not qrx:
            self.warnings.append("thermal.heating = ['NUCLEAR']: no reaction has a qeff")
        w("// nuclear energy release [erg cm^-3 s^-1]: sum of flux * qeff "
          "(MeV)")
        w("double inline heating_nuclear(double Tgas, const double* n, "
          "const double* k, const Params& p) {")
        w("  (void)Tgas; (void)n; (void)k; (void)p;")
        w("  double h = 0.0;")
        for rx in qrx:
            fl = "*".join(["k[%d]" % (rx["id"] - 1)]
                          + [self.nden(r["sp"]) for r in rx["reactants"]
                             if r["kinetic"]])
            w("  h += %s * (%s);  // %s" % (fl, self.cexpr(rx["qeff"]),
                                             rx["eq"]))
        w("  return h * constants::eV_to_erg * 1e6;  // MeV -> erg")
        w("}")
        w("// qeff of every reaction [MeV] (0: none)")
        w("inline void reaction_qeff(double Tgas, const double* n, "
          "const Params& p, double* q) {")
        w("  (void)Tgas; (void)n; (void)p;")
        w("  for (int i = 0; i < NREA; ++i) q[i] = 0.0;")
        for rx in qrx:
            w("  q[%d] = %s;" % (rx["id"] - 1, self.cexpr(rx["qeff"])))
        w("}")
        w("constexpr bool HAS_CUSTOM_COOLING = %s;"
          % ("true" if net.cooling else "false"))
        w("constexpr bool HAS_CUSTOM_HEATING = %s;"
          % ("true" if net.heating else "false"))
        w("")
        self.emit_gamma(w)

        # cosmic-ray heating (legacy -heating=CR): CR reactions are those
        # whose rate references the crate param (the KROME-X analog of the legacy
        # @CR blocks). Q_H = Q_He = 4.3 eV (Glassgold+2012), Q_H2(n_H2) from
        # the arXiv:1502.03380 fit, 20 eV for anything else (legacy Hfact).
        cr_rx = [rx for rx in net.reactions
                 if rx["rate"]["kind"] == "expr"
                 and any(nd.get("kind") == "param" and nd.get("s") == "crate"
                         for nd in self.deep_walk(rx["rate"]["ast"]))
                 and len(rx["reactants"]) == 1]
        # Q_H2 [eV] in log n_H2 (legacy's branches): the quartic crc for
        # 1e5 <= n_H2 <= 1e10, its value at n_H2 = 1e10 (log 10) above;
        # the line cr_lin below 1e5, down to minH2 where it reaches 0;
        # 0 below that
        crc = [9.1462198642562, -0.443120145068494, 0.601271617560071,
               -0.0710284101055109, 0.00242079494592405]
        cr_lin = (1.4510849135088, 7.23277032061136)
        minH2 = 10.0 ** (-cr_lin[1] / cr_lin[0])
        evalMax = sum(c * 10.0 ** i for i, c in enumerate(crc))
        w("double inline heating_CR(const double* n, const double* k, "
          "const Params& p) {")
        w("  (void)n; (void)k; (void)p;")
        w("  double h = 0.0;")
        if cr_rx:
            w("  const double ev2erg = constants::eV_to_erg;")
            def h2ion(r):
                a = r["reactants"][0]["sp"]
                return (a in self.H2_FORMS and sorted(
                    p0["sp"] for p0 in r["products"])
                        == sorted([a + "+", "e-"]))
            uses_H2 = any(h2ion(r) for r in cr_rx)
            if uses_H2:
                w("  // Q_H2 [eV] vs n_H2, fit from arXiv:1502.03380")
                w("  double QH2;")
                w("  const double logH2 = "
                  "std::log10(std::fmax(nH2_(n), 1e-40));")
                w("  if (nH2_(n) > 1e10) QH2 = %s;" % self.cnum(evalMax))
                w("  else if (nH2_(n) >= %s && nH2_(n) < 1e5)"
                  % self.cnum(minH2))
                w("    QH2 = %s*logH2 + %s;" % (self.cnum(cr_lin[0]),
                                                self.cnum(cr_lin[1])))
                w("  else if (nH2_(n) < %s) QH2 = 0.0;" % self.cnum(minH2))
                w("  else QH2 = %s;"
                  % " + ".join("%s*std::pow(logH2, %d.0)"
                               % (self.cnum(c), i) if i > 1 else
                               ("%s*logH2" % self.cnum(c) if i == 1
                                else self.cnum(c))
                               for i, c in enumerate(crc)))
                w("  QH2 *= ev2erg;")
            for rx in cr_rx:
                rea = rx["reactants"][0]["sp"]
                prods = sorted(p0["sp"] for p0 in rx["products"])
                if rea == "H" and prods == ["H+", "e-"]:
                    q = "4.3*ev2erg"
                elif rea == "He" and prods == ["He+", "e-"]:
                    q = "4.3*ev2erg  /* = Q_H, Glassgold+2012 */"
                elif h2ion(rx):
                    q = "QH2"
                else:
                    q = "20.0*ev2erg"
                w("  h += k[%d] * n[IDX_%s] * (%s);  // %s"
                  % (rx["id"] - 1, self.cid(rea), q, rx["eq"]))
        w("  return h;")
        w("}")
        w("")

        # chemical heating/cooling (legacy -heating=CHEM / -cooling=CHEM):
        # the H2 chemistry energy ledger — formation channels heat by their
        # exothermicity times the thermalized fraction h2heatfac (HM79
        # critical density), dissociation channels cool by 4.48 eV, and the
        # Solomon process pumps (18.7 h2heatfac + 0.4) eV per dissociation.
        chem_terms = []   # (sign_expr_eV, rx, use_ffac)
        CHEM_RPK = [
            (["H", "H", "H"], ["H", "H2"], "4.48", True),
            (["H", "H", "H2"], ["H2", "H2"], "4.48", True),
            (["H", "H-"], ["H2", "e-"], "3.53", True),
            (["H", "H2+"], ["H+", "H2"], "1.83", True),
            (["H", "H2"], ["H", "H", "H"], "-4.48", False),
            (["H2", "e-"], ["H", "H", "e-"], "-4.48", False),
            (["H2", "H2"], ["H", "H", "H2"], "-4.48", False),
        ]
        solomon = []
        dust_h2 = []  # H + {H} -> H2 with an H2-on-dust rate
        for rx in net.reactions:
            R = sorted(r["sp"] for r in rx["reactants"])
            P = sorted(p0["sp"] for p0 in rx["products"])
            for rr, pp, q, ffac in CHEM_RPK:
                if R == rr and P == pp:
                    chem_terms.append((q, rx, ffac))
                    break
            ast0 = rx["rate"].get("ast")
            if R == ["H", "H"] and P == ["H2"] and ast0 is not None and any(
                    nd.get("f") in H2DUST_EXTERNS
                    for nd in self.deep_walk(ast0)):
                dust_h2.append(rx)
            if self.is_solomon(rx):
                solomon.append(rx)
        # legacy: -heating=CHEM = the formation heating (+ Solomon pumping),
        # -cooling=CHEM / DISS = the collisional dissociation cooling
        heat_terms = [x for x in chem_terms if not x[0].startswith("-")]
        diss_terms = [x for x in chem_terms if x[0].startswith("-")]
        w("double inline cooling_chem(double Tgas, const double* n, "
          "const double* k, const Params& p) {")
        w("  (void)Tgas; (void)n; (void)k; (void)p;")
        w("  double c = 0.0;")
        for q, rx, ffac in diss_terms:
            prod = "*".join("n[IDX_%s]" % self.cid(r["sp"])
                            for r in rx["reactants"])
            w("  c += k[%d] * %s * %s;  // %s"
              % (rx["id"] - 1, q[1:], prod, rx["eq"]))
        if diss_terms:
            w("  c *= constants::eV_to_erg;  // eV -> erg")
        w("  return c;  // collisional H2 dissociation, 4.48 eV each")
        w("}")
        chem_terms = heat_terms
        w("double inline heating_chem(double Tgas, const double* n, "
          "const double* k, const Params& p) {")
        w("  (void)Tgas; (void)n; (void)k; (void)p;")
        w("  double h = 0.0;")
        if chem_terms or solomon or dust_h2:
            w("  const double ev2erg = constants::eV_to_erg;")
            w("  // thermalized fraction 1/(1 + ncr/n) (HM79)")
            w("  const double dd = Hnuclei_(n);")
            w("  const double ncr = 1.0e6 / std::sqrt(Tgas)")
            w("      / (1.6 * std::exp(-std::pow(4.0e2 / Tgas, 2.0)) "
              "* n[IDX_H] / dd")
            w("         + 1.4 * std::exp(-1.2e4 / (Tgas + 1.2e3)) "
              "* nH2_(n) / dd);")
            w("  const double ffac = 1.0 / (1.0 + ncr / dd);")
            for q, rx, ffac in chem_terms:
                prod = "*".join("n[IDX_%s]" % self.cid(r["sp"])
                                for r in rx["reactants"])
                qexpr = "%s*ffac" % q if ffac else q
                w("  h += k[%d] * (%s) * %s;  // %s"
                  % (rx["id"] - 1, qexpr, prod, rx["eq"]))
            for rx in solomon:
                w("  h += k[%d] * (18.7*ffac + 0.4) * n[IDX_%s];"
                  "  // Solomon: %s" % (rx["id"] - 1,
                                       self.cid(rx["reactants"][0]["sp"]),
                                       rx["eq"]))
            for rx in dust_h2:
                prod = "*".join(["k[%d]" % (rx["id"] - 1)]
                                + ["n[IDX_%s]" % self.cid(r["sp"])
                                   for r in rx["reactants"] if r["kinetic"]])
                w("  h += %s * (4.2*ffac + 0.2);  // H2 on dust: %s"
                  % (prod, rx["eq"]))
            w("  h *= ev2erg;")
        w("  return h;  // H2 formation heating + Solomon pumping")
        w("}")
        w("")
        self.emit_h2_extra(w)

        # bremsstrahlung / free-free cooling (legacy -cooling=BSS): Z^2-
        # weighted sum over ATOMIC ions only (legacy brem_ions pragma;
        # molecular ions excluded), Gaunt-averaged 1.42e-27 * gff prefactor
        # with gff = 1.5 (legacy choice).
        ff_ions = [(s, net.species[s]["charge"]) for s in sp_sorted
                   if net.species[s]["charge"] > 0
                   and net.species[s]["composition"]
                   and sum(net.species[s]["composition"].values()) == 1]
        w("double inline cooling_ff(double Tgas, const double* n) {")
        w("  (void)Tgas; (void)n;")
        if ff_ions and "e-" in net.species:
            w("  const double bms_ions =")
            w("      %s;" % "\n      + ".join(
                ("%d.0*n[IDX_%s]" % (z * z, self.cid(s))) if z > 1
                else "n[IDX_%s]" % self.cid(s) for s, z in ff_ions))
            w("  return std::fmax(1.42e-27 * 1.5 * std::sqrt(Tgas)"
              " * bms_ions * n[IDX_E], 0.0);")
        else:
            w("  return 0.0;  // no atomic ions / no electrons in network")
        w("}")
        w("")

        self.emit_mol_cooling(w)
        self.emit_metal_cooling(w)
        self.emit_dust_tables(w)
        self.emit_dust_heating(w)
        self.emit_dh_thermal(w)
        self.emit_recomb_cooling(w)
        w("#undef Tdust")

        # ==== jacobian section ===============================================
        w = sec["jacobian"].append
        w(self._tdust_define())
        w("")

        # rhs_scale(): per-species sum of |flux| contributions; an upper
        # bound on cancellation in dn (used for FD validation noise floors
        # and useful for solver error weighting)
        w("void inline rhs_scale(const double* n, const double* k, "
          "double* s) {")
        w("  for (int i = 0; i < NSPEC; ++i) s[i] = 0.0;")
        w("  double flux;")
        w("  (void)n; (void)k; (void)flux;  // an empty network")
        for rx in net.reactions:
            i = rx["id"] - 1
            fac = ["k[%d]" % i] + [self.nden(r["sp"])
                                   for r in rx["reactants"] if r["kinetic"]]
            w("  flux = std::fabs(%s);" % "*".join(fac))
            rows = [r["sp"] for r in rx["reactants"]] \
                + [pr["sp"] for pr in rx["products"]]
            if self.gas_total:  # totals also see the exchange noise floor
                rows += [self.total_of[sp] for sp in rows
                         if sp in self.total_of]
            for sp in rows:
                w("  s[IDX_%s] += flux;" % self.cid(sp))
        w("}")
        w("")

        # numerical dk/dn helpers
        w("// difference quotient of one rate along species j: central where")
        w("// positivity allows (strongly curved rates, e.g. dissH2_Martin96),")
        w("// else one-sided; a difference at the roundoff of k is noise, not a")
        w("// derivative (~4 ulps: 0)")
        w("template <class F>")
        w("static inline double fd_dk_(F rate, double* np_, const double* n, "
          "double k0, int j) {")
        w("  const double h_ = std::fmax(1e-7 * std::fabs(n[j]), 1e-30);")
        w("  np_[j] = n[j] + h_;")
        w("  const double kp_ = rate(np_);")
        w("  double km_ = k0, den_ = h_;")
        w("  if (n[j] > h_) {")
        w("    np_[j] = n[j] - h_;")
        w("    km_ = rate(np_);")
        w("    den_ = 2.0 * h_;")
        w("  }")
        w("  np_[j] = n[j];")
        w("  if (std::fabs(kp_ - km_) <= 1e-15 * std::fabs(k0)) return 0.0;")
        w("  return (kp_ - km_) / den_;")
        w("}")
        w("// dk/d(aggregate) c[a] for the aggregates in mask (bit a): the")
        w("// difference quotients g_q along the dominant carrier of each")
        w("// aggregate (largest |w_aj n_j| outside the rate's direct species,")
        w("// which moves the aggregate well above roundoff) satisfy")
        w("// g_q = sum_a c_a w_a,j_q; solved by elimination (m <= 3)")
        w("template <class F>")
        w("static inline void fd_dk_agg_(F rate, double* np_, const double* n,"
          " double k0, int mask, const int* dir, int ndir, double* c) {")
        w("  int am[NAGG], pj[NAGG], m = 0;")
        w("  for (int a = 0; a < NAGG; ++a) {")
        w("    c[a] = 0.0;")
        w("    if (mask & (1 << a)) am[m++] = a;")
        w("  }")
        w("  double M[NAGG][NAGG + 1];")
        w("  for (int q = 0; q < m; ++q) {")
        w("    int best = -1;")
        w("    double bv = 0.0;")
        w("    for (int j = 0; j < NSPEC; ++j) {")
        w("      bool skip = false;")
        w("      for (int t = 0; t < ndir; ++t) skip = skip || dir[t] == j;")
        w("      for (int t = 0; t < q; ++t) skip = skip || pj[t] == j;")
        w("      const double v = std::fabs(agg_w_(am[q], j) * n[j]);")
        w("      if (!skip && v > bv) bv = v, best = j;")
        w("    }")
        w("    if (best < 0) return;  // no carrier: no dependence to resolve")
        w("    pj[q] = best;")
        w("    for (int r = 0; r < m; ++r) M[q][r] = agg_w_(am[r], best);")
        w("    M[q][m] = fd_dk_(rate, np_, n, k0, best);")
        w("  }")
        w("  for (int col = 0; col < m; ++col) {  // partial pivoting")
        w("    int piv = col;")
        w("    for (int r = col + 1; r < m; ++r)")
        w("      if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;")
        w("    if (M[piv][col] == 0.0) return;")
        w("    for (int t = 0; t <= m; ++t) {  // (std::swap is host-only)")
        w("      const double x = M[piv][t];")
        w("      M[piv][t] = M[col][t];")
        w("      M[col][t] = x;")
        w("    }")
        w("    for (int r = col + 1; r < m; ++r) {")
        w("      const double f = M[r][col] / M[col][col];")
        w("      for (int t = col; t <= m; ++t) M[r][t] -= f * M[col][t];")
        w("    }")
        w("  }")
        w("  double x[NAGG];")
        w("  for (int r = m - 1; r >= 0; --r) {")
        w("    double v = M[r][m];")
        w("    for (int t = r + 1; t < m; ++t) v -= M[r][t] * x[t];")
        w("    x[r] = v / M[r][r];")
        w("  }")
        w("  for (int r = 0; r < m; ++r) c[am[r]] = x[r];")
        w("}")
        w("// the aggregate part sum_a c_a w_aj of a total derivative along j")
        w("static inline double agg_part_(const double* c, int j) {")
        w("  double s = 0.0;")
        w("  for (int a = 0; a < NAGG; ++a) s += c[a] * agg_w_(a, j);")
        w("  return s;")
        w("}")
        w("")
        # jacobian_split(): the direct part J[i*NSPEC+j] and the aggregate
        # columns U[i*NAGG+a]; jacobian() = J + U W^T
        w("// analytic Jacobian split: d dn_i / d n_j = J[i*NSPEC+j] + sum_a")
        w("// U[i*NAGG+a] agg_w_(a, j): J the direct dependence (the network's")
        w("// sparsity), U the coupling through the aggregates ntot, Hnuclei,")
        w("// rho (dense in j but of rank <= NAGG)")
        w("// J element (i, j) at J[i*jrs_ + j*jcs_]: dense row-major by default,")
        w("// or a cell's block of an interleaved matrix in its own layout")
        w("void inline jacobian_split(double Tgas, const double* n, "
          "const Params& p, const double* k, double* J, double* U, "
          "long jrs_ = NSPEC, long jcs_ = 1) {")
        w("  (void)Tgas; (void)n; (void)p; (void)U;")
        w("  for (int i = 0; i < NSPEC; ++i)")
        w("    for (int j = 0; j < NSPEC; ++j) J[i*jrs_ + j*jcs_] = 0.0;")
        w("  for (int i = 0; i < NSPEC*NAGG; ++i) U[i] = 0.0;")
        w("  double d;")
        w("  (void)k; (void)d;  // an empty network")
        w("  // stoichiometric part: k[i] * d(prod n)/dn_j; k already")
        w("  // carries the piecewise selection and clamping")
        if self.rhs_loop():  # --rhs=loop: over the index tables
            w("  for (int r = 0; r < NREA; ++r)")
            w("    for (int t = lt_jt_off[r]; t < lt_jt_off[r + 1]; ++t) {")
            w("      d = lt_jt_m[t] > 1 ? static_cast<double>(lt_jt_m[t]) * k[r]")
            w("                         : k[r];")
            w("      for (int q = lt_jt_foff[t]; q < lt_jt_foff[t + 1]; ++q)")
            w("        d *= KROME_LT_DEN(lt_jf_sp[q], lt_jf_pr[q]);")
            w("      for (int c = lt_jt_coff[t]; c < lt_jt_coff[t + 1]; ++c)")
            w("        for (int q = lt_sc_off[r]; q < lt_sc_off[r + 1]; ++q)")
            w("          KROME_LT_ADD(J[lt_sc_row[q]*jrs_ + lt_jc_col[c]*jcs_],")
            w("                       lt_sc_c[q] * lt_jc_sgn[c], d);")
            w("    }")
        else:
            for rx in net.reactions:
                for line in self.jac_stoich_lines(rx):
                    w("  " + line)
        w("  // rate-derivative part: (dk/dn_j) * prod n, gated and clamped")
        w("  // exactly like the k evaluation itself")
        self.emit_grouped(w, self.jac_dk_lines)
        w("}")
        w("// the full analytic Jacobian J[i*NSPEC+j] = d dn_i / d n_j")
        w("void inline jacobian(double Tgas, const double* n, "
          "const Params& p, const double* k, double* J) {")
        w("  std::vector<double> U(static_cast<size_t>(NSPEC) * NAGG);")
        w("  jacobian_split(Tgas, n, p, k, J, U.data());")
        w("  for (int i = 0; i < NSPEC; ++i)")
        w("    for (int a = 0; a < NAGG; ++a) {")
        w("      const double u = U[static_cast<size_t>(i) * NAGG + a];")
        w("      if (u == 0.0) continue;")
        w("      for (int j = 0; j < NSPEC; ++j)")
        w("        J[static_cast<size_t>(i) * NSPEC + j] += u * agg_w_(a, j);")
        w("    }")
        w("}")
        self.emit_sparse_tables(w)
        w("#undef Tdust")
        core = sec["core"]
        k = core.index(self.H2OP_SLOT)
        core[k:k + 1] = (self.h2op_lines()
                         if getattr(self, "_need_h2op", False) else [])
        self.emit_params(sec)
        return sec

    def load_ebind(self):
        """data/Ebare_ice.dat -> {species: (Ebare [K], Eice [K])}."""
        out = {}
        f = self.data_dir.parent / "Ebare_ice.dat"
        if f.exists():
            for line in f.read_text().splitlines():
                s = line.split("#")[0].split()
                if len(s) == 3:
                    out[s[0]] = (float(s[1]), float(s[2]))
        return out

    def load_table(self, fname):
        """Tabulated file -> (kind, data). 1D: two columns (x, v).
        2D (legacy @tabvar format): 'nx, ny' size header then row-major
        x y v rows on a rectilinear grid."""
        path = Path(fname)
        if not path.exists():
            path = self.data_dir.parent.parent / fname
        if not path.exists():
            raise FileNotFoundError("table file %s not found" % fname)
        rows = []
        size = None
        for line in path.read_text().splitlines():
            s = line.split("#")[0].replace(",", " ").split()
            if not s:
                continue
            v = [float(x) for x in s]
            if size is None and len(v) == 2 and v[0] == int(v[0]) \
                    and v[1] == int(v[1]) and not rows:
                size = (int(v[0]), int(v[1]))
                continue
            rows.append(v)
        if size:
            nx, ny = size
            xs = sorted({r[0] for r in rows})
            ys = sorted({r[1] for r in rows})
            if len(xs) != nx or len(ys) != ny or len(rows) != nx * ny:
                raise ValueError("table %s is not a %dx%d grid"
                                 % (fname, nx, ny))
            vals = {(r[0], r[1]): r[2] for r in rows}
            grid = [[vals[(x, y)] for y in ys] for x in xs]
            return ("2d", (xs, ys, grid))
        return ("1d", ([r[0] for r in rows], [r[1] for r in rows]))

    def load_mol_cool(self, coolant):
        """Legacy molecular cooling table (data/cool<X>.dat): lines of
        'i j k  log10(T) log10(nH+nH2) log10(Ncol) log10(cool)' on a
        uniform grid in each log variable (1-based Fortran indices)."""
        fname = self.MOL_COOLANTS[coolant]
        for base in (self.data_dir.parent, self.data_dir,
                     self.data_dir.parent.parent):
            path = base / fname
            if path.exists():
                break
        else:
            raise FileNotFoundError("cooling table %s not found" % fname)
        x1, x2, x3, y = {}, {}, {}, {}
        for line in path.read_text().splitlines():
            s = line.split("#")[0].split()
            if len(s) < 7:
                continue
            i, j, k = int(s[0]), int(s[1]), int(s[2])
            x1[i], x2[j], x3[k] = float(s[3]), float(s[4]), float(s[5])
            y[(i, j, k)] = float(s[6])
        n1, n2, n3 = max(x1), max(x2), max(x3)
        if len(y) != n1 * n2 * n3:
            raise ValueError("table %s is not a full %dx%dx%d grid"
                             % (fname, n1, n2, n3))
        ax1 = [x1[i] for i in range(1, n1 + 1)]
        ax2 = [x2[j] for j in range(1, n2 + 1)]
        ax3 = [x3[k] for k in range(1, n3 + 1)]
        # flattened C order: y[((i*n2)+j)*n3+k], 0-based
        flat = [y[(i + 1, j + 1, k + 1)] for i in range(n1)
                for j in range(n2) for k in range(n3)]
        return (ax1, ax2, ax3, flat)

    def emit_mol_cooling(self, w):
        """cooling_mol(): legacy -cooling=CO/OH/H2O/HCN, Neufeld & Kaufman
        1993 eqn. 5 — trilinear interpolation of log10 W(T, nH+nH2, Ncol),
        cooling = 10^logW * (nH+nH2) * n_coolant [erg cm^-3 s^-1]. Ported
        verbatim from legacy cooling_CO & init_coolingCO (krome_cooling.f90);
        the num2col column-density estimator is a codegen option exactly
        like the legacy -columnDensityMethod flag."""
        active = []
        for c in self.mol_cooling:
            if c not in self.MOL_COOLANTS:
                self.warnings.append("unknown molecular coolant %s" % c)
            elif c not in self.net.species:
                self.warnings.append("molecular coolant %s not in network, "
                                     "skipped" % c)
            elif not ("H" in self.net.species and self.h2_forms()):
                self.warnings.append("molecular coolant %s needs H and H2 "
                                     "colliders, skipped" % c)
            else:
                active.append(c)
        if active:
            self._emit_interp3_body(w)

        for c in active:
            data = self.load_mol_cool(c)
            self._emit_cool3d_table(
                w, "molcool_%s" % self.cid(c), data,
                "%s cooling table %s (%dx%dx%d), K. Omukai / NK93"
                % ((c, self.MOL_COOLANTS[c]) + tuple(len(a)
                                                     for a in data[:3])))
        if active:
            w("#define KROME_HAS_COOL_MOL 1")
        w("double inline cooling_mol(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if active:
            w("  double L = 0.0;")
            w("  const double cH = n[IDX_H] + nH2_(n);")
            w("  if (cH <= 0.0) return 0.0;")
            w("  const double lT = std::log10(Tgas);")
            w("  const double lcH = std::log10(cH);")
            for c in active:
                cc = self.cid(c)
                w("  L += std::pow(10.0, mol_cool_interp_(molcool_%s_tab, "
                  "lT, lcH," % cc)
                w("       std::log10(num2col_(n[IDX_%s], Tgas, n, p)))) "
                  "* cH * n[IDX_%s];" % (cc, cc))
            w("  return L;")
        else:
            w("  return 0.0;  // no molecular coolants configured")
        w("}")
        # the coolants selected for the cooling floor (--cooling-floor)
        fl = [c for c in active if c in self.mol_floor]
        w("// molecular coolants under the cooling floor: %s"
          % (", ".join(fl) or "none"))
        w("double inline cooling_mol_floor(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if fl:
            w("  double L = 0.0;")
            w("  const double cH = n[IDX_H] + nH2_(n);")
            w("  if (cH <= 0.0) return 0.0;")
            w("  const double lT = std::log10(Tgas);")
            w("  const double lcH = std::log10(cH);")
            for c in fl:
                cc = self.cid(c)
                w("  L += std::pow(10.0, mol_cool_interp_(molcool_%s_tab, "
                  "lT, lcH," % cc)
                w("       std::log10(num2col_(n[IDX_%s], Tgas, n, p)))) "
                  "* cH * n[IDX_%s];" % (cc, cc))
            w("  return L;")
        else:
            w("  return 0.0;")
        w("}")
        w("")

    def _emit_interp3_body(self, w):
        """the shared uniform-grid trilinear interpolator (legacy
        cooling_CO / cooling_Z_CIE structure), emitted once."""
        if getattr(self, "_interp3_done", False):
            return
        self._interp3_done = True
        if True:
            w("struct MolCoolTab {")
            w("  int n2, n3;  // y is flattened [i1][i2][i3]")
            w("  const double *x1, *x2, *x3, *ixd1, *ixd2, *ixd3, *y;")
            w("  double x1min, x1max, x2min, x2max, x3min, x3max;")
            w("  double dvn1, dvn2, dvn3;")
            w("};")
            w("// returns log10(W); -1e99 (i.e. 0 cooling) below the table")
            w("inline double mol_cool_interp_(const MolCoolTab& t, "
              "double v1, double v2, double v3) {")
            w("  const double eps = 1e-5;")
            w("  if (v1 >= t.x1max) v1 = t.x1max * (1.0 - eps);")
            w("  if (v2 >= t.x2max) v2 = t.x2max * (1.0 - eps);")
            w("  if (v3 >= t.x3max) v3 = t.x3max * (1.0 - eps);")
            w("  if (v1 < t.x1min || v2 < t.x2min || v3 < t.x3min) "
              "return -1e99;")
            w("  const int i = (int)((v1 - t.x1min) * t.dvn1);")
            w("  const int j = (int)((v2 - t.x2min) * t.dvn2);")
            w("  const int k = (int)((v3 - t.x3min) * t.dvn3);")
            w("  const double prev1 = (v1 - t.x1[i]) * t.ixd1[i];")
            w("  const double prev2 = (v2 - t.x2[j]) * t.ixd2[j];")
            w("  const double* y = t.y;")
            w("  const int n3 = t.n3, n2 = t.n2;")
            w("  #define YT_(kk, jj, ii) y[((ii)*n2 + (jj))*n3 + (kk)]")
            w("  const double vv1 = prev1 * (YT_(k, j, i+1) - YT_(k, j, i))"
              " + YT_(k, j, i);")
            w("  const double vv2 = prev1 * (YT_(k, j+1, i+1) - "
              "YT_(k, j+1, i)) + YT_(k, j+1, i);")
            w("  const double vv12 = prev2 * (vv2 - vv1) + vv1;")
            w("  const double vv3 = prev1 * (YT_(k+1, j, i+1) - "
              "YT_(k+1, j, i)) + YT_(k+1, j, i);")
            w("  const double vv4 = prev1 * (YT_(k+1, j+1, i+1) - "
              "YT_(k+1, j+1, i)) + YT_(k+1, j+1, i);")
            w("  const double vv34 = prev2 * (vv4 - vv3) + vv3;")
            w("  #undef YT_")
            w("  return (v3 - t.x3[k]) * t.ixd3[k] * (vv34 - vv12) + vv12;")
            w("}")

    def _emit_cool3d_table(self, w, name, data, comment):
        """embed one uniform-grid 3D table as <name>_* arrays plus the
        <name>_tab MolCoolTab descriptor."""
        ax1, ax2, ax3, flat = data
        n2, n3 = len(ax2), len(ax3)
        w("// " + comment)
        for nm, ax in (("x1", ax1), ("x2", ax2), ("x3", ax3)):
            w("inline const double %s_%s[%d] = {%s};"
              % (name, nm, len(ax), ", ".join(self.cnum(v) for v in ax)))
            ixd = [1.0 / (ax[i + 1] - ax[i]) for i in range(len(ax) - 1)]
            w("inline const double %s_i%s[%d] = {%s};"
              % (name, nm, len(ixd), ", ".join(self.cnum(v) for v in ixd)))
        self.emit_data_array(w, "%s_y" % name, flat, comment)
        # the descriptor as a device-safe accessor (it points to the
        # accessor-backed arrays, so it is built where it is used)
        w("KROME_HD inline MolCoolTab %s_tab_get_() {" % name)
        w("  return {%d, %d," % (n2, n3))
        w("    %s_x1, %s_x2, %s_x3," % (name, name, name))
        w("    %s_ix1, %s_ix2, %s_ix3, %s_y," % (name, name, name, name))
        w("    %s, %s, %s, %s, %s, %s,"
          % tuple(self.cnum(v) for v in (min(ax1), max(ax1), min(ax2),
                                         max(ax2), min(ax3), max(ax3))))
        w("    %s, %s, %s};"
          % tuple(self.cnum((len(a) - 1) / (max(a) - min(a)))
                  for a in (ax1, ax2, ax3)))
        w("}")
        w("#define %s_tab (::krome::%s_tab_get_())" % (name, name))

    # -- metal line cooling (legacy -cooling=Z) --------------------------------

    COOLZ_FUNCS = ("exp", "log", "log10", "sqrt", "min", "max", "tanh",
                   "abs", "where")

    @staticmethod
    def _godunov_spline(x, y):
        """legacy fspline coefficients (src/krome_fit.f90 spline(): A.
        Godunov's cubic spline with its end conditions): b, c, d with
        s = y_i + b_i dx + c_i dx^2 + d_i dx^3 on [x_i, x_i+1]"""
        n = len(x)
        b, c, d = [0.0] * n, [0.0] * n, [0.0] * n
        if n < 3:
            b[0] = b[1] = (y[1] - y[0]) / (x[1] - x[0])
            return b, c, d
        gap = n - 1
        d[0] = x[1] - x[0]
        c[1] = (y[1] - y[0]) / d[0]
        for i in range(1, gap):
            d[i] = x[i + 1] - x[i]
            b[i] = 2.0 * (d[i - 1] + d[i])
            c[i + 1] = (y[i + 1] - y[i]) / d[i]
            c[i] = c[i + 1] - c[i]
        b[0], b[n - 1] = -d[0], -d[n - 2]
        c[0] = c[n - 1] = 0.0
        if n != 3:
            c[0] = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0])
            c[n - 1] = (c[n - 2] / (x[n - 1] - x[n - 3])
                        - c[n - 3] / (x[n - 2] - x[n - 4]))
            c[0] = c[0] * d[0] ** 2 / (x[3] - x[0])
            c[n - 1] = -c[n - 1] * d[n - 2] ** 2 / (x[n - 1] - x[n - 4])
        for i in range(1, n):
            h = d[i - 1] / b[i - 1]
            b[i] = b[i] - h * d[i - 1]
            c[i] = c[i] - h * c[i - 1]
        c[n - 1] = c[n - 1] / b[n - 1]
        for j in range(1, gap + 1):
            i = n - 1 - j
            c[i] = (c[i] - d[i] * c[i + 1]) / b[i]
        b[n - 1] = ((y[n - 1] - y[gap - 1]) / d[gap - 1]
                    + d[gap - 1] * (c[gap - 1] + 2.0 * c[n - 1]))
        for i in range(gap):
            b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i])
            d[i] = (c[i + 1] - c[i]) / d[i]
            c[i] = 3.0 * c[i]
        c[n - 1] = 3.0 * c[n - 1]
        d[n - 1] = d[n - 2]
        return b, c, d

    def _coolz_tables(self, text):
        """fspline((/x.../), (/y.../), ARG) and flin(...) of a coolant
        rate -> cftab_<k>(ARG): generated tables (self._cftabs), the
        spline coefficients computed here once (legacy rebuilt the spline
        at every call)"""
        out, pos = [], 0
        rx = re.compile(r"\b(fspline|flin)\s*\(", re.I)
        tabs = self.__dict__.setdefault("_cftabs", [])
        while True:
            m = rx.search(text, pos)
            if not m:
                out.append(text[pos:])
                return "".join(out)
            out.append(text[pos:m.start()])
            depth, k = 1, m.end()
            args, start = [], m.end()
            while depth:
                ch = text[k]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        args.append(text[start:k])
                elif ch == "," and depth == 1:
                    args.append(text[start:k])
                    start = k + 1
                k += 1
            if len(args) != 3:
                raise ValueError("%s: needs 3 arguments in %r"
                                 % (m.group(1), text))

            def arr(a):
                a = a.strip()
                if not (a.startswith("(/") and a.endswith("/)")):
                    raise ValueError("%s: a (/.../) array expected, got %r"
                                     % (m.group(1), a))
                return [float(d2e(v)) for v in a[2:-2].split(",")
                        if v.strip()]
            X, Y = arr(args[0]), arr(args[1])
            if len(X) != len(Y) or len(X) < 2 or any(
                    X[i + 1] <= X[i] for i in range(len(X) - 1)):
                raise ValueError("%s: bad table in %r" % (m.group(1), text))
            kind = m.group(1).lower()
            if kind == "fspline":
                B, C, D = self._godunov_spline(X, Y)
                vals = X + Y + B + C + D
            else:
                vals = X + Y
            name = "cftab_%d" % len(tabs)
            tabs.append((name, kind, len(X), vals, False))
            out.append("%s(%s)" % (name, args[2].strip()))
            pos = k

    def _emit_cftabs(self, w):
        """the coolant tables registered since the last call"""
        tabs = self.__dict__.setdefault("_cftabs", [])
        if tabs and not getattr(self, "_cftab_helpers", False):
            self._cftab_helpers = True
            w("// coolant-file rate tables (legacy fspline: Godunov cubic")
            w("// spline, t = X, Y, B, C, D; flin: linear, t = X, Y), the")
            w("// end values outside the table")
            w("KROME_HD inline double cftab_spline_(const double* t, int n,"
              " double u) {")
            w("  const double* X = t;")
            w("  if (u <= X[0]) return t[n];")
            w("  if (u >= X[n - 1]) return t[2 * n - 1];")
            w("  int i = 0, j = n;")
            w("  while (j > i + 1) {")
            w("    const int k = (i + j) / 2;")
            w("    if (u < X[k]) j = k; else i = k;")
            w("  }")
            w("  const double dx = u - X[i];")
            w("  return t[n + i] + dx * (t[2 * n + i] + dx * (t[3 * n + i]"
              " + dx * t[4 * n + i]));")
            w("}")
            w("KROME_HD inline double cftab_lin_(const double* t, int n, "
              "double u) {")
            w("  const double* X = t;")
            w("  const double x = std::fmin(std::fmax(u, X[0]), X[n - 1]);")
            w("  for (int i = 1; i < n; ++i)")
            w("    if (x <= X[i])")
            w("      return (t[n + i] - t[n + i - 1]) / (X[i] - X[i - 1])"
              " * (x - X[i - 1]) + t[n + i - 1];")
            w("  return t[2 * n - 1];")
            w("}")
        for q, (name, kind, n, vals, done) in enumerate(tabs):
            if done:
                continue
            self.emit_data_array(w, name + "_t", vals,
                                 "coolant %s table" % kind)
            w("KROME_HD inline double %s(double u) { return cftab_%s_(%s_t, "
              "%d, u); }" % (name, "spline" if kind == "fspline" else "lin",
                             name, n))
            tabs[q] = (name, kind, n, vals, True)

    def _coolz_expr(self, text, localnames):
        """Fortran-ish coolZ.dat expression -> classified knet AST."""
        s = d2e(self._coolz_tables(text.strip()))
        s = re.sub(r"\b(EXP|LOG10|LOG|SQRT|MIN|MAX|TANH|ABS)\s*\(",
                   lambda m: m.group(1).lower() + "(", s, flags=re.I)
        ast = knet.parse_expr(s, 0)
        for nd in knet.walk(ast):
            if "s" in nd and "kind" not in nd:
                sym = nd["s"]
                if sym == "Tgas" or sym in localnames:
                    nd["kind"] = "state"  # locals print as bare names
                elif sym in knet.CONSTANTS:
                    nd["kind"] = "const"
                else:
                    raise ValueError("unknown symbol %r in coolZ expression "
                                     "%r" % (sym, text))
            elif "f" in nd and "kind" not in nd:
                if nd["f"].startswith("cftab_"):
                    nd["kind"] = "tfunc"  # a generated table
                elif nd["f"] not in self.COOLZ_FUNCS:
                    raise ValueError("unknown function %r in coolZ "
                                     "expression %r" % (nd["f"], text))
        return ast

    def load_coolz(self):
        """Legacy data/coolZ.dat, the embedded metal coolants (see
        parse_coolfile)."""
        for base in (self.data_dir.parent, self.data_dir,
                     self.data_dir.parent.parent):
            path = base / "coolZ.dat"
            if path.exists():
                break
        else:
            raise FileNotFoundError("coolZ.dat not found")
        return self.parse_coolfile(path)

    def load_cool_files(self, paths):
        """the coolants of the user coolant files (thermal.coolant_files), in file
        order: each with its species (None if not in the network), its
        legacy cooling flag, its file's @TgasUpperLimit and @var list"""
        out, seen = [], {}
        for path in paths:
            path = resolve_data_path(path, self.data_dir)
            limit, gvars, metals = self.parse_coolfile(path)
            for name, md in metals.items():
                sp = self.match_species_(name)
                md["name"], md["species"] = name, sp
                md["gvars"], md["file"] = gvars, str(path)
                md["flag"] = self.cool_flag_(sp or name)
                f = md["flag"].upper()
                if f in self.COOLING_NAMES:
                    raise ValueError(
                        "%s: coolant %s is named like the built-in cooling "
                        "%s (rename the metal: line)" % (path, name, f))
                if f in seen:
                    raise ValueError("coolant %s defined twice (%s, %s)"
                                     % (md["flag"], seen[f], path))
                seen[f] = str(path)
                out.append(md)
        return out

    def cool_flag_(self, sp):
        """the legacy cooling name of a coolant: atoms and atomic ions in
        Roman numerals (C -> CI, Fe+ -> FeII, C- -> CmII), molecules as
        the species (O2)"""
        base = sp.rstrip("+-")
        info = self.net.species.get(sp)
        if info is None:  # not in the network: parse the name (CHIANTI
            # upper case: FE+ -> Fe+)
            for cand in (sp, base.capitalize() + sp[len(base):]):
                try:
                    info = knet.parse_species(cand)
                except Exception:
                    continue
                if info.get("composition"):
                    sp, base = cand, cand.rstrip("+-")
                    break
        comp = (info or {}).get("composition") or {}
        atom = len(comp) == 1 and list(comp.values())[0] == 1
        if not atom:
            return sp
        q = sp.count("+") - sp.count("-")

        def roman(k):
            r = ""
            for v, t in ((10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
                         (1, "I")):
                while k >= v:
                    r += t
                    k -= v
            return r
        return base + (roman(q + 1) if q >= 0 else "m" + roman(-q + 1))

    def match_species_(self, name):
        """the network species of a coolant-file name: exact, else
        case-insensitive (CHIANTI writes NE+++++, FE+); 'e' / 'E' the
        electrons; None when absent"""
        if name in ("e", "E"):
            name = "e-"
        if name in self.net.species:
            return name
        hits = [s for s in self.net.species if s.upper() == name.upper()]
        return hits[0] if len(hits) == 1 else None

    def parse_coolfile(self, path):
        """A legacy coolant file (data/coolZ.dat, or a user -coolFile:
        LAMDA / CHIANTI conversions, tests/Cinterface/coolX.dat): per
        coolant the levels (energy [K], g), the radiative transitions
        (up -> low, Aij[, wavelength in Angstrom]) and the per-collider
        de-excitation rate expressions with their 'if(COND): EXPR'
        overrides; '@var:' shortcuts and '@TgasUpperLimit:'. Legacy
        grammar, with these differences: an override stays attached to
        its own rate (legacy attached it to the last kept one, wrong
        after -coolLevels), informational 'key: value' lines are skipped
        (para/ortho, ortho/para, escape), 'E' is the electron collider,
        and thermal.levels = L drops the levels above L with their lines."""
        limit, gvars, metals, cur = 1e99, [], {}, None
        L = self.cool_levels
        for raw in Path(path).read_text().splitlines():
            s = raw.strip()
            if not s or s.startswith("#") or s.startswith("//"):
                continue
            low = s.lower()
            if low.startswith("@tgasupperlimit:"):
                limit = float(d2e(s.split(":", 1)[1].strip()))
            elif low.startswith("@var:"):
                name, expr = s[5:].split("=", 1)
                gvars.append((name.strip(), expr.strip()))
            elif low.startswith("metal:"):
                cur = {"name": s.split(":", 1)[1].strip(), "levels": {},
                       "rad": [], "coll": [], "file": str(path)}
            elif low.startswith("endmetal") or low.startswith("end metal"):
                metals[cur["name"]] = cur
                cur = None
            elif cur is None:
                continue  # outside a coolant block (legacy: ignored)
            elif low.startswith("level"):
                head, rest = s.split(":", 1)
                idx = int(head.split()[1])
                e, g = (float(d2e(x)) for x in rest.split(",")[:2])
                if L is None or idx <= L:
                    cur["levels"][idx] = (e, g)
            elif low.startswith("if("):
                cond, expr = s.split(":", 1)
                if cur["coll"] and cur["coll"][-1] is not None:
                    cur["coll"][-1][3].append((cond[2:].strip(),
                                               expr.strip()))
            elif "->" in s:
                parts = d2e(s.replace("->", ",")).split(",")
                u, d = int(parts[0]), int(parts[1])
                if L is None or (u <= L and d <= L):
                    cur["rad"].append((u, d, float(parts[2]),
                                       float(parts[3]) if len(parts) > 3
                                       else None))
            elif re.match(r"^[A-Za-z][\w/ ]*:", s) and s.count(",") < 3:
                continue  # informational (para/ortho, ortho/para, escape)
            else:
                parts = s.split(",")
                if len(parts) < 4:
                    raise ValueError("%s: bad coolant line %r" % (path, s))
                coll = parts[0].strip()
                if coll == "E":
                    coll = "e"
                u, d = int(parts[1]), int(parts[2])
                if L is not None and (u > L or d > L):
                    cur["coll"].append(None)  # dropped, with its overrides
                    continue
                cur["coll"].append([coll, u, d,
                                    [(None, ",".join(parts[3:]).strip())]])
        for md in metals.values():
            md["coll"] = [c for c in md["coll"] if c is not None]
            md["limit"] = limit
        return limit, gvars, metals

    def emit_metal_cooling(self, w):
        """cooling_Z(): legacy -cooling=Z (Maio et al. 2007 approach) —
        per-metal level populations from a steady-state linear system
        (collisional up/down + radiative decay, ground-level row replaced
        by number conservation), cooling = kB sum B_up Aul dE_K. Ported
        from the legacy coolZ.dat-driven generator; deviations: rates are
        evaluated directly instead of interpolated from the init-time
        log-T table, and negative populations are clamped instead of
        aborting (a solver library cannot stop())."""
        active, limit, gvars, metals = [], 1e99, [], {}
        if self.metal_cooling:
            limit, gvars, metals = self.load_coolz()
            want = list(metals) if self.metal_cooling == "auto" \
                else self.metal_cooling
            fsp = {md["species"] for md in self.file_active}
            for m in want:
                if m in fsp:
                    continue  # from its coolant file
                if m not in metals:
                    self.warnings.append("no coolZ.dat data for metal %s"
                                         % m)
                elif m not in self.net.species:
                    if self.metal_cooling != "auto":
                        self.warnings.append("metal coolant %s not in "
                                             "network, skipped" % m)
                else:
                    active.append(m)
        files = getattr(self, "file_active", [])
        fsp = {md["species"] for md in files}
        # a coolant a user file defines is taken from the file
        active = [m for m in active if m not in fsp]
        if active or files:
            self._emit_coolz_solver(w)
        # the coolant-file coolants (--cool-file) selected by --cooling,
        # FILE or --metal-cooling: their own cooling terms
        localnames = {"Tgas", "invT", "logTgas"} | {nm for nm, _ in gvars}
        for m in active:
            self._emit_coolant(w, metals[m], m, "coolingZ_%s" % self.cid(m),
                               gvars, localnames)
        for md in files:
            fl = {"Tgas", "invT", "logTgas"} | {nm for nm, _ in md["gvars"]}
            self._emit_coolant(w, md, md["species"],
                               "coolingF_%s" % self.cid(md["species"]),
                               md["gvars"], fl)
        w("")
        if active:
            w("#define KROME_HAS_COOL_Z 1")
        w("double inline cooling_Z(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if active:
            w("  // above the coolZ.dat @TgasUpperLimit all collisional")
            w("  // rates are cut in legacy, so line cooling vanishes")
            w("  if (Tgas >= %s) return 0.0;" % self.cnum(limit))
            w("  double c = 0.0;")
            for m in active:
                w("  c += coolingZ_%s(Tgas, n, p);" % self.cid(m))
            w("  return c * constants::boltzmann_erg;  // K -> erg")
        else:
            w("  return 0.0;  // no metal line cooling configured")
        w("}")
        w("")
        # one cooling term per coolant-file coolant (solve_cell.hpp slots after
        # the fixed terms: C_FILE0 + i), each zero above its file's
        # @TgasUpperLimit
        w("constexpr int NCOOL_FILE = %d;" % len(files))
        w("#define KROME_COOLFILE_NAMES %s" % "".join(
            ', "%s"' % md["flag"] for md in files))
        w("#define KROME_COOLFILE_ACTIVE %s" % (", true" * len(files)))
        w("KROME_HD inline void cooling_files_(double Tgas, const double* n, "
          "const Params& p, double* out) {")
        w("  (void)Tgas; (void)n; (void)p; (void)out;")
        for q, md in enumerate(files):
            w("  out[%d] = (Tgas >= %s) ? 0.0 : coolingF_%s(Tgas, n, p) "
              "* constants::boltzmann_erg;  // %s (%s)" % (q, self.cnum(md["limit"]),
                                                self.cid(md["species"]),
                                                md["flag"],
                                                Path(md["file"]).name))
        w("}")
        w("")
        self._emit_coolant_table(w, [(self.cool_flag_(m), m, metals[m],
                                      "coolingZ_%s" % self.cid(m), limit)
                                     for m in active]
                                 + [(md["flag"], md["species"], md,
                                     "coolingF_%s" % self.cid(md["species"]),
                                     md["limit"]) for md in files])
        self.emit_zcie(w)

    def _emit_coolant_table(self, w, cl):
        """the line coolants' metadata and evaluator (the query API:
        get_coolant_cooling / get_coolant_populations, legacy
        krome_cooling<metal> and krome_popcool_dump): cl = [(name,
        species, parsed coolant, function, @TgasUpperLimit)], the coolZ.dat
        metals then the thermal.coolant_files coolants. The levels are the file's
        (after thermal.levels), energies [K] and degeneracies."""
        nc = len(cl)
        lev0, E, g, nl = [], [], [], []
        for _, _, md, _, _ in cl:
            n_ = max(md["levels"]) + 1
            lev0.append(len(E))
            nl.append(n_)
            for l in range(n_):
                e_, g_ = md["levels"].get(l, (0.0, 0.0))
                E.append(e_)
                g.append(g_)
        one = lambda v: v if v else ["0"]  # noqa: E731 (no empty arrays)
        w("// the line coolants: the coolZ.dat metals (-cooling=Z), then the")
        w("// thermal.coolant_files ones; their levels as in the file (thermal.levels")
        w("// applied): energies [K], degeneracies")
        w("constexpr int NCOOLANT = %d;" % nc)
        w("inline const char* coolant_names[%d] = {%s};"
          % (max(nc, 1), ", ".join(one(['"%s"' % c[0] for c in cl]))))
        w("constexpr int coolant_species[%d] = {%s};"
          % (max(nc, 1), ", ".join(one(["IDX_%s" % self.cid(c[1])
                                        for c in cl]))))
        w("constexpr int coolant_nlev[%d] = {%s};"
          % (max(nc, 1), ", ".join(one([str(v) for v in nl]))))
        w("constexpr int coolant_lev0[%d] = {%s};"
          % (max(nc, 1), ", ".join(one([str(v) for v in lev0]))))
        w("inline constexpr double coolant_level_E[%d] = {%s};"
          % (max(len(E), 1), ", ".join(one([self.cnum(v) for v in E]))))
        w("inline constexpr double coolant_level_g[%d] = {%s};"
          % (max(len(g), 1), ", ".join(one([self.cnum(v) for v in g]))))
        w("// coolant i at Tgas: its line cooling [erg cm^-3 s^-1] (the")
        w("// terms cooling_Z / cooling_files_ sum), and its level")
        w("// populations into pop[coolant_nlev[i]] unless nullptr (above")
        w("// the file's @TgasUpperLimit: no cooling, all in the ground)")
        w("inline double coolant_eval_(int i, double Tgas, const double* n,")
        w("                            const Params& p, double* pop) {")
        w("  (void)Tgas; (void)n; (void)p; (void)pop;")
        w("  switch (i) {")
        for i, (name, sp, md, fn, lim) in enumerate(cl):
            w("    case %d:  // %s" % (i, name))
            w("      if (Tgas >= %s) {" % self.cnum(lim))
            w("        if (pop != nullptr) {")
            w("          for (int l = 0; l < %d; ++l) pop[l] = 0.0;" % nl[i])
            w("          pop[0] = std::fmax(n[IDX_%s], 0.0);" % self.cid(sp))
            w("        }")
            w("        return 0.0;")
            w("      }")
            w("      return %s(Tgas, n, p, pop) * constants::boltzmann_erg;" % fn)
        w("    default:")
        w("      return 0.0;")
        w("  }")
        w("}")
        w("")

    def _emit_coolz_solver(self, w):
        """the level-population solver of the metal / file coolants"""
        w("// metal line cooling: level-population solvers (legacy")
        w("// mylin2/mylin3 closed forms; partial-pivot Gaussian")
        w("// elimination replaces LAPACK dgesv for >3 levels)")
        w("template <int N>")
        w("inline void coolz_solve_(double (&A)[N][N], double (&B)[N])"
          " {")
        w("  if constexpr (N == 2) {")
        w("    const double iab = B[0] / (A[1][1] - A[1][0]);")
        w("    B[0] = A[1][1] * iab;")
        w("    B[1] = -A[1][0] * iab;")
        w("  } else if constexpr (N == 3) {")
        w("    const double iab = B[0] / (A[1][0]*(A[2][2] - A[2][1])")
        w("        + A[1][1]*(A[2][0] - A[2][2]) "
          "+ A[1][2]*(A[2][1] - A[2][0]));")
        w("    const double c0 = (A[1][2]*A[2][1] - A[1][1]*A[2][2]) "
          "* iab;")
        w("    const double c1 = -(A[1][2]*A[2][0] - A[1][0]*A[2][2]) "
          "* iab;")
        w("    const double c2 = (A[2][0]*A[1][1] - A[1][0]*A[2][1]) "
          "* iab;")
        w("    B[0] = c0; B[1] = c1; B[2] = c2;")
        w("  } else {")
        w("    for (int c = 0; c < N; ++c) {")
        w("      int piv = c;")
        w("      for (int r = c + 1; r < N; ++r)")
        w("        if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) "
          "piv = r;")
        w("      if (piv != c) {")
        w("        for (int j = 0; j < N; ++j) "
          "{ const double t = A[c][j]; A[c][j] = A[piv][j]; "
          "A[piv][j] = t; }")
        w("        const double t = B[c]; B[c] = B[piv]; B[piv] = t;")
        w("      }")
        w("      if (A[c][c] == 0.0) {  // singular: no cooling")
        w("        for (int r = 0; r < N; ++r) B[r] = 0.0;")
        w("        return;")
        w("      }")
        w("      for (int r = c + 1; r < N; ++r) {")
        w("        const double f = A[r][c] / A[c][c];")
        w("        for (int j = c; j < N; ++j) A[r][j] -= f * A[c][j];")
        w("        B[r] -= f * B[c];")
        w("      }")
        w("    }")
        w("    for (int r = N - 1; r >= 0; --r) {")
        w("      double s = B[r];")
        w("      for (int j = r + 1; j < N; ++j) s -= A[r][j] * B[j];")
        w("      B[r] = s / A[r][r];")
        w("    }")
        w("  }")
        w("}")

    def _coolant_reachable(self, md, m):
        """the coolant restricted to the levels its steady state can
        populate: those reachable from the ground level through collisional
        excitation / de-excitation, radiative decay and (with pumping)
        absorption. The others have zero population; a level with no
        transition at all makes the system singular (O2 of
        data/cooling/coolO2.dat: levels 36-47 have no collisional data).
        Legacy dropped every level above the last one excited directly from
        the ground (its nmax), which also cut levels fed by cascades (Fe+ of
        coolZ.dat)."""
        lev = md["levels"]
        edges = {k: set() for k in lev}
        for c in md["coll"]:
            _, u, d, _ = c
            if u in lev and d in lev:
                edges[u].add(d)
                edges[d].add(u)
        pump = self.photo_induced or self.cmb_pumping
        for (u, d, A, wvl) in md["rad"]:
            if u in lev and d in lev:
                edges[u].add(d)
                if pump and (wvl is None or wvl > 0.0):
                    edges[d].add(u)
        seen, todo = {0}, [0]
        while todo:
            k = todo.pop()
            for j in edges.get(k, ()):
                if j not in seen:
                    seen.add(j)
                    todo.append(j)
        keep = sorted(seen)
        if keep == sorted(lev):
            return md
        drop = sorted(set(lev) - seen)
        self.warnings.append("coolant %s: levels %s dropped (no transition "
                             "connects them to the ground level: zero "
                             "population)" % (m, ",".join(map(str, drop))))
        new = {old: i for i, old in enumerate(keep)}
        out = dict(md)
        out["orig_index"] = keep
        out["levels"] = {new[k]: lev[k] for k in keep}
        out["rad"] = [(new[u], new[d], A, wvl) for (u, d, A, wvl)
                      in md["rad"] if u in new and d in new]
        out["coll"] = [[c[0], new[c[1]], new[c[2]], c[3]] for c in md["coll"]
                       if c[1] in new and c[2] in new]
        return out

    def _emit_coolant(self, w, md, m, fname, gvars, localnames):
        """one coolant's level-population cooling function fname
        (Tgas, n, p) [K cm^-3 s^-1]: md the parsed coolant, m its
        network species"""
        kb_eV = constants.boltzmann_eV
        nfile = max(md["levels"]) + 1  # the file's levels (--cool-levels)
        md = self._coolant_reachable(md, m)
        nlev = max(md["levels"]) + 1
        # solved level i is the file's level omap[i] (the unreachable ones
        # are dropped from the system: population 0)
        omap = md.get("orig_index", list(range(nlev)))
        assert sorted(md["levels"]) == list(range(nlev)), \
            "non-contiguous levels for metal %s" % m
        cc = self.cid(m)
        # symbols used by this metal's rate expressions
        used = set()
        asts = {}
        for r, (coll, up, dn, exprs) in enumerate(md["coll"]):
            lst = []
            for cond, expr in exprs:
                ca = self._coolz_expr(cond, localnames) if cond else None
                ea = self._coolz_expr(expr, localnames)
                for a in ([ca] if ca else []) + [ea]:
                    used |= {nd["s"] for nd in knet.walk(a)
                             if nd.get("kind") == "state"}
                lst.append((ca, ea))
            asts[r] = lst
        # the @var shortcuts used, with the ones they use in turn
        gast = {nm: self._coolz_expr(expr, localnames) for nm, expr in gvars}
        changed = True
        while changed:
            changed = False
            for nm, a in gast.items():
                if nm in used:
                    dep = {nd["s"] for nd in knet.walk(a)
                           if nd.get("kind") == "state"} - used
                    if dep:
                        used |= dep
                        changed = True
        self._emit_cftabs(w)  # the tables these rates use
        w("")
        w("// %s: %d levels, %d radiative, %d collisional transitions%s"
          % (m, nlev, len(md["rad"]), len(md["coll"]),
             " (%s)" % Path(md["file"]).name if "file" in md else ""))
        w("// pop: nullptr, or the %d level populations [cm^-3] in the "
          "file's" % nfile)
        w("// level numbering (get_coolant_populations)")
        w("double inline %s(double Tgas, const double* n, "
          "const Params& p, double* pop = nullptr) {" % fname)
        w("  (void)Tgas; (void)p;")
        w("  if (pop != nullptr)")
        w("    for (int l = 0; l < %d; ++l) pop[l] = 0.0;" % nfile)
        w("  if (n[IDX_%s] < 1e-15) {" % cc)
        w("    if (pop != nullptr) pop[0] = std::fmax(n[IDX_%s], 0.0);" % cc)
        w("    return 0.0;")
        w("  }")
        w("  const double invT = 1.0 / Tgas;")
        if "logTgas" in used:
            w("  const double logTgas = std::log10(Tgas);")
        for nm, expr in gvars:
            if nm in used:
                w("  const double %s = %s;" % (nm, self.cexpr(gast[nm])))
        # colliders (clamped >= 0 as legacy)
        colliders = []
        for coll, up, dn, exprs in md["coll"]:
            if coll not in colliders:
                colliders.append(coll)
        for coll in colliders:
            cv = self.cid(coll if coll not in ("H2or", "H2pa") else coll)
            if coll.lower() in ("h2or", "h2pa"):
                # the network's own ortho / para species, else H2 split by
                # the parameter opratio (nH2o_ / nH2p_)
                src = None
                if self.h2_forms():
                    self._need_h2op = True
                    src = "%s(n, p)" % ("nH2o_" if coll.lower() == "h2or"
                                        else "nH2p_")
            else:
                sp = self.match_species_(coll)
                src = ("n[IDX_%s]" % self.cid(sp)
                       if sp is not None else None)
            if src is None:
                self.warnings.append("coolZ collider %s not in network "
                                     "(density 0 in cooling_Z_%s)"
                                     % (coll, m))
                w("  const double coll_%s = 0.0;  "
                  "// collider not in network" % cv)
            else:
                w("  const double coll_%s = std::fmax(%s, 0.0);"
                  % (cv, src))
        nR = len(md["coll"])
        w("  double kz[%d];" % (2 * nR))
        for r, (coll, up, dn, exprs) in enumerate(md["coll"]):
            for ci, (ca, ea) in enumerate(asts[r]):
                if ci == 0:
                    w("  kz[%d] = %s;  // %d->%d %s"
                      % (r, self.cexpr(ea), up, dn, coll))
                else:
                    w("  if (%s) kz[%d] = %s;"
                      % (self.cexpr(ca), r, self.cexpr(ea)))
            w("  kz[%d] = std::fmax(kz[%d], 0.0);" % (r, r))
        w("  // reverse (excitation) rates by detailed balance")
        for r, (coll, up, dn, exprs) in enumerate(md["coll"]):
            gu = md["levels"][up][1]
            gd = md["levels"][dn][1]
            dE = abs(md["levels"][up][0] - md["levels"][dn][0])
            w("  kz[%d] = kz[%d] * %s * std::exp(-%s * invT);"
              "  // %d->%d %s" % (nR + r, r, self.cnum(gu / gd),
                                  self.cnum(dE), dn, up, coll))
        # rate index lookup: (from, to, collider) -> kz index
        ridx = {}
        for r, (coll, up, dn, exprs) in enumerate(md["coll"]):
            ridx.setdefault((up, dn), []).append((r, coll))
            ridx.setdefault((dn, up), []).append((nR + r, coll))
        radmap, radidx = {}, {}
        for li, (u, d, A, wvl) in enumerate(md["rad"]):
            radmap[(u, d)] = A
            radidx[(u, d)] = li
        # --cool-escape=LVG: every radiative term times the line's escape
        # probability bet[li] (the level solve iterated on it); off, the
        # constants as they are
        esc = getattr(self, "cool_escape", "") == "LVG" and bool(md["rad"])

        def arad(u, d):
            A = self.cnum(radmap[(u, d)])
            return "%s * bet[%d]" % (A, radidx[(u, d)]) if esc else A

        def pterm(nm, li):
            return "%s%d * bet[%d]" % (nm, li, li) if esc else \
                "%s%d" % (nm, li)
        # --photo-induced (legacy -usePhotoInduced, corrected): each
        # line with a wavelength gets stimulated emission B_ud J
        # (downward, with A) and absorption B_du J (upward); legacy
        # left the absorption out of the level balance and weighted
        # the absorbed energy with the upper level
        pump = {}
        if self.photo_induced or self.cmb_pumping:
            h_eV, c_cm = constants.planck_eV, constants.clight
            for li, (u, d, A, wvl) in enumerate(md["rad"]):
                # line energy: kB dE of the levels, or h c / lambda
                # when coolZ.dat gives a wavelength [Angstrom] (none
                # for a non-positive one), as legacy
                if wvl is None:
                    E = kb_eV * abs(md["levels"][u][0]
                                    - md["levels"][d][0])
                elif wvl > 0.0:
                    E = h_eV * c_cm / (wvl * 1e-8)
                else:
                    continue
                Bud = A * 0.5 * (h_eV * c_cm) ** 2 / E ** 3
                gu = md["levels"][u][1]
                gd = md["levels"][d][1]
                # stimulated emission / absorption per line: the bin
                # field (--photo-induced) and the CMB blackbody at the
                # per-cell Tcmb, occupation 1/(exp(E/kTcmb) - 1)
                # (--cmb-pumping): B_ud J_cmb = A n_cmb
                sd, ab = [], []
                if self.photo_induced:
                    w("  const double Jl%d = photo_J_(p, %s);  // %d->%d "
                      "line" % (li, self.cnum(E), u, d))
                    sd.append("%s * Jl%d" % (self.cnum(Bud), li))
                    ab.append("%s * Jl%d" % (self.cnum(Bud * gu / gd),
                                             li))
                if self.cmb_pumping:
                    w("  const double oc%d = 1.0 / std::expm1(%s / "
                      "p.Tcmb);  // CMB occupation, %d->%d line"
                      % (li, self.cnum(E / kb_eV), u, d))
                    sd.append("%s * oc%d" % (self.cnum(A), li))
                    ab.append("%s * oc%d" % (self.cnum(A * gu / gd), li))
                w("  const double sd%d = %s, ab%d = %s;"
                  % (li, " + ".join(sd), li, " + ".join(ab)))
                pump[(u, d)] = li
        w("  double A[%d][%d], B[%d];" % (nlev, nlev, nlev))
        if esc:
            nrad = len(md["rad"])
            w("  // thermal.escape = 'LVG': the Sobolev optical depth of each line,")
            w("  // tau = A lambda^3 gu / (8 pi gl) (N / (b n)) (n_l - n_u gl /")
            w("  // gu), N the species column (num2col), b the line width")
            w("  // sqrt(2 kB T / m + vturb^2); escape probability (1 -")
            w("  // exp(-tau)) / tau (1 for tau <= 0, masers). The level solve")
            w("  // is iterated with damped populations Bd (0.3 old + 0.7")
            w("  // new, after 5 plain steps, as RADEX) until they change by")
            w("  // < 1e-10 of the species density, at most 100 times; the")
            w("  // damped populations and their bet[] are the result")
            w("  const double Ncol = num2col_(n[IDX_%s], Tgas, n, p);" % cc)
            w("  const double bD = std::sqrt(2.0 * constants::boltzmann_erg * Tgas / "
              "species_mass_g[IDX_%s] + p.vturb * p.vturb);" % cc)
            w("  const double tpre = Ncol / (bD * n[IDX_%s]);" % cc)
            w("  double bet[%d], Bd[%d];" % (nrad, nlev))
            w("  for (int l = 0; l < %d; ++l) bet[l] = 1.0;" % nrad)
            w("  for (int it = 0; it < 100; ++it) {")
        w("  for (int r = 0; r < %d; ++r) { B[r] = 0.0; "
          "for (int c = 0; c < %d; ++c) A[r][c] = 0.0; }"
          % (nlev, nlev))
        w("  for (int c = 0; c < %d; ++c) A[0][c] = 1.0;  "
          "// conservation row" % nlev)
        w("  B[0] = n[IDX_%s];" % cc)

        def coll_var(coll):
            return "coll_" + self.cid(coll)

        for k in range(1, nlev):
            for i in range(nlev):
                terms = []
                if i == k:  # depopulation of level k
                    for j in range(nlev):
                        if j == k:
                            continue
                        for (r, coll) in ridx.get((k, j), []):
                            terms.append("- kz[%d]*%s"
                                         % (r, coll_var(coll)))
                        if (k, j) in radmap:
                            terms.append("- %s" % arad(k, j))
                            if (k, j) in pump:
                                terms.append("- %s"
                                             % pterm("sd", pump[(k, j)]))
                        if (j, k) in pump:  # absorption k -> j
                            terms.append("- %s" % pterm("ab", pump[(j, k)]))
                else:  # population of level k from level i
                    for (r, coll) in ridx.get((i, k), []):
                        terms.append("+ kz[%d]*%s" % (r, coll_var(coll)))
                    if (i, k) in radmap:
                        terms.append("+ %s" % arad(i, k))
                        if (i, k) in pump:
                            terms.append("+ %s" % pterm("sd", pump[(i, k)]))
                    if (k, i) in pump:  # absorption i -> k
                        terms.append("+ %s" % pterm("ab", pump[(k, i)]))
                if terms:
                    expr = " ".join(terms)
                    if expr.startswith("+ "):
                        expr = expr[2:]
                    w("  A[%d][%d] = %s;" % (k, i, expr))
        w("  coolz_solve_<%d>(A, B);" % nlev)
        w("  // sanitize small negative populations (legacy: 1e-40);")
        w("  // large ones are clamped too instead of aborting")
        w("  for (int r = 0; r < %d; ++r)" % nlev)
        w("    if (B[r] < 0.0)")
        w("      B[r] = (std::fabs(B[r] / n[IDX_%s]) > 1e-10) "
          "? 0.0 : 1e-40;" % cc)
        if esc:
            w("  double chg = 0.0;  // damped populations, then the bet[]")
            w("  for (int r = 0; r < %d; ++r) {" % nlev)
            w("    if (it >= 5) {")
            w("      chg = std::fmax(chg, std::fabs(B[r] - Bd[r]));")
            w("      Bd[r] = 0.3 * Bd[r] + 0.7 * B[r];")
            w("    } else {")
            w("      Bd[r] = B[r];")
            w("    }")
            w("  }")
            for li, (u, d, A, wvl) in enumerate(md["rad"]):
                dEK = abs(md["levels"][u][0] - md["levels"][d][0])
                if dEK <= 0.0:
                    continue  # no line: stays 1
                gu, gd = md["levels"][u][1], md["levels"][d][1]
                lam = constants.hc_k / dEK  # [cm], hc / (kB dE)
                C = A * lam ** 3 * gu / (8.0 * math.pi * gd)
                w("  {  // %d->%d" % (u, d))
                w("    const double tau = %s * tpre * (Bd[%d] - Bd[%d] * %s);"
                  % (self.cnum(C), d, u, self.cnum(gd / gu)))
                w("    bet[%d] = tau > 1e-5 ? -std::expm1(-tau) / tau"
                  " : (tau > 0.0 ? 1.0 - 0.5 * tau : 1.0);" % li)
                w("  }")
            w("  if (it >= 5 && chg < 1e-10 * n[IDX_%s]) break;" % cc)
            w("  }")
            w("  for (int r = 0; r < %d; ++r) B[r] = Bd[r];" % nlev)
        if omap == list(range(nlev)):
            w("  if (pop != nullptr)")
            w("    for (int r = 0; r < %d; ++r) pop[r] = B[r];" % nlev)
        else:
            w("  if (pop != nullptr) {  // back to the file's numbering")
            w("    static constexpr int omap[%d] = {%s};"
              % (nlev, ", ".join(map(str, omap))))
            w("    for (int r = 0; r < %d; ++r) pop[omap[r]] = B[r];" % nlev)
            w("  }")
        cool = []
        for (u, d, A, wvl) in md["rad"]:
            # legacy roundtrips dE through eV with its kB value
            de_K = kb_eV * abs(md["levels"][u][0]
                               - md["levels"][d][0]) / kb_eV
            b = " * bet[%d]" % radidx[(u, d)] if esc else ""
            if (u, d) in pump:
                li = pump[(u, d)]
                cool.append("(B[%d] * (%s + sd%d) - B[%d] * ab%d)%s * %s"
                            % (u, self.cnum(A), li, d, li, b,
                               self.cnum(de_K)))
            else:
                cool.append("B[%d] * %s%s * %s"
                            % (u, self.cnum(A), b, self.cnum(de_K)))
        w("  return %s;  // [K cm^-3 s^-1]" % (" + ".join(cool)))
        w("}")

    # -- ZCIE: Shen et al. 2013 metal cooling/heating tables -------------------

    def load_zcie_3d(self):
        """coolZ_CIE2012.dat (CLOUDY 10 + Haardt & Madau 2012 UV): lines of
        'i j k  log10(T) log10(nH) z  log10(cool) log10(heat)' on a uniform
        grid; returns axes plus both flattened tables."""
        path = extract_cie_tables(self.data_dir.parent) / "coolZ_CIE2012.dat"
        x1, x2, x3, yc, yh = {}, {}, {}, {}, {}
        for line in path.read_text().splitlines():
            s = line.split("#")[0].split()
            if len(s) < 8:
                continue
            i, j, k = int(s[0]), int(s[1]), int(s[2])
            x1[i], x2[j], x3[k] = float(s[3]), float(s[4]), float(s[5])
            yc[(i, j, k)] = float(s[6])
            yh[(i, j, k)] = float(s[7])
        n1, n2, n3 = max(x1), max(x2), max(x3)
        if len(yc) != n1 * n2 * n3:
            raise ValueError("coolZ_CIE2012.dat is not a full grid")
        ax = ([x1[i] for i in range(1, n1 + 1)],
              [x2[j] for j in range(1, n2 + 1)],
              [x3[k] for k in range(1, n3 + 1)])
        flats = []
        for y in (yc, yh):
            flats.append([y[(i + 1, j + 1, k + 1)] for i in range(n1)
                          for j in range(n2) for k in range(n3)])
        return ax[0], ax[1], ax[2], flats[0], flats[1]

    def load_u2d(self, relpath):
        """legacy init_anytab2D file: an 'n1, n2' header (in either order:
        legacy takes the sizes from the Fortran arrays, and escape_H2.dat
        lists them reversed), then one block of '(x, y, value)' rows per x;
        returns (xs, ys, flat values [x][y])."""
        path = self.data_dir.parent / relpath
        if relpath.startswith("coolZ_CIE"):
            path = extract_cie_tables(self.data_dir.parent) / relpath
        head = None
        xs, ys, zs = [], [], []
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if head is None:
                head = sorted(int(v) for v in s.split(","))
                continue
            v = [float(tk) for tk in s.split()]
            if not xs or v[0] != xs[-1]:
                xs.append(v[0])
                zs.append([])
            if len(xs) == 1:
                ys.append(v[1])
            zs[-1].append(v[2])
        if sorted((len(xs), len(ys))) != head or \
                any(len(r) != len(ys) for r in zs):
            raise ValueError("%s is not a %sx%s grid" % (relpath, *head))
        return xs, ys, [v for row in zs for v in row]

    def load_zcie_2d(self):
        """coolZ_CIE2012NOUV.dat (CLOUDY 10, no UV): 'nx, ny' header, then
        nx blocks of ny '(log10 nH, log10 T, log10 cool)' rows."""
        return self.load_u2d("coolZ_CIE2012NOUV.dat")

    def emit_zcie(self, w):
        """cooling_ZCIE()/heating_ZCIE(): legacy -cooling=Z_CIE /
        Z_CIENOUV — Shen et al. 2013 (ApJ 765, 89) CLOUDY tables scaled as
        10^tab * n_Hnuclei^2 * total_Z. The 2012 variant (with the HM12 UV
        background) is trilinear in (log T, log nH, z) with a heating
        table; the NOUV variant is bilinear in (log nH, log T) with the
        legacy fit_anytab2D edge behavior (linear extrapolation)."""
        mode = self.zcie
        if mode == "2012":
            ax1, ax2, ax3, fc, fh = self.load_zcie_3d()
            self._emit_interp3_body(w)
            self._emit_cool3d_table(
                w, "zcie_cool", (ax1, ax2, ax3, fc),
                "ZCIE cooling, coolZ_CIE2012.dat (%dx%dx%d), Shen+13"
                % (len(ax1), len(ax2), len(ax3)))
            self._emit_cool3d_table(
                w, "zcie_heat", (ax1, ax2, ax3, fh),
                "ZCIE heating, coolZ_CIE2012.dat (%dx%dx%d), Shen+13"
                % (len(ax1), len(ax2), len(ax3)))
            w("#define KROME_HAS_COOL_ZCIE 1")
            for nm, tab in (("cooling", "zcie_cool"), ("heating",
                                                       "zcie_heat")):
                w("double inline %s_ZCIE(double Tgas, const double* n, "
                  "const Params& p) {" % nm)
                w("  const double cH = Hnuclei_(n);")
                w("  if (cH < 1e-20) return 0.0;")
                w("  const double x = mol_cool_interp_(%s_tab," % tab)
                w("      std::log10(Tgas), std::log10(cH), p.zredshift);")
                w("  return std::pow(10.0, x) * cH * cH * p.total_Z;")
                w("}")
        elif mode == "nouv":
            xs, ys, flat = self.load_zcie_2d()
            w("// ZCIE cooling (no UV), coolZ_CIE2012NOUV.dat (%dx%d "
              "in log nH, log T), Shen+13" % (len(xs), len(ys)))
            w("inline const double zcie_x[%d] = {%s};"
              % (len(xs), ", ".join(self.cnum(v) for v in xs)))
            w("inline const double zcie_t[%d] = {%s};"
              % (len(ys), ", ".join(self.cnum(v) for v in ys)))
            self.emit_data_array(w, "zcie_z", flat,
                                 "ZCIE NOUV table (log10 cool)")
            w("// legacy fit_anytab1D/2D: uniform grids, clamped cell")
            w("// index, linear EXTRAPOLATION beyond the edges")
            w("inline double zcie_fit1d_(const double* x, const double* z, "
              "int nx, double xmul, double xx) {")
            w("  int i = (int)((xx - x[0]) * xmul);")
            w("  i = (i < 0) ? 0 : ((i > nx - 2) ? nx - 2 : i);")
            w("  const double pf = (xx - x[i]) * xmul;")
            w("  return pf * (z[i+1] - z[i]) + z[i];")
            w("}")
            w("#define KROME_HAS_COOL_ZCIE 1")
            w("double inline cooling_ZCIE(double Tgas, const double* n, "
              "const Params& p) {")
            nx, ny = len(xs), len(ys)
            xmul = 1.0 / (xs[1] - xs[0])
            ymul = 1.0 / (ys[1] - ys[0])
            w("  const double cH = Hnuclei_(n);")
            w("  if (cH < 1e-20) return 0.0;")
            w("  const double xx = std::log10(cH), yy = std::log10(Tgas);")
            w("  const double ymul = %s;" % self.cnum(ymul))
            w("  int j = (int)((yy - zcie_t[0]) * ymul);")
            w("  j = (j < 0) ? 0 : ((j > %d) ? %d : j);" % (ny - 2, ny - 2))
            w("  double zcol_l[%d], zcol_r[%d];" % (nx, nx))
            w("  for (int i = 0; i < %d; ++i) {" % nx)
            w("    zcol_l[i] = zcie_z[i*%d + j];" % ny)
            w("    zcol_r[i] = zcie_z[i*%d + j + 1];" % ny)
            w("  }")
            w("  const double vl = zcie_fit1d_(zcie_x, zcol_l, %d, %s, xx);"
              % (nx, self.cnum(xmul)))
            w("  const double vr = zcie_fit1d_(zcie_x, zcol_r, %d, %s, xx);"
              % (nx, self.cnum(xmul)))
            w("  const double x = (yy - zcie_t[j]) * ymul * (vr - vl) + vl;")
            w("  return std::pow(10.0, x) * cH * cH * p.total_Z;")
            w("}")
            w("double inline heating_ZCIE(double, const double*, "
              "const Params&) { return 0.0; }  // NOUV has no heating table")
        if not mode:
            w("double inline cooling_ZCIE(double, const double*, "
              "const Params&) { return 0.0; }")
            w("double inline heating_ZCIE(double, const double*, "
              "const Params&) { return 0.0; }")
        w("")

    # -- dust tables (legacy -useDustTabs) -------------------------------------

    def load_u3d(self, relpath):
        """legacy init_anytab3D file: an 'nx, ny, nz' header, then rows
        '(z, x, y, value)' with z outermost (blank lines ignored); returns
        (xs, ys, zs, flat values [z][x][y])."""
        path = self.data_dir.parent / relpath
        head, rows = None, []
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if head is None:
                head = [int(v) for v in s.split(",")]
                continue
            rows.append([float(v) for v in s.split()])
        zs = sorted({r[0] for r in rows})
        xs = sorted({r[1] for r in rows})
        ys = sorted({r[2] for r in rows})
        if [len(xs), len(ys), len(zs)] != head or len(rows) != \
                len(xs) * len(ys) * len(zs):
            raise ValueError("%s is not a %sx%sx%s grid" % (relpath, *head))
        val = {(r[0], r[1], r[2]): r[3] for r in rows}
        flat = [val[(z, x, y)] for z in zs for x in xs for y in ys]
        return xs, ys, zs, flat

    def _emit_dtab_helpers(self, w):
        w("// legacy fit_anytab2D: uniform-grid bilinear, clamped cell")
        w("// index = linear extrapolation beyond every edge")
        w("KROME_HD inline double dtab_fit2d_(const double* xs, "
          "const double* ys, const double* z, int nx, int ny,")
        w("                          double xmul, double ymul, "
          "double xx, double yy) {")
        w("  int i = (int)((xx - xs[0]) * xmul);")
        w("  i = i < 0 ? 0 : (i > nx - 2 ? nx - 2 : i);")
        w("  int j = (int)((yy - ys[0]) * ymul);")
        w("  j = j < 0 ? 0 : (j > ny - 2 ? ny - 2 : j);")
        w("  const double pf = (xx - xs[i]) * xmul;")
        w("  const double zl = pf * (z[(i+1)*ny + j] - z[i*ny + j])"
          " + z[i*ny + j];")
        w("  const double zr = pf * (z[(i+1)*ny + j+1] - "
          "z[i*ny + j+1]) + z[i*ny + j+1];")
        w("  return (yy - ys[j]) * ymul * (zr - zl) + zl;")
        w("}")
        w("// sign-aware geometric interpolation for linear-stored tables")
        w("// (legacy fit_anytab1D/2D_linlog, ported verbatim)")
        w("KROME_HD inline double dtab_fit1d_ll_(const double* xs, "
          "const double* z, int nx, int ny, int j,")
        w("                             double xmul, double xx) {")
        w("  int i = (int)((xx - xs[0]) * xmul);")
        w("  i = i < 0 ? 0 : (i > nx - 2 ? nx - 2 : i);")
        w("  const double pf = (xx - xs[i]) * xmul;")
        w("  double z1 = z[i*ny + j], z2 = z[(i+1)*ny + j];")
        w("  if (z1 < 0.0 && z2 < 0.0) {")
        w("    z1 = std::log10(-z1); z2 = std::log10(-z2);")
        w("    return -std::pow(10.0, pf * (z2 - z1) + z1);")
        w("  }")
        w("  if (z1 > 0.0 && z2 > 0.0) {")
        w("    z1 = std::log10(z1); z2 = std::log10(z2);")
        w("    return std::pow(10.0, pf * (z2 - z1) + z1);")
        w("  }")
        w("  return pf * (z2 - z1) + z1;")
        w("}")
        w("KROME_HD inline double dtab_fit2d_ll_(const double* xs, "
          "const double* ys, const double* z, int nx, int ny,")
        w("                             double xmul, double ymul, "
          "double xx, double yy) {")
        w("  int j = (int)((yy - ys[0]) * ymul);")
        w("  j = j < 0 ? 0 : (j > ny - 2 ? ny - 2 : j);")
        w("  const double zl = dtab_fit1d_ll_(xs, z, nx, ny, j, xmul, xx);")
        w("  const double zr = dtab_fit1d_ll_(xs, z, nx, ny, j + 1, xmul, "
          "xx);")
        w("  return (yy - ys[j]) * ymul * (zr - zl) + zl;")
        w("}")

    def emit_dust_tables(self, w):
        """table_Tdust() / cooling_dustTab() / krate_H2dust_table():
        legacy dust tables for a fixed radiation field.
          HM2012, MC  2D tables in (log10 ntot, log10 Tgas)
                      (data/dust_tables/dust_table_<field>_*.dat)
          3D          3D tables in (log10 Av, log10 ntot, log10 Tgas)
                      (data/dust_tables/dust_table_3D/), Av the per-cell
                      parameter (legacy -dustTabs=...,3D, whose file names
                      never matched these files)
          PHOTO3D     the 3D tables at the Av that reproduces the dust
                      absorption of the photon-bin field (legacy
                      -dustTabs=...,Photo3D, which never worked: the setup
                      was never called, and its Av weight was not
                      normalized by the grid step, fixed here)
        Tdust = 10^fit and the H2 rate mu 10^fit ntot [s^-1] (log10
        stored, as legacy); cooling = mu coolFit ntot^2. The cool tables:
        MC and 3D store linear coefficients (legacy sign-aware linlog in
        ntot, linear in Tgas and log Av); the HM2012 2D table stores log10
        values, inconsistent with the linlog reader legacy applies to it
        (user-confirmed), so it is read as 10^fit. Without the option the
        functions are stubs (table_Tdust falls back to p.Tdust)."""
        mode = self.dust_tables
        uses_h2 = "krate_H2dust_table" in self.net.externs
        if uses_h2 and not mode:
            raise ValueError("krate_H2dust_table() needs dust.tables")
        if mode in ("HM2012", "MC"):
            self._emit_dust_tables_2d(w, mode, uses_h2)
        elif mode in ("3D", "PHOTO3D"):
            self._emit_dust_tables_3d(w, mode, uses_h2)
        else:
            w("// no dust tables configured (dust.tables): the Tdust")
            w("// parameter, or the binned-dust mean when a hook is live")
            w("// (the macro is suspended: its expansion IS this function)")
            w("#undef Tdust")
            w("template <class P>")
            w("KROME_HD double table_Tdust(double, const double*, const P& p) {")
            w("  if (p.hk.dust.nb > 0) return dustbins_Td_mean_(p);")
            w("  return p.Tdust;")
            w("}")
            w(self._tdust_define())
            w("KROME_HD double inline cooling_dustTab(double, const double*, "
              "const Params&) { return 0.0; }")
        w("")

    def _emit_dust_tables_2d(self, w, mode, uses_h2):
        xs, ys, td = self.load_u2d(
            "dust_tables/dust_table_%s_Tdust.dat" % mode)
        xs2, ys2, cl = self.load_u2d(
            "dust_tables/dust_table_%s_cool.dat" % mode)
        assert xs == xs2 and ys == ys2, "dust table grids differ"
        h2 = None
        if uses_h2:
            xs3, ys3, h2 = self.load_u2d(
                "dust_tables/dust_table_%s_H2.dat" % mode)
            assert xs == xs3 and ys == ys3, "dust H2 table grid differs"
        nx, ny = len(xs), len(ys)
        xmul = 1.0 / (xs[1] - xs[0])
        ymul = 1.0 / (ys[1] - ys[0])
        self._emit_dtab_helpers(w)
        w("// dust tables %s (%dx%d in log10 ntot, log10 Tgas)"
          % (mode, nx, ny))
        arrs = [("x", xs), ("t", ys), ("td", td), ("cool", cl)]
        if h2 is not None:
            arrs.append(("h2", h2))
        for nm, arr in arrs:
            self.emit_data_array(w, "dusttab_%s" % nm, arr,
                                 "dust table %s" % mode)
        fit = ("%d, %d, %s, %s, std::log10(nt), std::log10(Tgas)"
               % (nx, ny, self.cnum(xmul), self.cnum(ymul)))
        w("#define KROME_HAS_DUST_TABLES 1")
        w("template <class P>")
        w("KROME_HD double table_Tdust(double Tgas, const double* n, const P& p) {")
        w("  (void)p;")
        w("  if (p.hk.dust.nb > 0) return dustbins_Td_mean_(p);")
        w("  const double nt = ntot_(n);")
        w("  if (nt <= 0.0) return 1.0;  // legacy default")
        w("  return std::pow(10.0, dtab_fit2d_(dusttab_x, dusttab_t, "
          "dusttab_td, %s));" % fit)
        w("}")
        log_stored = max(abs(v) for v in cl) > 1.0
        w("KROME_HD double inline cooling_dustTab(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)p;")
        w("  const double nt = ntot_(n);")
        w("  if (nt <= 0.0) return 0.0;")
        w("  const double mu = rho_(n) / nt / constants::p_mass;")
        if log_stored:
            w("  // %s cool table stores log10 values -> 10^fit" % mode)
            w("  return mu * std::pow(10.0, dtab_fit2d_(dusttab_x, dusttab_t,"
              " dusttab_cool, %s)) * nt * nt;" % fit)
        else:
            w("  // %s cool table stores linear values -> legacy linlog"
              % mode)
            w("  return mu * dtab_fit2d_ll_(dusttab_x, dusttab_t, "
              "dusttab_cool, %s) * nt * nt;" % fit)
        w("}")
        if h2 is not None:
            self._emit_h2dust_table(w, "dtab_fit2d_(dusttab_x, dusttab_t, "
                                       "dusttab_h2, %s)" % fit)

    def _emit_h2dust_table(self, w, fit_expr):
        w("// H2 formation on dust from the tables (legacy -dustTabs=H2):")
        w("// nH2dust = mu n_H 10^fit ntot, the rate of H + {H} -> H2")
        w("KROME_HD inline double krate_H2dust_table(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)p;")
        w("  const double nt = ntot_(n);")
        w("  if (nt <= 0.0) return 0.0;")
        w("  const double mu = rho_(n) / nt / constants::p_mass;")
        w("  return mu * std::pow(10.0, %s) * nt;" % fit_expr)
        w("}")

    def _emit_dust_tables_3d(self, w, mode, uses_h2):
        base = "dust_tables/dust_table_3D/dust_table_%s_3D.dat"
        xs, ys, zs, td = self.load_u3d(base % "Tdust")
        g2 = self.load_u3d(base % "cool")
        assert g2[:3] == (xs, ys, zs), "3D dust table grids differ"
        cl = g2[3]
        h2 = None
        if uses_h2:
            g3 = self.load_u3d(base % "H2")
            assert g3[:3] == (xs, ys, zs), "3D dust H2 table grid differs"
            h2 = g3[3]
        nx, ny, nz = len(xs), len(ys), len(zs)
        xmul = 1.0 / (xs[1] - xs[0])
        ymul = 1.0 / (ys[1] - ys[0])
        zmul = 1.0 / (zs[1] - zs[0])
        self._emit_dtab_helpers(w)
        w("// 3D dust tables (%dx%dx%d in log10 Av, log10 ntot, log10 Tgas;"
          % (nz, nx, ny))
        w("// data/dust_tables/dust_table_3D), flat [Av][ntot][Tgas]")
        arrs = [("x", xs), ("t", ys), ("z", zs), ("td", td), ("cool", cl)]
        if h2 is not None:
            arrs.append(("h2", h2))
        for nm, arr in arrs:
            self.emit_data_array(w, "dust3_%s" % nm, arr, "3D dust table")
        if mode == "3D":
            w("// log10 Av: the per-cell parameter; legacy fit_anytab3D")
            w("// extrapolates linearly in log Av beyond the grid")
            w("template <class P>")
            w("KROME_HD inline double dust3_lav_(const P& p) {")
            w("  return std::log10(std::fmax(p.Av, 1e-40));")
            w("}")
            clampz = False
        else:
            self._emit_photo3d_av(w)
            clampz = True
        w("// interpolation in log Av between the two neighbouring slices")
        w("// (bilinear, or linlog for the linear cool table, in each)")
        w("template <bool LL>")
        w("KROME_HD inline double dust3_fit_(const double* f, double xx, double yy, "
          "double zz) {")
        w("  int k = (int)((zz - dust3_z[0]) * %s);" % self.cnum(zmul))
        w("  k = k < 0 ? 0 : (k > %d ? %d : k);" % (nz - 2, nz - 2))
        w("  double wz = (zz - dust3_z[k]) * %s;" % self.cnum(zmul))
        if clampz:
            w("  wz = wz < 0.0 ? 0.0 : (wz > 1.0 ? 1.0 : wz);  // Photo3D "
              "clamps")
        w("  const double* f0 = f + k * %d;" % (nx * ny))
        w("  const double* f1 = f0 + %d;" % (nx * ny))
        args = ("%d, %d, %s, %s, xx, yy" % (nx, ny, self.cnum(xmul),
                                            self.cnum(ymul)))
        w("  const double a = LL ? dtab_fit2d_ll_(dust3_x, dust3_t, f0, %s)"
          % args)
        w("                     : dtab_fit2d_(dust3_x, dust3_t, f0, %s);"
          % args)
        w("  const double b = LL ? dtab_fit2d_ll_(dust3_x, dust3_t, f1, %s)"
          % args)
        w("                     : dtab_fit2d_(dust3_x, dust3_t, f1, %s);"
          % args)
        w("  return a + wz * (b - a);")
        w("}")
        lav = "dust3_lav_(p)" if mode == "3D" else "dust3_lav_photo_(p)"
        w("#define KROME_HAS_DUST_TABLES 1")
        w("template <class P>")
        w("KROME_HD double table_Tdust(double Tgas, const double* n, const P& p) {")
        w("  (void)p;")
        w("  if (p.hk.dust.nb > 0) return dustbins_Td_mean_(p);")
        w("  const double nt = ntot_(n);")
        w("  if (nt <= 0.0) return 1.0;  // legacy default")
        w("  return std::pow(10.0, dust3_fit_<false>(dust3_td, "
          "std::log10(nt), std::log10(Tgas), %s));" % lav)
        w("}")
        w("KROME_HD double inline cooling_dustTab(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)p;")
        w("  const double nt = ntot_(n);")
        w("  if (nt <= 0.0) return 0.0;")
        w("  const double mu = rho_(n) / nt / constants::p_mass;")
        w("  return mu * dust3_fit_<true>(dust3_cool, std::log10(nt), "
          "std::log10(Tgas), %s) * nt * nt;" % lav)
        w("}")
        if h2 is not None:
            self._emit_h2dust_table(
                w, "dust3_fit_<false>(dust3_h2, std::log10(nt), "
                   "std::log10(Tgas), %s)" % lav)

    def _emit_photo3d_av(self, w):
        """legacy get_int_JQabs_tab + convert_to_Av: the dust absorption of
        the photon field sum_b J_b I_b (I_b the absorption table integrated
        over bin b with overlap weights), mapped to Av by the Av-Lambda
        table (log-linear in the absorption, clamped at the ends)"""
        d = "dust_tables/dust_table_3D/"
        rows = [[float(v) for v in ln.split()] for ln in
                (self.data_dir.parent / (d + "dust_table_absorption.dat"))
                .read_text().splitlines() if ln.strip()]
        avl = [[float(v) for v in ln.split()] for ln in
               (self.data_dir.parent / (d + "dust_table_Av_Lambda.dat"))
               .read_text().splitlines() if ln.strip()]
        self.emit_data_array(w, "dust3_abs_el", [r[1] for r in rows],
                             "absorption table: bin lower edge [erg]")
        self.emit_data_array(w, "dust3_abs_eu", [r[2] for r in rows],
                             "upper edge [erg]")
        self.emit_data_array(w, "dust3_abs_g", [r[3] for r in rows],
                             "tGamma")
        self.emit_data_array(w, "dust3_av_tab", [r[0] for r in avl],
                             "Av of the Av-Lambda table")
        self.emit_data_array(w, "dust3_lam_tab", [r[1] for r in avl],
                             "absorption of the Av-Lambda table")
        nr, na = len(rows), len(avl)
        w("// per-bin absorption integrand, cached per bin layout")
        w("struct Dust3AbsCache_ {")
        w("  int nb = -1;")
        w("  double e0 = 0.0, e1 = 0.0;")
        w("  std::vector<double> I;")
        w("};")
        w("#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)")
        w("inline thread_local Dust3AbsCache_ dust3_abs_cache_;")
        w("#endif")
        w("KROME_HD inline double dust3_lav_photo_(const Params& p) {")
        w("  double lam = 0.0;")
        w("  const int B = p.hk.nbins;")
        w("#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)")
        w("  // device code: the same integrands, computed in place")
        w("  for (int b = 0; b < B; ++b) {")
        w("    const double eL = p.hk.elo[b] * constants::eV_to_erg;")
        w("    const double eR = p.hk.ehi[b] * constants::eV_to_erg;")
        w("    double Ib = 0.0;")
        w("    for (int j = 0; j < %d; ++j) {" % nr)
        w("      const double el = dust3_abs_el[j], eu = dust3_abs_eu[j];")
        w("      if (el > eR || eu < eL) continue;")
        w("      const double lo = el > eL ? el : eL;")
        w("      const double hi = eu < eR ? eu : eR;")
        w("      Ib += dust3_abs_g[j] * (hi - lo) / (eu - el);")
        w("    }")
        w("    lam += photo_Jbin_(p, b) * Ib;")
        w("  }")
        w("  if (B > 0) lam *= constants::eV_to_erg;  // eV -> erg")
        w("#else")
        w("  if (B > 0) {")
        w("    Dust3AbsCache_& C = dust3_abs_cache_;")
        w("    if (C.nb != B || C.e0 != p.hk.elo[0] ||")
        w("        C.e1 != p.hk.ehi[B - 1]) {")
        w("      C.nb = B; C.e0 = p.hk.elo[0]; "
          "C.e1 = p.hk.ehi[B - 1];")
        w("      C.I.assign(B, 0.0);")
        w("      for (int b = 0; b < B; ++b) {")
        w("        const double eL = p.hk.elo[b] * constants::eV_to_erg;")
        w("        const double eR = p.hk.ehi[b] * constants::eV_to_erg;")
        w("        for (int j = 0; j < %d; ++j) {" % nr)
        w("          const double el = dust3_abs_el[j], eu = dust3_abs_eu[j];")
        w("          if (el > eR || eu < eL) continue;")
        w("          const double lo = el > eL ? el : eL;")
        w("          const double hi = eu < eR ? eu : eR;")
        w("          C.I[b] += dust3_abs_g[j] * (hi - lo) / (eu - el);")
        w("        }")
        w("      }")
        w("    }")
        w("    for (int b = 0; b < B; ++b) lam += photo_Jbin_(p, b) * C.I[b];")
        w("    lam *= constants::eV_to_erg;  // eV -> erg")
        w("  }")
        w("#endif")
        w("  double Av;")
        w("  if (lam >= dust3_lam_tab[0]) Av = dust3_av_tab[0];")
        w("  else if (lam <= dust3_lam_tab[%d]) Av = dust3_av_tab[%d];"
          % (na - 1, na - 1))
        w("  else {")
        w("    int i = 1;")
        w("    while (!(lam > dust3_lam_tab[i])) ++i;")
        w("    const double wgt = (std::log(lam) - std::log(dust3_lam_tab[i]))")
        w("        / (std::log(dust3_lam_tab[i - 1]) - "
          "std::log(dust3_lam_tab[i]));")
        w("    Av = dust3_av_tab[i] * (1.0 - wgt) + dust3_av_tab[i - 1] * wgt;")
        w("  }")
        w("  return std::log10(Av);")
        w("}")

    # -- enthalpy heating/cooling (legacy -cooling=DH / -heating=DH) -----------

    def load_thermo(self):
        """data/thermo_burcat.toml (Burcat's NASA 7-coefficient polynomials
        of the networks' species, tools/maintainer/tables/
        gen_thermo_burcat.py): species name (as in the networks) -> the 17
        numbers [Tlow, Tcommon, Thigh, 7 upper-range, 7 lower-range
        coefficients]"""
        if not hasattr(self, "_thermo"):
            import tomlio
            for base in (self.data_dir.parent, self.data_dir,
                         self.data_dir.parent.parent):
                path = base / "thermo_burcat.toml"
                if path.exists():
                    break
            else:
                raise FileNotFoundError("thermo_burcat.toml not found")
            sp = tomlio.load(path).get("species", {})
            self._thermo = {k: list(v["T"]) + list(v["high"]) + list(v["low"])
                            for k, v in sp.items()}
        return self._thermo

    def recomb_reactions(self):
        """the electron recombinations of atomic ions, X^q+ + e- ->
        X^(q-1)+ (radiative and dielectronic, as the network writes
        them): the processes whose captured electrons the ATOMIC
        recombination cooling counts. Not: molecular ions (dissociative
        recombination), three-body recombination, recombination on dust
        grains, photo kinds."""
        net, out = self.net, []
        for rx in net.reactions:
            if rx["rate"]["kind"] != "expr":
                continue
            R, P = rx["reactants"], rx["products"]
            if len(R) != 2 or len(P) != 1:
                continue
            ions = [r["sp"] for r in R if r["sp"] != "e-"]
            if len(ions) != 1 or len(R) - len(ions) != 1:
                continue
            ion, prod = ions[0], P[0]["sp"]
            si, sp = net.species[ion], net.species[prod]
            if si["opaque"] or sp["opaque"] or si["phase"] != "gas" \
                    or sp["phase"] != "gas" or self._natoms(ion) != 1:
                continue
            if si["composition"] != sp["composition"] \
                    or sp["charge"] != si["charge"] - 1 or si["charge"] < 1:
                continue
            if any(nd.get("kind") in ("extern", "table")
                   and "dust" in str(nd.get("f", "")).lower()
                   or nd.get("s") == "Tdust"
                   for nd in self.deep_walk(rx["rate"]["ast"])):
                continue
            out.append(rx)
        return out

    def emit_recomb_cooling(self, w):
        """cooling_recomb_(): the recombination cooling of thermal.cooling = ['ATOMIC']
        from the network's own recombination rates (recomb_reactions).
        Captured electrons take their kinetic energy with them: for a
        Maxwellian rate k(T) = <sigma v>, <sigma v E> = kB T (3/2 k +
        T dk/dT) (exact for any electron-capture cross section; for a
        dielectronic fit k ~ T^-3/2 exp(-E0/kT) it is E0 k). T dk/dT by a
        central difference (relative step 1e-4) of the piece rates()
        selects at T, 0 where rates() clamps the temperature, at the rate
        temperature clamp(T, Tcmb, 1e9). Each term is floored at 0 (a
        fit steeper than T^-3/2 would give a negative energy). Legacy
        used the Cen (1992) fits for H+, He+ (with dielectronic) and He++,
        whatever the network's rates."""
        rxs = self.recomb_reactions() if self.thermal.get("cool_atomic") \
            else []
        idx = {rx["id"]: j for j, rx in enumerate(rxs)}
        if rxs:
            w("// the recombination rates of the ATOMIC cooling at T, and")
            w("// T dk/dT of the piece active at T (0 where the rate is")
            w("// clamped): a difference across a piece boundary would turn")
            w("// the jump between two fits into a spike")
            w("inline void recomb_rates_(double Tgas, const double* n, "
              "const Params& p, double* kr, double* tdk) {")
            w("  (void)n; (void)p;")
            w("  const double h = 1e-4;")

            def body(rx, tv):
                j = idx[rx["id"]]

                def e(t):
                    if self.has_extern(rx["rate"]["ast"]):
                        return "rate_%d(%s, n, p, Tgas)" % (rx["id"] - 1, t)
                    return self.cexpr(rx["rate"]["ast"], t)
                out = ["kr[%d] = %s;  // %s" % (j, e(tv), rx["eq"])]
                if tv == "Tgas":
                    out.append("{ const double Tp_ = Tgas * (1.0 + h), "
                               "Tm_ = Tgas * (1.0 - h);")
                else:
                    out.append("if (%s == Tgas) { const double Tp_ = %s * "
                               "(1.0 + h), Tm_ = %s * (1.0 - h);"
                               % (tv, tv, tv))
                out.append("  tdk[%d] = ((%s) - (%s)) / (2.0 * h); }"
                           % (j, e("Tp_"), e("Tm_")))
                return out
            w("  for (int j = 0; j < %d; ++j) kr[j] = tdk[j] = 0.0;"
              % len(rxs))
            self.emit_grouped(w, body, select=lambda rx: rx["id"] in idx)
            w("}")
        w("// recombination cooling [erg cm^-3 s^-1] from the network's "
          "rates")
        w("// (thermal.cooling = ['ATOMIC']): n_e n_X kB T (3/2 k + T dk/dT)")
        w("inline double cooling_recomb_(double Tgas, const double* n, "
          "const Params& p) {")
        w("  (void)Tgas; (void)n; (void)p;")
        if not rxs:
            w("  return 0.0;%s" % ("" if self.thermal.get("cool_atomic")
                                   else "  // no thermal.cooling = ['ATOMIC']"))
            w("}")
            w("")
            return
        nr = len(rxs)
        w("  const double T = std::fmin(std::fmax(Tgas, p.Tcmb), 1e9);")
        w("  double k0[%d], tdk[%d];" % (nr, nr))
        w("  recomb_rates_(T, n, p, k0, tdk);")
        w("  double c = 0.0;")
        for j, rx in enumerate(rxs):
            ion = [r["sp"] for r in rx["reactants"] if r["sp"] != "e-"][0]
            w("  c += n[IDX_E] * %s * std::fmax(0.0, 1.5 * k0[%d] + tdk[%d]);"
              "  // %s" % (self.nden(ion), j, j, rx["eq"]))
        w("  return constants::boltzmann_erg * T * c;")
        w("}")
        w("")

    def emit_dh_thermal(self, w):
        """cooling_dH()/heating_dH(): legacy -cooling=DH / -heating=DH, the
        reaction enthalpies at the gas temperature: dH(T) = sum H_reac -
        sum H_prod, H = kB T (H/RT) of Burcat's NASA polynomials
        (the lower block below Tcommon, the upper above, T clamped to the
        data range; electrons 0). Each reaction contributes
        k * prod(n_reactants) * |dH(T)| to the cooling (dH < 0) or to the
        heating (dH > 0). Legacy took every species' enthalpy at 300 K
        from the first-read (upper-range, T > 1000 K) block. Photo-kind
        rates are excluded as legacy (photoheating is the photo module's);
        reactions with a species without data are skipped with a
        warning (legacy silently used H = 0)."""
        rxs, species = [], set()
        if self.dh:
            thermo = self.load_thermo()
            for rx in self.net.reactions:
                if rx["rate"]["kind"] != "expr":
                    continue
                sps = [q["sp"] for q in rx["reactants"] + rx["products"]]
                miss = [x for x in sps
                        if x != "e-" and x not in thermo]
                if miss:
                    self.warnings.append("dH: no thermochemical data for %s in "
                                         "'%s', skipped"
                                         % (", ".join(sorted(set(miss))),
                                            rx["eq"]))
                    continue
                rxs.append(rx)
                species.update(x for x in sps if x != "e-")
            if species:
                self.burcat_used = True
                for ln in BURCAT_COMMENT:
                    w(ln)
            for sp in sorted(species):
                c = thermo[sp]
                w("// H/k [K] of %s (NASA, Burcat; T clamped to the data)"
                  % sp)
                w("static inline double hK_%s(double T) {" % self.cid(sp))
                w("  T = std::fmin(std::fmax(T, %s), %s);"
                  % (self.cnum(c[0]), self.cnum(c[2])))
                for tag, a in (("if (T <= %s) {" % self.cnum(c[1]),
                                c[10:17]), ("} else {", c[3:10])):
                    w("  " + tag)
                    w("    const double a[6] = {%s};"
                      % ", ".join(self.cnum(x) for x in a[:6]))
                    w("    return T * (a[0] + a[1]*T/2.0 + a[2]*T*T/3.0 "
                      "+ a[3]*T*T*T/4.0 + a[4]*T*T*T*T/5.0) + a[5];")
                w("  }")
                w("}")
        for nm, sign in (("cooling_dH", -1.0), ("heating_dH", 1.0)):
            w("inline double %s(double Tgas, const double* n, "
              "const double* k, const Params& p) {" % nm)
            w("  (void)Tgas; (void)n; (void)k; (void)p;")
            if not rxs:
                w("  return 0.0;%s" % ("" if self.dh
                                       else "  // --dh not enabled"))
                w("}")
                continue
            w("  const double kB = constants::boltzmann_erg;")
            w("  double s = 0.0, dh;")
            for rx in rxs:
                h = " + ".join(["hK_%s(Tgas)" % self.cid(q["sp"])
                                for q in rx["reactants"] if q["sp"] != "e-"]
                               or ["0.0"])
                hp = "".join([" - hK_%s(Tgas)" % self.cid(q["sp"])
                              for q in rx["products"] if q["sp"] != "e-"])
                prod = "*".join([self.nden(r["sp"]) for r in
                                 rx["reactants"] if r["kinetic"]])
                w("  dh = %skB * (%s%s);  // %s" % (
                    "-" if sign < 0 else "", h, hp, rx["eq"]))
                w("  if (dh > 0.0) s += k[%d] * %s * dh;" % (rx["id"] - 1,
                                                           prod))
            w("  return s;")
            w("}")
        w("")

    # -- dust photoelectric heating (legacy -heating=PHOTODUST/...) -----------

    def load_solar(self):
        """data/asplund.dat (Asplund+2009 Table 1, fixed-width): element ->
        log10(n/nH) = epsilon - 12 (photospheric)."""
        for base in (self.data_dir.parent, self.data_dir,
                     self.data_dir.parent.parent):
            path = base / "asplund.dat"
            if path.exists():
                break
        else:
            raise FileNotFoundError("asplund.dat not found")
        out = {}
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            # Z  A  Name  Solar  [ErrSolar  Meteor  ErrMeteor]
            m = re.match(r"\s*\d+\s+[\d.]+\s+([A-Za-z]+)\s+(-?[\d.]+)",
                         line)
            if m:
                out[m.group(1)] = float(m.group(2)) - 12.0
        return out

    def load_h2pd(self):
        """data/H2pdB.dat: X levels i (ground energy EX [eV]) x B levels
        j, the transition energy dE [eV] and the prefactor pre (pi e^2 /
        (m_e c) f eps / dE: k = sum pre J(dE) z_i, legacy kpd_H2)"""
        if getattr(self, "_h2pd", None) is None:
            f = self.data_dir.parent / "H2pdB.dat"
            rows = []
            for line in f.read_text().splitlines():
                t = line.split("#")[0].split()
                if len(t) == 5:
                    rows.append((int(t[0]), int(t[1]), float(t[2]),
                                 float(t[3]), float(t[4])))
            nx = max(r[0] for r in rows) + 1
            EX = [0.0] * nx
            for ix, _, ex, _, _ in rows:
                EX[ix] = ex
            self._h2pd = (nx, EX, [r for r in rows if r[4] != 0.0])
        return self._h2pd

    def h2pd_levels(self):
        """(lowest, highest) line energy of each X level, None when the
        level has no line"""
        nx, _, lines = self.load_h2pd()
        out = []
        for j in range(nx):
            e = [r[3] for r in lines if r[0] == j]
            out.append((min(e), max(e)) if e else (None, None))
        return out

    def emit_h2pd(self, w, photo):
        """the photo channels given per bin (h2pd): photo_binned[r] (the X
        level, -1 for the others), photo_bin_xsec() and the level
        populations h2pd_pop()"""
        n = len(photo)
        w("// channels given per photon bin (h2pd lines): the X level, or -1")
        w("inline const int photo_binned[%d] = {%s};"
          % (max(n, 1), ", ".join(str(p.get("binned", -1)) for p in photo)
             or "-1"))
        binned = [p for p in photo if "binned" in p]
        if not binned:
            w("constexpr bool PHOTO_BINNED = false;")
            w("inline void photo_bin_xsec(int, int, const double*, "
              "const double*, double*) {}")
            return
        nx, EX, lines = self.load_h2pd()
        w("constexpr bool PHOTO_BINNED = true;")
        w("// H2pdB.dat (legacy kpd_H2): X-level energies [eV] and the lines")
        w("constexpr int H2PD_NX = %d, H2PD_NL = %d;" % (nx, len(lines)))
        w("inline const double h2pd_EX[%d] = {%s};"
          % (nx, ", ".join(self.cnum(x) for x in EX)))
        self.emit_data_array(w, "h2pd_dE", [r[3] for r in lines],
                             "line energies [eV]")
        self.emit_data_array(w, "h2pd_pre", [r[4] for r in lines],
                             "line prefactors pi e^2 f eps / (m_e c dE)")
        w("inline const int h2pd_X[%d] = {%s};"
          % (len(lines), ", ".join(str(r[0]) for r in lines)))
        w("// bin cross sections of the lines of X level j: legacy k = pre")
        w("// J(dE), with J constant across the bin (legacy photoBinJ), so")
        w("// sigma_b = sum_lines pre h / (4 pi ln(Ehi/Elo)) and k = sum_b")
        w("// sigma_b Phi_b. A line on a shared edge goes to the upper bin")
        w("// (legacy: the last match); a line outside every bin absorbs")
        w("// nothing (legacy stopped)")
        w("inline void photo_bin_xsec(int r, int B, const double* Elo, "
          "const double* Ehi, double* sig) {")
        w("  const int j = photo_binned[r];")
        w("  for (int b = 0; b < B; ++b) sig[b] = 0.0;")
        w("  if (j < 0) return;")
        w("  for (int l = 0; l < H2PD_NL; ++l) {")
        w("    if (h2pd_X[l] != j) continue;")
        w("    int bl = -1;")
        w("    for (int b = 0; b < B; ++b)")
        w("      if (Elo[b] <= h2pd_dE[l] && Ehi[b] >= h2pd_dE[l]) bl = b;")
        w("    if (bl >= 0)")
        w("      sig[bl] += h2pd_pre[l] * constants::planck_eV / (4.0 * "
          "constants::pi * std::log(Ehi[bl] / Elo[bl]));")
        w("  }")
        w("}")
        w("// normalized Boltzmann population of X level j at Tgas (legacy")
        w("// partitionH2_vib; energies relative to the lowest level, so it")
        w("// does not underflow to 0/0 at low T as legacy did)")
        w("KROME_HD inline double h2pd_pop(double jf, double Tgas, "
          "const double* n, const Params& p) {")
        w("  (void)n; (void)p;")
        w("  const double b = 1.0 / (constants::boltzmann_eV * Tgas);")
        w("  double z = 0.0;")
        w("  for (int i = 0; i < H2PD_NX; ++i) z += std::exp(-(h2pd_EX[i] - "
          "h2pd_EX[0]) * b);")
        w("  return std::exp(-(h2pd_EX[static_cast<int>(jf)] - h2pd_EX[0]) "
          "* b) / z;")
        w("}")

    def is_solomon(self, rx):
        """H2 -> H + H by Lyman-Werner line pumping (the Solomon process):
        a J_LW rate, a rate through H2_solomonLW, or a photo() rate whose
        window lies inside the LW band. These get the Solomon heating
        (dissociation 0.4 eV + UV pumping, legacy CHEM / PHOTOAV) and no
        photoelectron deposit. The direct photodissociation continuum
        (H2_sigmaLW, above 14 eV) is not one: legacy CHEM counted every
        photo H2 dissociation, and H2_sigmaLW rates, as Solomon pumping."""
        R = sorted(r["sp"] for r in rx["reactants"])
        P = sorted(p0["sp"] for p0 in rx["products"])
        if len(R) != 1 or R[0] not in self.H2_FORMS or P != ["H", "H"]:
            return False
        r = rx["rate"]
        if r["kind"] == "photo" and r["xsec"]["kind"] == "h2pd":
            return True  # line absorption followed by dissociation
        if r["kind"] == "photo":
            return (r.get("Emax") is not None
                    and r["Emax"] <= knet.LW_BAND[1] + 1e-9)
        ast = r.get("ast")
        return bool(rx.get("lw")) or (ast is not None and any(
            nd.get("f") == "H2_solomonLW" for nd in self.deep_walk(ast)))

    def is_uv_h2diss(self, rx):
        """H2 -> H + H by an analytic interstellar-UV rate: not a photo()
        rate, with an expression that reads the UV field (the parameters
        Av or G0, e.g. 5.6e-11*exp(-3.74*Av)*H2self) and no cosmic-ray
        parameter. Legacy PHOTOAV took the first H2 -> H + H reaction; its
        CHEM term counted only the photo() ones, so these join the PHOTOAV
        rates only (not CHEM's Solomon set)."""
        R = sorted(r["sp"] for r in rx["reactants"])
        P = sorted(p0["sp"] for p0 in rx["products"])
        r = rx["rate"]
        if len(R) != 1 or R[0] not in self.H2_FORMS or P != ["H", "H"] \
                or r["kind"] == "photo":
            return False
        ast = r.get("ast")
        if ast is None:
            return False
        names = {nd.get("s") for nd in self.deep_walk(ast)
                 if nd.get("kind") == "param"}
        return bool(names & {"Av", "G0"}) and not names & {"crate",
                                                            "crflux"}

    def emit_wd01(self, w, have_e):
        """thermal.photoelectric = 'WD01': PHOTODUST = the Weingartner & Draine 2001
        photoelectric heating, eq. (44), and PHOTODUSTNET = eq. (44) minus
        their grain recombination cooling, eq. (45), for the thermal.wd01 row of
        their Tables 2-3, scaled by d2g / 0.01 (the fits are for Milky Way
        dust). psi = G sqrt(T) / n_e with G the Habing flux of the solver
        (the fits hold for 10 <= T <= 1e4 K and 1e2 <= psi <= 1e6; eq. 45
        for T >= 1e3 K: psi is clamped at 1e2 there, where its D1 / ln psi
        diverges at psi = 1). PHOTOAV stays the BT94-era formula."""
        (c0, c1, c2, c3, c4, c5, c6), (d0, d1, d2, d3, d4) = \
            WD01_TABLE[self.wd01]
        cn = self.cnum
        tag = "R_V %g, b_C %g, case %s, %s" % self.wd01
        if not have_e:
            self.warnings.append("thermal.photoelectric = 'WD01' stubbed: no electrons")
        w("// Weingartner & Draine 2001 photoelectric heating, eq. (44),")
        w("// Table 2 row %s [erg cm^-3 s^-1], scaled by d2g / 0.01" % tag)
        w("double inline heating_photoDust(double Tgas, const double* n, "
          "const Params& p, double Ghab) {")
        w("  (void)Tgas; (void)n; (void)p; (void)Ghab;")
        if have_e:
            w("  if (n[IDX_E] <= 0.0 || Ghab <= 0.0) return 0.0;")
            w("  const double psi = Ghab * std::sqrt(Tgas) / n[IDX_E];")
            w("  return 1e-26 * Ghab * Hnuclei_(n) * (%s + %s * std::pow("
              "Tgas, %s))" % (cn(c0), cn(c1), cn(c4)))
            w("      / (1.0 + %s * std::pow(psi, %s) * (1.0 + %s * "
              "std::pow(psi, %s)))" % (cn(c2), cn(c5), cn(c3), cn(c6)))
            w("      * (p.d2g / 1e-2);")
        else:
            w("  return 0.0;")
        w("}")
        w("// net: eq. (44) minus the grain recombination cooling, eq. (45),")
        w("// Table 3 row %s" % tag)
        w("double inline heating_netPhotoDust(double Tgas, const double* n, "
          "const Params& p, double Ghab) {")
        w("  (void)Tgas; (void)n; (void)p; (void)Ghab;")
        if have_e:
            w("  if (n[IDX_E] <= 0.0) return 0.0;")
            w("  const double psi = Ghab * std::sqrt(Tgas) / n[IDX_E];")
            w("  const double x = std::log(std::fmax(psi, 1e2));")
            w("  const double cool = 1e-28 * n[IDX_E] * Hnuclei_(n) * std::pow("
              "Tgas, %s + %s / x)" % (cn(d0), cn(d1)))
            w("      * std::exp(%s + %s * x - %s * x * x) * (p.d2g / 1e-2);"
              % (cn(d2), cn(d3), cn(d4)))
            w("  return heating_photoDust(Tgas, n, p, Ghab) - cool;")
        else:
            w("  return 0.0;")
        w("}")
        self.emit_photoAv(w)

    def emit_dust_heating(self, w):
        """heating_photoDust / heating_netPhotoDust / heating_photoAv:
        legacy -heating=PHOTODUST / PHOTODUSTNET / PHOTOAV, ported from
        krome_heating.f90. The Habing flux is an explicit argument computed
        by the thermal solver (p.G0 * exp(-2.5 Av), or the photobin field's
        Habing integral when radiation is present) — legacy hardcoded 1.69
        (Draine78) in PHOTODUST unless the user supplied @var G0/Av; KROME-X
        makes the field explicit. PHOTODUSTNET scales by D/D_sol =
        p.d2g/0.01 (legacy dust2gas_ratio)."""
        net = self.net
        th = self.thermal
        have_e = "e-" in net.species
        # metallicity reference: first of Fe, C, O, Si whose neutral atom
        # is a network species (legacy #KROME_photoDustZ), scaled to the
        # Asplund+09 solar log10(n/nH); the sum runs over every species
        # containing that atom, weighted by multiplicity (legacy zGets)
        zel = next((zz for zz in ("Fe", "C", "O", "Si")
                    if zz in net.species), None)
        eps_lines = [
            "  const double eps = 4.9e-2 / (1.0 + 4e-3 * "
            "std::pow(psi, 0.73))",
            "      + 3.7e-2 * std::pow(Tgas * 1e-4, 0.7) / "
            "(1.0 + 2e-4 * psi);  // BT94 eq. 43 / W03"]
        if zel and th["heat_photoDust"]:
            solar = self.load_solar()[zel]
            terms = []
            for s in sorted(net.species):
                mult = (net.species[s]["composition"] or {}).get(zel, 0)
                if mult:
                    terms.append(("%d.0*" % mult if mult > 1 else "")
                                 + "n[IDX_%s]" % self.cid(s))
            w("// metallicity wrt Asplund+09 solar, reference element %s"
              " (legacy get_metallicity%s)" % (zel, zel))
            w("static inline double zscale_photoDust_(const double* n) {")
            w("  const double nH = Hnuclei_(n);")
            w("  const double zx = std::fmax(%s, 0.0);"
              % " + ".join(terms))
            w("  return std::pow(10.0, std::log10(zx / nH + 1e-40) - (%s));"
              % self.cnum(solar))
            w("}")
        if self.photoelectric == "WD01" and (th["heat_photoDust"]
                                             or th["heat_netPhotoDust"]):
            return self.emit_wd01(w, have_e)
        if zel and have_e and th["heat_photoDust"]:
            w("// Bakes & Tielens 1994 photoelectric heating (Wolfire+03")
            w("// absorbed-energy tweak), local-flux form without grain")
            w("// recombination (legacy heat_photoDust); Ghab in Habing")
            w("double inline heating_photoDust(double Tgas, "
              "const double* n, const Params& p, double Ghab) {")
            w("  (void)p;")
            w("  const double nH = Hnuclei_(n);")
            w("  double psi = 1e99;")
            w("  if (n[IDX_E] > 0.0) psi = Ghab * std::sqrt(Tgas) "
              "/ n[IDX_E];")
            for ln in eps_lines:
                w(ln)
            w("  return 1.3e-24 * eps * Ghab * nH * "
              "zscale_photoDust_(n);")
            w("}")
        else:
            if th["heat_photoDust"]:
                self.warnings.append("heating_photoDust stubbed: %s"
                                     % ("no Fe/C/O/Si species for "
                                        "metallicity" if not zel
                                        else "no electrons"))
            w("double inline heating_photoDust(double, const double*, "
              "const Params&, double) { return 0.0; }")
        if have_e and "H" in net.species and th["heat_netPhotoDust"]:
            w("// net photoelectric heating: BT94 eqs. 42+44 including the")
            w("// grain recombination cooling (legacy heat_netPhotoDust),")
            w("// scaled by D/D_sol = d2g/0.01")
            w("double inline heating_netPhotoDust(double Tgas, "
              "const double* n, const Params& p, double Ghab) {")
            w("  const double nH = Hnuclei_(n);")
            w("  const double bet = 0.735 * std::pow(Tgas, -0.068);")
            w("  double psi = 0.0;")
            w("  if (n[IDX_E] > 0.0) psi = Ghab * std::sqrt(Tgas) "
              "/ n[IDX_E];")
            w("  const double recomb_cool = 4.65e-30 * "
              "std::pow(Tgas, 0.94)")
            w("      * std::pow(psi, bet) * n[IDX_E] * n[IDX_H];")
            for ln in eps_lines:
                w(ln)
            w("  return (1.3e-24 * eps * Ghab * nH - recomb_cool) "
              "* (p.d2g / 1e-2);")
            w("}")
        else:
            if th["heat_netPhotoDust"]:
                self.warnings.append("heating_netPhotoDust stubbed: needs "
                                     "e- and H")
            w("double inline heating_netPhotoDust(double, const double*, "
              "const Params&, double) { return 0.0; }")
        self.emit_photoAv(w)

    def emit_photoAv(self, w):
        """heating_photoAv (legacy -heating=PHOTOAV)"""
        net, th = self.net, self.thermal
        # PHOTOAV: the network's Solomon-process rates, and its analytic
        # interstellar-UV H2 dissociation rates (legacy took the first
        # H2 -> H + H reaction, whatever its process)
        solc = [rx for rx in net.reactions if self.is_solomon(rx)]
        if solc and th["heat_photoAv"] and th["heat_chem"]:
            raise ValueError("thermal.heating = ['CHEM'] and PHOTOAV both add the "
                             "Solomon-process heating: choose one")
        sol = [rx for rx in net.reactions
               if self.is_solomon(rx) or self.is_uv_h2diss(rx)]
        if sol and "H" in net.species and self.h2_forms() \
                and th["heat_photoAv"]:
            w("// H2 photodissociation (0.4 eV) + UV pumping (2.7e-11 erg")
            w("// * thermalized fraction) heating, legacy heat_photoAv;")
            w("// Rdiss = the Solomon / UV dissociation rates: %s"
              % ", ".join("k[%d]" % (rx["id"] - 1) for rx in sol))
            w("double inline heating_photoAv(double Tgas, const double* n, "
              "const double* k, const Params& p) {")
            w("  (void)p;")
            w("  const double dd = Hnuclei_(n);")
            w("  const double ncr = 1.0e6 / std::sqrt(Tgas)")
            w("      / (1.6 * std::exp(-std::pow(4.0e2 / Tgas, 2.0)) "
              "* n[IDX_H] / dd")
            w("         + 1.4 * std::exp(-1.2e4 / (Tgas + 1.2e3)) "
              "* nH2_(n) / dd);")
            w("  const double h2heatfac = 1.0 / (1.0 + ncr / dd);")
            if all(rx["reactants"][0]["sp"] == "H2" for rx in sol):
                w("  const double Rdiss = %s;"
                  % " + ".join("k[%d]" % (rx["id"] - 1) for rx in sol))
                w("  return 6.4e-13 * Rdiss * n[IDX_H2]")
                w("      + 2.7e-11 * Rdiss * h2heatfac * n[IDX_H2];")
            else:
                # nuclear-spin H2 species: each rate times its own reactant
                w("  const double Rn = %s;" % " + ".join(
                    "k[%d] * n[IDX_%s]" % (rx["id"] - 1,
                                           self.cid(rx["reactants"][0]["sp"]))
                    for rx in sol))
                w("  return 6.4e-13 * Rn + 2.7e-11 * Rn * h2heatfac;")
            w("}")
        else:
            if th["heat_photoAv"]:
                self.warnings.append("heating_photoAv stubbed: needs a "
                                     "Solomon-process or UV-field (Av, G0) "
                                     "H2 -> H + H reaction")
            w("double inline heating_photoAv(double, const double*, "
              "const double*, const Params&) { return 0.0; }")
        w("")

    def emit_tables(self, w):
        if self.net.tables:
            w("static inline double tab_interp1(double x, const double* tx, "
              "const double* tv, int n) {")
            w("  if (x <= tx[0]) return tv[0];")
            w("  if (x >= tx[n-1]) return tv[n-1];")
            w("  int lo = 0, hi = n - 1;")
            w("  while (hi - lo > 1) { int m = (lo+hi)/2; "
              "if (tx[m] <= x) lo = m; else hi = m; }")
            w("  const double f = (x - tx[lo]) / (tx[lo+1] - tx[lo]);")
            w("  return tv[lo] + f * (tv[lo+1] - tv[lo]);")
            w("}")
            w("static inline double tab_interp2(double x, double y, "
              "const double* tx, const double* ty, const double* tv, "
              "int nx, int ny) {")
            w("  // bilinear on a rectilinear grid, clamped at the edges")
            w("  if (x < tx[0]) x = tx[0];")
            w("  if (x > tx[nx-1]) x = tx[nx-1];")
            w("  if (y < ty[0]) y = ty[0];")
            w("  if (y > ty[ny-1]) y = ty[ny-1];")
            w("  int i = 0, j = 0, hi;")
            w("  for (hi = nx - 1; hi - i > 1;) { int m = (i+hi)/2; "
              "if (tx[m] <= x) i = m; else hi = m; }")
            w("  for (hi = ny - 1; hi - j > 1;) { int m = (j+hi)/2; "
              "if (ty[m] <= y) j = m; else hi = m; }")
            w("  const double fx = (x - tx[i]) / (tx[i+1] - tx[i]);")
            w("  const double fy = (y - ty[j]) / (ty[j+1] - ty[j]);")
            w("  const double v00 = tv[i*ny + j], v01 = tv[i*ny + j + 1];")
            w("  const double v10 = tv[(i+1)*ny + j], "
              "v11 = tv[(i+1)*ny + j + 1];")
            w("  return v00*(1-fx)*(1-fy) + v10*fx*(1-fy) + v01*(1-fx)*fy "
              "+ v11*fx*fy;")
            w("}")
        for name, fname in sorted(self.net.tables.items()):
            try:
                kind, data = self.load_table(fname)
            except (FileNotFoundError, ValueError) as exc:
                self.warnings.append("table %s: %s (stubbed to 0)"
                                     % (name, exc))
                w("template <class... A> static inline double %s(A...) "
                  "{ return 0.0; }  // STUB" % name)
                continue
            if name in getattr(self.net, "table_log", ()):
                # `table ... log`: interpolated in log10 of the axes, and
                # called with log10 of the axis quantities
                def lg(axis):
                    if min(axis) <= 0.0:
                        raise ValueError(
                            "table %s (%s) declared log has an axis value "
                            "<= 0" % (name, fname))
                    return [math.log10(v) for v in axis]
                if kind == "1d":
                    data = (lg(data[0]), data[1])
                else:
                    data = (lg(data[0]), lg(data[1]), data[2])
                w("// table %s: interpolated in log10 of its axes (called "
                  "with log10 values)" % name)
            if kind == "1d":
                X, V = data
                w("// table %s: %s (%d points)" % (name, fname, len(X)))
                # device-safe (emit_data_array -> devicify accessors), as
                # every table a KROME_HD function reads
                self.emit_data_array(w, "tab_%s_x" % name, X,
                                     "table %s x" % name)
                self.emit_data_array(w, "tab_%s_v" % name, V,
                                     "table %s values" % name)
                w("static inline double %s(double x) {" % name)
                w("  return tab_interp1(x, tab_%s_x, tab_%s_v, %d);"
                  % (name, name, len(X)))
                w("}")
            else:
                xs, ys, grid = data
                w("// table %s: %s (%dx%d grid)"
                  % (name, fname, len(xs), len(ys)))
                self.emit_data_array(w, "tab_%s_x" % name, xs,
                                     "table %s x" % name)
                self.emit_data_array(w, "tab_%s_y" % name, ys,
                                     "table %s y" % name)
                self.emit_data_array(w, "tab_%s_v" % name,
                                     [v for row in grid for v in row],
                                     "table %s values" % name)
                w("static inline double %s(double x, double y) {" % name)
                w("  return tab_interp2(x, y, tab_%s_x, tab_%s_y, "
                  "tab_%s_v, %d, %d);" % (name, name, name, len(xs), len(ys)))
                w("}")
        w("")

    def resolve_photo(self, rx):
        """Resolve one photo() rate to C++ xsec code + energy window."""
        r = rx["rate"]
        i = rx["id"] - 1
        out = {"i": i, "rx": rx, "code": None,
               "reactant": rx["reactants"][0]["sp"],
               "Emin": r.get("Emin", 0.0), "Emax": r.get("Emax")}
        x = r["xsec"]
        if x["kind"] == "expr":
            out["code"] = self.cexpr(x["ast"])
        elif x["kind"] == "auto":
            reacts = sorted(rr["sp"] for rr in rx["reactants"])
            prods = sorted(p["sp"] for p in rx["products"])
            segs = [e for e in self.photo_db
                    if e.get("type") == "photoion"
                    and e.get("reacts") == reacts and e.get("prods") == prods]
            segs.sort(key=lambda e: e["limits"][0])
            if segs:
                # piecewise-in-E: each segment applies from its own threshold
                # up to the next segment's threshold (the outer-shell fit,
                # then the sum over the shells above the inner-shell edge);
                # the LAST segment is extrapolated above its upper limit
                code = None
                for e in segs:
                    expr = re.sub(r"\benergy_ev\b", "E_eV", d2e(e["rate"]),
                                  flags=re.I)
                    code = expr if code is None else \
                        "(E_eV >= %s ? %s : %s)" \
                        % (self.cnum(e["limits"][0]), expr, code)
                out["code"] = code
                if "Emin" not in r:
                    out["Emin"] = segs[0]["limits"][0]
                # no Emax from the database: extrapolate (photo() args win)
                db_top = segs[-1]["limits"][1]
                if db_top < 4.9e4:
                    self.warnings.append(
                        "photo auto %s: database covers only up to %g eV "
                        "(inner shells missing); outermost-shell fit "
                        "extrapolated above" % (rx["eq"], db_top))
            else:
                self.warnings.append("photo auto: no database entry for %s"
                                     % rx["eq"])
        else:
            ref = (x.get("ref") or "").lower()
            tab = None
            try:
                if ref == "leiden":
                    tab = self.load_leiden(rx)
                elif ref == "swri":
                    tab = self.load_swri(rx)
            except FileNotFoundError as exc:
                self.warnings.append("photo %s: %s" % (rx["eq"], exc))
            if tab:
                out["table"] = tab
                if "Emin" not in r:
                    out["Emin"] = tab[0][0]
                if "Emax" not in r:
                    out["Emax"] = tab[0][-1]
            elif ref not in ("leiden", "swri"):
                self.warnings.append("photo xsec source %r not implemented "
                                     "yet (%s)" % (x.get("ref"), rx["eq"]))
        return out

    def load_leiden(self, rx):
        return leiden_table(self.data_dir, rx)

    def load_swri(self, rx):
        return swri_table(self.data_dir, rx, self.warnings)

    def emit_photo(self, w, photo):
        if photo:
            w("// Verner+96 cross-section fit [cm^2] (as in legacy KROME)")
            w("static inline double sigma_v96(double E, double E0, "
              "double s0, double ya, double P, double yw, double y0, "
              "double y1) {")
            w("  const double x = E/E0 - y0;")
            w("  const double y = std::sqrt(x*x + y1*y1);")
            w("  const double Fy = ((x - 1.0)*(x - 1.0) + yw*yw)")
            w("      * std::pow(y, 0.5*P - 5.5)")
            w("      * std::pow(1.0 + std::sqrt(y/ya), -P);")
            w("  return 1e-18 * s0 * Fy;")
            w("}")
            w("")
            w("// Verner & Yakovlev 1995 fit of one shell [cm^2], 0 below its")
            w("// threshold Eth; above the inner-shell edge the database sums")
            w("// them over the shells (Verner's phfit2.f)")
            w("static inline double sigma_vy95(double E, double Eth, "
              "double E0, double s0, double ya, double P, double yw, "
              "double l) {")
            w("  if (E < Eth) return 0.0;")
            w("  const double y = E/E0;")
            w("  return 1e-18 * s0 * ((y - 1.0)*(y - 1.0) + yw*yw)")
            w("      * std::pow(y, 0.5*P - 5.5 - l)")
            w("      * std::pow(1.0 + std::sqrt(y/ya), -P);")
            w("}")
            w("")
        if any("table" in ph for ph in photo):
            w("// linear interpolation on a tabulated cross section;")
            w("// zero outside the tabulated range")
            w("static inline double xsec_interp(double E, const double* tE, "
              "const double* tS, int n) {")
            w("  if (E < tE[0] || E > tE[n-1]) return 0.0;")
            w("  int lo = 0, hi = n - 1;")
            w("  while (hi - lo > 1) {")
            w("    int m = (lo + hi) / 2;")
            w("    if (tE[m] <= E) lo = m; else hi = m;")
            w("  }")
            w("  const double f = (E - tE[lo]) / (tE[lo+1] - tE[lo]);")
            w("  return tS[lo] + f * (tS[lo+1] - tS[lo]);")
            w("}")
            w("// exact integral of that interpolant over [a, b] [cm^2 eV]")
            w("// (the Doppler broadening of PhotoModule::init)")
            w("static inline double xsec_int(double a, double b, "
              "const double* tE, const double* tS, int n) {")
            w("  a = std::fmax(a, tE[0]);")
            w("  b = std::fmin(b, tE[n-1]);")
            w("  if (!(b > a)) return 0.0;")
            w("  int i = 0, hi = n - 1;")
            w("  while (hi - i > 1) {")
            w("    int m = (i + hi) / 2;")
            w("    if (tE[m] <= a) i = m; else hi = m;")
            w("  }")
            w("  double s = 0.0;")
            w("  for (; i < n - 1 && tE[i] < b; ++i) {")
            w("    const double l = std::fmax(a, tE[i]);")
            w("    const double r = std::fmin(b, tE[i+1]);")
            w("    if (!(r > l)) continue;")
            w("    const double g = (tS[i+1] - tS[i]) / (tE[i+1] - tE[i]);")
            w("    s += 0.5 * (2.0 * tS[i] + g * (l + r - 2.0 * tE[i])) "
              "* (r - l);")
            w("  }")
            w("  return s;")
            w("}")
            w("")
        for ph in photo:
            if "table" in ph:
                E, S = ph["table"]
                w("// tabulated xsec of %s: %d points, %.4g-%.4g eV"
                  % (ph["rx"]["eq"], len(E), E[0], E[-1]))
                # device-safe: nvcc compiles every KROME_HD body, so a
                # plain host array here fails even if no kernel calls it
                self.emit_data_array(w, "xsec_tab_%d_E" % ph["ch"], E,
                                     "xsec %d energies [eV]" % ph["ch"])
                self.emit_data_array(w, "xsec_tab_%d_S" % ph["ch"], S,
                                     "xsec %d [cm^2]" % ph["ch"])
                w("static inline double photo_xsec_%d(double E_eV) {"
                  % ph["ch"])
                w("  return xsec_interp(E_eV, xsec_tab_%d_E, xsec_tab_%d_S, "
                  "%d);" % (ph["ch"], ph["ch"], len(E)))
                w("}")
                w("static inline double photo_xsec_int_%d(double a, double b)"
                  " {" % ph["ch"])
                w("  return xsec_int(a, b, xsec_tab_%d_E, xsec_tab_%d_S, %d);"
                  % (ph["ch"], ph["ch"], len(E)))
                w("}")
                continue
            if "binned" in ph:
                w("// %s, H2pdB.dat lines of X level v=%d: given per bin"
                  % (ph["rx"]["eq"], ph["binned"]))
                w("static inline double photo_xsec_%d(double E_eV) {"
                  % ph["ch"])
                w("  (void)E_eV;")
                w("  return 0.0;  // lines: see photo_bin_xsec")
                w("}")
                continue
            w("// xsec of %s [cm^2]" % ph["rx"]["eq"])
            w("static inline double photo_xsec_%d(double E_eV) {" % ph["ch"])
            w("  (void)E_eV;")
            w("  return %s;" % (ph["code"] if ph["code"] else
                                "0.0;  // STUB (unresolved xsec source)"))
            w("}")
        n = len(photo)
        w("")
        w("constexpr int NPHOTO = %d;" % n)
        w("typedef double (*XsecFn)(double);")
        w("inline const int photo_reaction[%d] = {%s};"
          % (max(n, 1), ", ".join(str(p["i"]) for p in photo) or "-1"))
        w("inline const int photo_reactant[%d] = {%s};"
          % (max(n, 1), ", ".join("IDX_" + self.cid(p["reactant"])
                                  for p in photo) or "-1"))
        w("inline const double photo_Emin[%d] = {%s};"
          % (max(n, 1), ", ".join(self.cnum(p["Emin"]) for p in photo)
             or "0.0"))
        w("inline const double photo_Emax[%d] = {%s};"
          % (max(n, 1),
             ", ".join(self.cnum(p["Emax"] if p["Emax"] is not None
                                 else math.inf) for p in photo) or "0.0"))
        w("inline const XsecFn photo_xsec[%d] = {%s};"
          % (max(n, 1), ", ".join("photo_xsec_%d" % p["ch"] for p in photo)
             or "nullptr"))
        w("// exact integrals of the tabulated (file) cross sections, the")
        w("// ones legacy -useBroadening broadened; nullptr for the others")
        w("typedef double (*XsecIntFn)(double, double);")
        w("inline const XsecIntFn photo_xsec_int[%d] = {%s};"
          % (max(n, 1), ", ".join(("photo_xsec_int_%d" % p["ch"])
                                  if "table" in p else "nullptr"
                                  for p in photo) or "nullptr"))
        self.emit_h2pd(w, photo)
        w("")

    # shielding externs whose first argument is the absorber column: a
    # factor containing exactly one of them gets a column form (the
    # curve-of-growth optical depth of get_opacity)
    SHIELD_COL_EXTERNS = ("fselfH2", "fselfH2_WG11", "fselfH2_R14")

    def col_form(self, ast):
        """copy of ast with the column argument of its single shielding
        extern call replaced by Ncol, or None"""
        calls = [nd for nd in knet.walk(ast)
                 if nd.get("f") in self.SHIELD_COL_EXTERNS]
        if len(calls) != 1 or not calls[0].get("a"):
            return None
        calls[0]["_mark"] = True
        c = copy.deepcopy(ast)
        del calls[0]["_mark"]
        for nd in knet.walk(c):
            if nd.pop("_mark", False):
                nd["a"][0] = {"s": "Ncol", "kind": "ncol"}
        return c

    def emit_rad_channels(self, w, photo):
        """network-format.md 5.4: shielded photo reactions
        (photo(shield=...)), the photons absorbed per event (; pdiss=), and
        the J_LW channel. A J_LW rate k = khat(T, n) J_LW acts on the photon
        bins as a photo reaction with cross section khat w_b per bin (w_b the
        band weight, J_LW = sum_b w_b Phi_b); every channel removes
        rate/pdiss photons."""
        net = self.net
        n = len(photo)
        m1 = max(n, 1)
        w("// ---- radiative channels (network-format.md 5.4)")
        w("// photons absorbed per photo event: 1/pdiss")
        w("inline const double photo_pdiss[%d] = {%s};"
          % (m1, ", ".join(self.cnum(p["rx"].get("pdiss", 1.0))
                           for p in photo) or "1.0"))
        w("// photoelectron deposit (E - Eth) of each photo rate: none for")
        w("// the Solomon process (line pumping; its heating is the")
        w("// Solomon term of CHEM / PHOTOAV)")
        w("inline const bool photo_deposit[%d] = {%s};"
          % (m1, ", ".join("false" if self.is_solomon(p["rx"]) else "true"
                           for p in photo) or "true"))
        shielded = [p for p in photo if p.get("shield")]
        w("constexpr bool PHOTO_SHIELDED = %s;"
          % ("true" if shielded else "false"))
        w("inline const bool photo_shielded[%d] = {%s};"
          % (m1, ", ".join("true" if p.get("shield")
                           else "false" for p in photo) or "false"))
        w("// shielding factor s_r of each photo rate (1: unshielded)")
        w("inline void photo_shields(double Tgas, const double* n, "
          "const Params& p, double* s) {")
        w("  (void)Tgas; (void)n; (void)p; (void)s;")
        for q, ph in enumerate(photo):
            sh = ph.get("shield")
            w("  s[%d] = %s;  // %s" % (q, self.cexpr(sh) if sh else "1.0",
                                         ph["rx"]["eq"]))
        w("}")
        colf = []
        w("// s_r with its absorber column set to Ncol (column form)")
        w("inline double photo_shield_col(int r, double Ncol, double Tgas, "
          "const double* n, const Params& p) {")
        w("  (void)Ncol; (void)Tgas; (void)n; (void)p;")
        w("  switch (r) {")
        for q, ph in enumerate(photo):
            sh = ph.get("shield")
            c = self.col_form(sh) if sh else None
            colf.append(c is not None)
            if c is not None:
                w("    case %d: return %s;" % (q, self.cexpr(c)))
            elif sh and "binned" not in ph:  # a level population: no column
                self.warnings.append(
                    "shield of %s has no column form (one fselfH2* call "
                    "with the column as first argument): get_opacity uses "
                    "s(N) N" % ph["rx"]["eq"])
        w("    default: return 1.0;")
        w("  }")
        w("}")
        w("inline const bool photo_shield_colform[%d] = {%s};"
          % (m1, ", ".join("true" if c else "false" for c in colf)
             or "false"))
        lw = [rx for rx in net.reactions if rx.get("lw")]
        m2 = max(len(lw), 1)
        w("// J_LW rates: k = khat(T, n) J_LW")
        w("constexpr int NLW = %d;" % len(lw))
        w("inline const int lw_reaction[%d] = {%s};"
          % (m2, ", ".join(str(rx["id"] - 1) for rx in lw) or "-1"))
        w("inline const int lw_reactant[%d] = {%s};"
          % (m2, ", ".join("IDX_" + self.cid(next(
              q["sp"] for q in rx["reactants"] if q["kinetic"]))
              for rx in lw) or "-1"))
        w("inline const double lw_pdiss[%d] = {%s};"
          % (m2, ", ".join(self.cnum(rx["pdiss"]) for rx in lw) or "1.0"))
        self._jlw_unit = True
        w("inline void lw_unit_rates(double Tgas, const double* n, "
          "const Params& p, double* kh) {")
        w("  (void)Tgas; (void)n; (void)p; (void)kh;")
        for q, rx in enumerate(lw):
            w("  kh[%d] = %s;  // %s" % (q, self.cexpr(rx["rate"]["ast"]),
                                         rx["eq"]))
        w("}")
        colf = []
        w("inline double lw_unit_rate_col(int q, double Ncol, double Tgas, "
          "const double* n, const Params& p) {")
        w("  (void)Ncol; (void)Tgas; (void)n; (void)p;")
        w("  switch (q) {")
        for q, rx in enumerate(lw):
            c = self.col_form(rx["rate"]["ast"])
            colf.append(c is not None)
            if c is not None:
                w("    case %d: return %s;" % (q, self.cexpr(c)))
            else:
                self.warnings.append(
                    "J_LW rate %s has no column form: get_opacity uses "
                    "khat(N) N" % rx["eq"])
        w("    default: return 0.0;")
        w("  }")
        w("}")
        self._jlw_unit = False
        w("inline const bool lw_colform[%d] = {%s};"
          % (m2, ", ".join("true" if c else "false" for c in colf)
             or "false"))
        w("")

    RESERVED_LOCALS = {"T", "Tgas", "n", "p", "k", "d", "np_", "jj", "h_",
                       "dk_", "q", "flux", "J", "s"}

    def hoist_vars(self, w, rxs, tvar, indent):
        """Hoist the vars directly referenced by these rates into locals of
        the current scope (evaluated once at this scope's temperature)."""
        names = sorted({nd["s"] for rx in rxs
                        for nd in knet.walk(rx["rate"]["ast"])
                        if nd.get("kind") == "var"})
        self._hoisted = {s: (s if s not in self.RESERVED_LOCALS
                             else "vh_" + s) for s in names}
        for s in names:
            w("%sconst double %s = v_%s(%s, n, p);"
              % (indent, self._hoisted[s], s, tvar))
            # unused when every rate that needs it is an out-of-line
            # rate_<i>() function
            w("%s(void)%s;" % (indent, self._hoisted[s]))

    def emit_secion(self, w):
        """Per-species secondary ionization by fast photoelectrons
        (radiation.secondary_ionization = ['H', 'H2', 'He', 'C']=3.92,...): emits the channel
        tables consumed by the photo module's Quadri et al. 2025 energy
        split. Fit species (H, H2, He) use the paper's phi_i fits; metals
        take an explicit sigma_A/sigma_H ratio and rescale the H rate."""
        if not self.secion:
            return
        rows = []
        table = self.load_secion_ratios()
        for tok in self.secion:
            sp, _, ratio = tok.partition("=")
            sp = sp.strip()
            ion = sp + "+"
            for need in (sp, ion, "e-"):
                if need not in self.net.species:
                    raise ValueError("secondary ionization of %s needs "
                                     "species %r in the network" % (sp, need))
            if sp in self.Q25_FIT:
                if ratio:
                    raise ValueError("secondary ionization: %s is a fit "
                                     "species, no ratio expected" % sp)
                kind, eth = self.Q25_FIT[sp]
                rows.append((sp, ion, kind, eth, 1.0))
            else:
                # metal: ratio from the data table (an explicit X=ratio in
                # the option overrides it)
                r = float(ratio) if ratio else table.get(sp)
                if r is None:
                    raise ValueError(
                        "secondary ionization of %s: no entry in %s and no "
                        "explicit ratio given (use %s=<sigma_ei ratio>)"
                        % (sp, self.SECION_RATIOS_FILE.name, sp))
                rows.append((sp, ion, self.Q25_KIND_METAL, 0.0, r))
        if "H" not in self.net.species:
            raise ValueError("secondary ionization needs species H in the "
                             "network (mixture ratio n_H2/n_H and metal "
                             "scaling are defined per H atom)")
        w("// secondary ionization by fast photoelectrons, Quadri et al. "
          "2025 App. B")
        w("// (radiation.secondary_ionization=%s)" % ",".join(self.secion))
        w("#define KROME_HAS_SECION 1")
        w("constexpr int N_SECION = %d;" % len(rows))
        w("constexpr int SECION_KIND_METAL = %d;" % self.Q25_KIND_METAL)
        w("inline const int secion_target[N_SECION] = {%s};"
          % ", ".join("IDX_" + self.cid(r[0]) for r in rows))
        w("inline const int secion_ion[N_SECION] = {%s};"
          % ", ".join("IDX_" + self.cid(r[1]) for r in rows))
        w("inline const int secion_kind[N_SECION] = {%s};"
          % ", ".join(str(r[2]) for r in rows))
        w("inline const double secion_Eth_eV[N_SECION] = {%s};"
          % ", ".join(self.cnum(r[3]) for r in rows))
        w("inline const double secion_ratio[N_SECION] = {%s};"
          % ", ".join(self.cnum(r[4]) for r in rows))
        w("")

    @staticmethod
    def lit_num(nd):
        """the value of a numeric literal node, also a negated one (the
        parser writes x**(-0.5) as unary minus), else None"""
        if is_num(nd):
            return float(nd["n"])
        if nd.get("u") == "-" and is_num(nd.get("x", {})):
            return -float(nd["x"]["n"])
        return None

    @staticmethod
    def cheap_power(e):
        """exponents printed without std::pow (ipow_ / sqrt forms)"""
        return e is not None and abs(e) <= 32 and (
            e.is_integer() or ((2.0 * e).is_integer() and abs(e) <= 2.5))

    def rate_is_fixed(self, rx):
        """a rate of Tgas and the parameters only: no densities (ntot,
        Hnuclei, n(X)), externs, tables, J_LW or Tdust (the live
        table_Tdust / binned-dust hook). Constant during a fixed-Tgas
        solve (rates_fixed). A piecewise group is fixed only as a whole
        (its range conditions are emitted relative to its pieces)."""
        def one(r):
            ast = r["rate"]["ast"]
            if self.has_extern(ast) or self.ndep(ast):
                return False
            return not any(nd.get("s") == "Tdust"
                           or (nd.get("kind") == "extern"
                               and not self.state_free_extern(nd, True))
                           for nd in self.deep_walk(ast))
        g = rx.get("group")
        if g is None:
            return one(rx)
        cache = self.__dict__.setdefault("_group_fixed", {})
        if g not in cache:
            cache[g] = all(one(r) for r in self.net.reactions
                           if r.get("group") == g
                           and r["rate"]["kind"] == "expr")
        return cache[g]

    def emit_grouped(self, w, body, hoist=False, select=None):
        """Emit statements from body(rx, tvar) for every expr-rate reaction,
        wrapped in the piecewise group selection/clamp structure. rates()
        and the jacobian dk terms share this, so they always agree on which
        piece is active and at which (clamped) temperature it is evaluated.
        With hoist=True the vars used in each scope become local constants.
        select(rx) restricts the emission to a subset of the reactions."""
        singles, groups = [], {}
        for rx in self.net.reactions:
            if rx["rate"]["kind"] != "expr":
                continue
            if select is not None and not select(rx):
                continue
            if rx.get("group") is None:
                singles.append(rx)
            else:
                groups.setdefault(rx["group"], []).append(rx)
        if hoist and singles:
            self.hoist_vars(w, singles, "Tgas", "  ")
        for rx in singles:
            for line in body(rx, "Tgas"):
                w("  " + line)
        self._hoisted = {}
        for g, pieces in sorted(groups.items()):
            pieces.sort(key=lambda r: r["trange"])
            bodies = [body(rx, "T") for rx in pieces]
            if not any(bodies):
                continue
            glo, ghi = pieces[0]["trange"][0], pieces[-1]["trange"][1]
            w("  { // group %d: %s (%d pieces, outside=%s)"
              % (g, pieces[0]["eq"], len(pieces), self.outside))
            w("    double T = Tgas;")
            clamp_lo = self.outside == "clamp"
            zero_lo = self.outside in ("zero", "asym") and glo > 0.0
            clamp_hi = self.outside in ("clamp", "asym")
            if clamp_lo and glo > 0.0:
                w("    if (T < %s) T = %s;" % (self.cnum(glo),
                                               self.cnum(glo)))
            if clamp_hi and ghi != math.inf:
                w("    if (T > %s) T = %s;" % (self.cnum(ghi),
                                               self.cnum(ghi)))
            if hoist:
                self.hoist_vars(w, pieces, "T", "    ")
                bodies = [body(rx, "T") for rx in pieces]  # re-render hoisted
            cond0 = "T >= %s && " % self.cnum(glo) if zero_lo else ""
            for j, rx in enumerate(pieces):
                lo, hi = rx["trange"]
                last = j == len(pieces) - 1
                cmp_hi = "<=" if (last and hi != math.inf) else "<"
                conds = []
                if j > 0 or zero_lo:
                    conds.append("T >= %s" % self.cnum(lo))
                if hi != math.inf and not (last and clamp_hi):
                    conds.append("T %s %s" % (cmp_hi, self.cnum(hi)))
                cond = " && ".join(conds) if conds else "true"
                kw = "if" if j == 0 else "else if"
                w("    %s (%s%s) {" % (kw, cond0 if j == 0 else "", cond))
                for line in bodies[j]:
                    w("      " + line)
                w("    }")
                cond0 = ""
            w("  }")
            self._hoisted = {}

    def kin_counts(self, rx):
        targets = {}
        for r in rx["reactants"]:
            if r["kinetic"]:
                targets[r["sp"]] = targets.get(r["sp"], 0) + 1
        return targets

    def jac_updates(self, rx, tgt, slot=False):
        """row updates of the derivative 'd' of one reaction's flux with
        respect to the density of species tgt (slot=False), or to the
        state slot tgt directly (slot=True, used by the symbolic dk part
        whose diff() is already slot-based). In gas+total mode a surface
        member's density maps to two slots (+total, -gas) and the row set
        includes the pairs' total slots."""
        if not self.gas_total:
            col = "IDX_%s" % self.cid(tgt)
            lines = []
            for r in rx["reactants"]:
                lines.append("J[IDX_%s*jrs_ + %s*jcs_] -= d;"
                             % (self.cid(r["sp"]), col))
            for pr in rx["products"]:
                lines.append("J[IDX_%s*jrs_ + %s*jcs_] += d;"
                             % (self.cid(pr["sp"]), col))
            return lines
        lines = []
        for colsp, csign in ([(tgt, 1)] if slot else self.den_cols(tgt)):
            col = "IDX_%s" % self.cid(colsp)
            for sp, v in self.stoich_net(rx):
                v *= csign
                op = "+=" if v > 0 else "-="
                mag = "" if abs(v) == 1 else "%d.0*" % abs(v)
                lines.append("J[IDX_%s*jrs_ + %s*jcs_] %s %sd;"
                             % (self.cid(sp), col, op, mag))
        return lines

    RHS_LOOP_MIN = 500  # --rhs=auto: loops above this many reactions

    def rhs_loop(self):
        """the mass-action rhs / fluxes / Jacobian as loops over index
        tables (rates.rhs = 'loop', or auto above RHS_LOOP_MIN reactions)"""
        nrea = len(self.net.reactions)
        if nrea == 0:
            return False
        return self.rhs_mode == "loop" or (
            self.rhs_mode == "auto" and nrea > self.RHS_LOOP_MIN)

    def loop_tables(self):
        """the index tables of the loop rhs and Jacobian, built from the
        same helpers as the unrolled code (so that both perform the same
        operations in the same order). Per reaction r:
          flux factors  kf_sp/kf_pr[kf_off[r]..]: the kinetic reactants in
                        reaction order; a density is n[sp], or n[sp] -
                        n[pr] for a dust.gas_total surface slot (pr >= 0)
          scatter rows  sc_row/sc_c[sc_off[r]..]: dn[row] += c * flux
                        (reactants -1, products +1, duplicates kept; the
                        net increments with dust.gas_total)
          Jacobian      jt_*[jt_off[r]..]: one term per distinct kinetic
                        reactant (sorted by name): multiplicity jt_m, the
                        factors jf_sp/jf_pr[jt_foff[t]..] (sorted by name),
                        the columns jc_col/jc_sgn[jt_coff[t]..]; each column
                        takes the rows of the scatter with c * sign"""
        sp_idx = {x: i for i, x in enumerate(sorted(self.net.species))}

        def pr_of(sp):
            g = self.pair_gas.get(sp)
            return sp_idx[g] if g is not None else -1
        T = {k: [] for k in ("kf_off", "kf_sp", "kf_pr", "sc_off", "sc_row",
                             "sc_c", "jt_off", "jt_m", "jt_foff", "jf_sp",
                             "jf_pr", "jt_coff", "jc_col", "jc_sgn")}
        for rx in self.net.reactions:
            T["kf_off"].append(len(T["kf_sp"]))
            for r in rx["reactants"]:
                if r["kinetic"]:
                    T["kf_sp"].append(sp_idx[r["sp"]])
                    T["kf_pr"].append(pr_of(r["sp"]))
            T["sc_off"].append(len(T["sc_row"]))
            if self.gas_total:
                rows = [(sp_idx[x], v) for x, v in self.stoich_net(rx)]
            else:
                rows = ([(sp_idx[r["sp"]], -1) for r in rx["reactants"]]
                        + [(sp_idx[q["sp"]], 1) for q in rx["products"]])
            for row, c in rows:
                T["sc_row"].append(row)
                T["sc_c"].append(c)
            T["jt_off"].append(len(T["jt_m"]))
            targets = self.kin_counts(rx)
            for tgt, mult in sorted(targets.items()):
                T["jt_m"].append(mult)
                T["jt_foff"].append(len(T["jf_sp"]))
                for x, c in sorted(targets.items()):
                    for _ in range(c - 1 if x == tgt else c):
                        T["jf_sp"].append(sp_idx[x])
                        T["jf_pr"].append(pr_of(x))
                T["jt_coff"].append(len(T["jc_col"]))
                cols = self.den_cols(tgt) if self.gas_total else [(tgt, 1)]
                for colsp, sgn in cols:
                    T["jc_col"].append(sp_idx[colsp])
                    T["jc_sgn"].append(sgn)
        for k, tot in (("kf_off", "kf_sp"), ("sc_off", "sc_row"),
                       ("jt_off", "jt_m"), ("jt_foff", "jf_sp"),
                       ("jt_coff", "jc_col")):
            T[k].append(len(T[tot]))
        return T

    def emit_loop_tables(self, w):
        """the tables of loop_tables() as int arrays (device-safe: the
        generated headers' constant arrays become accessors on device)"""
        T = self.loop_tables()
        w("// rates.rhs = 'loop': the mass-action rhs, fluxes and Jacobian as loops")
        w("// over these index tables (Codegen.loop_tables)")
        for name, vals in T.items():
            vals = vals or [0]  # (no zero-length arrays)
            w("inline const int lt_%s[%d] = {" % (name, len(vals)))
            for q in range(0, len(vals), 16):
                w("  %s%s" % (", ".join(str(v) for v in vals[q:q + 16]),
                              "," if q + 16 < len(vals) else ""))
            w("};")
        w("// a density of the loop tables: n[sp], or n[sp] - n[pr] (gas+total)")
        w("#define KROME_LT_DEN(sp, pr) ((pr) < 0 ? n[sp] : (n[sp] - n[pr]))")
        w("// x += c v, written as the unrolled code does (+= v, -= v, += c.0*v)")
        w("#define KROME_LT_ADD(x, c, v) \\")
        w("  do { const int c_ = (c); \\")
        w("    if (c_ == 1) (x) += (v); else if (c_ == -1) (x) -= (v); \\")
        w("    else if (c_ > 0) (x) += static_cast<double>(c_) * (v); \\")
        w("    else (x) -= static_cast<double>(-c_) * (v); } while (0)")
        w("")

    JAC_EL_RE = re.compile(r"\bJ\[IDX_(\w+)\*jrs_ \+ (IDX_\w+|jj)\*jcs_\]")
    JAC_U_RE = re.compile(r"\bU\[IDX_(\w+)\*NAGG \+ \d+\]")

    def jac_pattern(self):
        """the species block of the analytic Jacobian as a set of (row,
        column) indices: every element the generated jacobian_split()
        writes (the mass-action and the rate-derivative parts, read from
        the same emitters), the rows with aggregate terms (U W^T) dense,
        and the diagonal (the Newton matrix I - gamma J). The entries the
        runtime adds by difference quotients (photo shielding, secondary
        ionization, dust exchange) are not in it (solve_cell.hpp
        sparse_exact_())."""
        names = sorted(self.net.species)
        idx = {self.cid(x): i for i, x in enumerate(names)}
        ns = len(names)
        pat = {(i, i) for i in range(ns)}
        lines = []
        for rx in self.net.reactions:
            lines += self.jac_stoich_lines(rx)
            if rx["rate"]["kind"] == "expr":
                lines += self.jac_dk_lines(rx, "Tgas")
        for ln in lines:
            for r, c in self.JAC_EL_RE.findall(ln):
                if c == "jj":
                    pat.update((idx[r], j) for j in range(ns))
                else:
                    pat.add((idx[r], idx[c[4:]]))
            for r in self.JAC_U_RE.findall(ln):
                pat.update((idx[r], j) for j in range(ns))
        return ns, pat

    @staticmethod
    def min_degree_order(ns, pat):
        """a minimum-degree elimination order of the symmetrized pattern
        (A + A^T; ties: the lowest index), on the explicit elimination
        graph (bit sets): perm[k] = the species eliminated k-th"""
        adj = [0] * ns
        for i, j in pat:
            if i != j:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
        left = set(range(ns))
        deg = [bin(a).count("1") for a in adj]
        perm = []
        while left:
            v = min(left, key=lambda u: (deg[u], u))
            left.discard(v)
            perm.append(v)
            nb = adj[v]
            bit = 1 << v
            u = nb
            while u:
                low = u & -u
                w = low.bit_length() - 1
                u ^= low
                adj[w] = (adj[w] | nb) & ~(bit | low)
                deg[w] = bin(adj[w]).count("1")
        return perm

    def sparse_lu_tables(self):
        """the static-pivot sparse LU of the species block of the Newton
        matrix (bdf_cell.hpp, Opts::sparse): the elimination order (the
        legacy DLSODES way: a minimum-degree ordering, no pivoting) and the
        rows of L + U with their fill, in pivot order:
          perm[k]   the species of pivot k
          rp[k]     row k at ci[rp[k] .. rp[k+1]), columns (pivot order)
                    ascending; L left of the diagonal, U from it on
          di[k]     the position of the diagonal in row k"""
        ns, pat = self.jac_pattern()
        perm = self.min_degree_order(ns, pat)
        inv = [0] * ns
        for k, s in enumerate(perm):
            inv[s] = k
        rows = [0] * ns
        for i, j in pat:
            rows[inv[i]] |= 1 << inv[j]
        for k in range(ns):  # symbolic elimination (no pivoting)
            up = rows[k] >> (k + 1) << (k + 1)
            colk = 1 << k
            for i in range(k + 1, ns):
                if rows[i] & colk:
                    rows[i] |= up
        rp, ci, di = [0], [], []
        for k in range(ns):
            r = rows[k]
            cols = []
            while r:
                low = r & -r
                cols.append(low.bit_length() - 1)
                r ^= low
            di.append(len(ci) + cols.index(k))
            ci += cols
            rp.append(len(ci))
        flops = 0
        for k in range(ns):
            nl = sum(1 for i in range(k + 1, ns) if rows[i] >> k & 1)
            flops += nl * (1 + 2 * (rp[k + 1] - di[k] - 1))
        return dict(perm=perm, rp=rp, ci=ci, di=di, nnz_a=len(pat),
                    flops=flops)

    def emit_sparse_tables(self, w):
        """the tables of sparse_lu_tables() (device-safe int arrays),
        SPARSE_NNZ and SPARSE_FLOPS (one factorization's)"""
        T = self.sparse_lu_tables()
        ns = len(T["perm"])
        w("// the static-pivot sparse LU of the species block of the Newton")
        w("// matrix (Codegen.sparse_lu_tables; LINSOL_SPARSE): %d nonzeros"
          % T["nnz_a"])
        w("// in I - gamma J, %d in L + U (%.1f%% of NSPEC^2), %d flops per"
          % (len(T["ci"]), 100.0 * len(T["ci"]) / max(ns * ns, 1),
             T["flops"]))
        w("// factorization (dense: %d)" % (2 * ns ** 3 // 3))
        w("constexpr int SPARSE_NNZ = %d;" % len(T["ci"]))
        w("constexpr double SPARSE_FLOPS = %d;" % T["flops"])
        for name in ("perm", "rp", "ci", "di"):
            vals = T[name] or [0]
            w("inline const int sp_%s[%d] = {" % (name, len(vals)))
            for q in range(0, len(vals), 16):
                w("  %s%s" % (", ".join(str(v) for v in vals[q:q + 16]),
                              "," if q + 16 < len(vals) else ""))
            w("};")
        w("")

    def jac_stoich_lines(self, rx):
        i = rx["id"] - 1
        targets = self.kin_counts(rx)
        lines = []
        for tgt, mult in sorted(targets.items()):
            rest = ["k[%d]" % i]
            for s, c in sorted(targets.items()):
                rest += [self.nden(s)] * (c - 1 if s == tgt else c)
            coef = ("%d.0*" % mult) if mult > 1 else ""
            lines.append("d = %s;  // %s wrt %s" % (coef + "*".join(rest),
                                                    rx["eq"], tgt))
            lines += self.jac_updates(rx, tgt)
        return lines

    def state_free_extern(self, nd, fixed=False):
        """an extern call whose value depends on the densities only
        through its explicit arguments, and these are n-independent
        (EXTERN_STATE); fixed=True also requires it to be constant during
        a fixed-Tgas solve ("field" only without radiation.evolve)."""
        if nd.get("kind") != "extern":
            return False
        cls = EXTERN_STATE.get(nd["f"])
        if cls is None or any(self._arg_ndep(x) for x in nd["a"]):
            return False
        return not (fixed and cls == "field" and self.evolve_rad)

    def _arg_ndep(self, a):
        """n-dependence of an extern argument: densities, ntot, Hnuclei,
        the live Tdust, n-dependent vars or externs (a bare species index,
        as in n_global(X), is a label, not a density)"""
        for nd in knet.walk(a):
            k = nd.get("kind")
            if k == "density" or nd.get("s") in ("ntot", "Hnuclei"):
                return True
            if nd.get("s") == "Tdust" and self.dust_tables:
                return True
            if k == "var" and self.ndep(self.net.vars[nd["s"]]):
                return True
            if k in ("extern", "table") and nd is not a \
                    and not self.state_free_extern(nd):
                return True
        return False

    def has_extern(self, ast):
        """the rate needs the numerical dk/dn path (and a callable
        rate_<i>): an extern/table call that may depend on the densities,
        the live Tdust of the dust tables, or J_LW"""
        return any((nd.get("kind") in ("extern", "table")
                    and not self.state_free_extern(nd))
                   or (nd.get("s") == "Tdust" and self.dust_tables)
                   or nd.get("s") == "J_LW"
                   for nd in self.deep_walk(ast))

    def flux_rows(self, rx, lhs):
        """row updates `lhs % row += / -= d` of a reaction's flux
        derivative d (net slot increments in gas+total mode)"""
        if self.gas_total:
            out = []
            for sp, v in self.stoich_net(rx):
                op = "+=" if v > 0 else "-="
                mag = "" if abs(v) == 1 else "%d.0*" % abs(v)
                out.append("%s %s %sd;" % (lhs % ("IDX_" + self.cid(sp)),
                                           op, mag))
            return out
        return (["%s -= d;" % (lhs % ("IDX_" + self.cid(r["sp"])))
                 for r in rx["reactants"]]
                + ["%s += d;" % (lhs % ("IDX_" + self.cid(pr["sp"])))
                   for pr in rx["products"]])

    AGG_INDEX = {"ntot": 0, "Hnuclei": 1, "rho": 2}

    def extern_deps(self, nd):
        """(direct species, aggregates) through which one extern/table
        call depends on the densities beyond its explicit arguments, or
        None when unknown (runtime-module externs reading n themselves,
        the live Tdust of the dust tables, dust.gas_total slots): the rate
        then takes the full numerical dk/dn loop over every species."""
        if nd.get("kind") == "table" or self.state_free_extern(nd):
            return set(), set()
        f = nd["f"]
        if not f.startswith("surf_") or self.gas_total or self.dust_tables:
            return None
        # generated grain-surface functions (emit_surface): the sites
        # follow rho (surf_sites_), the coverage the water ice
        # (surf_fice_); Tdust is the parameter or the binned-dust hook
        if f in ("surf_chem_PC", "surf_chem_CP"):
            return set(), set()
        # (a single grain population's sites are fixed: no rho)
        rho = set() if self.net.surface_grains is not None else {"rho"}
        if f in ("surf_adsorption", "surf_adsorption_ns", "surf_chem_PCCP",
                 "surf_chem_CC"):
            return set(), rho
        if f == "surf_desorption" or f.startswith("surf_2body_"):
            return set(getattr(self, "_surf_ice", ())), rho
        return None

    def rate_deps(self, ast):
        """(direct species, aggregates) of a rate's dependence on the
        densities, or None when an extern's dependence is unknown"""
        direct, aggs = set(), set()
        for nd in self.deep_walk(ast):
            if nd.get("kind") == "density":
                sp = nd["a"][0]["sp"]
                direct.add(sp)
                if sp in self.pair_gas:
                    direct.add(self.pair_gas[sp])
            elif nd.get("s") in ("ntot", "Hnuclei"):
                aggs.add(nd["s"])
            elif nd.get("s") == "Tdust" and self.dust_tables:
                return None
            elif nd.get("kind") in ("extern", "table"):
                dd = self.extern_deps(nd)
                if dd is None:
                    return None
                direct |= dd[0]
                aggs |= dd[1]
        return direct, aggs

    def jac_dk_lines(self, rx, tvar):
        """the rate-derivative terms of one reaction: dk/dn_j prod(n) into
        J (the direct dependence) and dk/d(aggregate) prod(n) into U (the
        coupling through ntot, Hnuclei and rho, J += U W^T)"""
        i = rx["id"] - 1
        ast = rx["rate"]["ast"]
        targets = self.kin_counts(rx)
        prod = "*".join(self.nden(s) for s, c in sorted(targets.items())
                        for _ in range(c)) or "1.0"
        # extern/table calls cannot be differentiated symbolically: their
        # dk/dn enters the Jacobian by finite differences of the rate
        # function itself (evaluated at the same clamped temperature)
        if self.has_extern(ast):
            deps = self.rate_deps(ast)
            if deps is None:
                return self.jac_dk_full_loop(rx, tvar, prod)
            direct, aggs = deps
            direct = sorted(direct)
            mask = sum(1 << self.AGG_INDEX[a] for a in aggs)
            lines = [
                "{ // numerical dk/dn: %s (direct: %s; aggregates: %s)"
                % (rx["eq"], ", ".join(direct) or "-",
                   ", ".join(sorted(aggs)) or "-"),
                "  double np_[NSPEC];",
                "  for (int q = 0; q < NSPEC; ++q) np_[q] = n[q];",
                "  auto rf_ = [&](const double* x) { return rate_%d(%s, x, p, "
                "Tgas); };" % (i, tvar),
                "  const double pr_ = %s;" % prod,
                "  double ca_[NAGG] = {0.0};",
            ]
            if mask:
                lines.append("  const int dir_[] = {%s};"
                             % ", ".join(["IDX_" + self.cid(x)
                                          for x in direct] or ["-1"]))
                lines.append("  fd_dk_agg_(rf_, np_, n, k[%d], %d, dir_, %d,"
                             " ca_);" % (i, mask, len(direct)))
            for sp in direct:
                col = "IDX_" + self.cid(sp)
                lines.append("  d = (fd_dk_(rf_, np_, n, k[%d], %s) - "
                             "agg_part_(ca_, %s)) * pr_;  // wrt %s"
                             % (i, col, col, sp))
                lines += ["  " + x for x in self.flux_rows(
                    rx, "J[%%s*jrs_ + %s*jcs_]" % col)]
            for a in sorted(aggs, key=self.AGG_INDEX.get):
                lines.append("  d = ca_[%d] * pr_;  // wrt %s"
                             % (self.AGG_INDEX[a], a))
                lines += ["  " + x for x in self.flux_rows(
                    rx, "U[%%s*NAGG + %d]" % self.AGG_INDEX[a])]
            lines.append("}")
            return lines
        if not self.ndep(ast):
            return []
        # targets are STATE SLOTS (diff() is slot-based): in gas+total
        # mode a surface density adds its gas slot
        dk_targets = set()
        aggs = set()
        for nd in self.deep_walk(ast):
            if nd.get("kind") == "density":
                sp = nd["a"][0]["sp"]
                dk_targets.add(sp)
                if sp in self.pair_gas:
                    dk_targets.add(self.pair_gas[sp])
            elif nd.get("s") in ("ntot", "Hnuclei"):
                aggs.add(nd["s"])
        lines = []
        for tgt in sorted(dk_targets):
            dk = self.diff(ast, tgt)
            if is_num(dk, 0.0):
                continue
            lines.append("d = (%s)*%s;  // dk of %s wrt %s"
                         % (self.cexpr(dk, tvar), prod, rx["eq"], tgt))
            lines += self.jac_updates(rx, tgt, slot=True)
        for a in sorted(aggs, key=self.AGG_INDEX.get):
            dk = self.diff(ast, "@" + a)
            if is_num(dk, 0.0):
                continue
            lines.append("d = (%s)*%s;  // dk of %s wrt %s"
                         % (self.cexpr(dk, tvar), prod, rx["eq"], a))
            lines += self.flux_rows(rx, "U[%%s*NAGG + %d]"
                                    % self.AGG_INDEX[a])
        return lines

    def jac_dk_full_loop(self, rx, tvar, prod):
        """numerical dk/dn of a rate whose density dependence is unknown:
        difference quotients along every species (the total derivative,
        aggregates included, straight into J)"""
        i = rx["id"] - 1
        lines = [
            "{ // numerical dk/dn: %s (extern/table in rate, all species)"
            % rx["eq"],
            "  double np_[NSPEC];",
            "  for (int q = 0; q < NSPEC; ++q) np_[q] = n[q];",
            "  auto rf_ = [&](const double* x) { return rate_%d(%s, x, p, "
            "Tgas); };" % (i, tvar),
            "  const double pr_ = %s;" % prod,
            "  for (int jj = 0; jj < NSPEC; ++jj) {",
            "    const double dk_ = fd_dk_(rf_, np_, n, k[%d], jj);" % i,
            "    if (dk_ == 0.0) continue;",
            "    d = dk_ * pr_;",
        ]
        lines += ["    " + x for x in self.flux_rows(rx, "J[%s*jrs_ + jj*jcs_]")]
        lines += ["  }", "}"]
        return lines

    def deep_walk(self, ast):
        """walk() through var references too."""
        for nd in knet.walk(ast):
            yield nd
            if nd.get("kind") == "var":
                yield from self.deep_walk(self.net.vars[nd["s"]])


# ---- numbers in data files --------------------------------------------------

def _fortran_real(v):
    """a Fortran real as written in legacy data files (1d-40, 1.D-6), as a
    Python float literal"""
    if v is None:
        raise TypeError
    t = v.strip().replace("d", "e").replace("D", "e")
    try:
        float(t)
    except ValueError:
        raise ValueError("%r is not a number" % v)
    return t


TESTS_DIR = paths.ROOT / "examples" / "tests"


def list_tests():
    return sorted(d.name for d in TESTS_DIR.iterdir()
                  if (d / "options.toml").is_file()) \
        if TESTS_DIR.is_dir() else []


def is_legacy_network(path):
    """a legacy react_* file: its first statement is not the KROME-X header"""
    try:
        for ln in Path(path).read_text().splitlines():
            t = ln.strip()
            if t and not t.startswith("#"):
                return not t.startswith("krome_format")
    except OSError:
        return False
    return False


def convert_legacy_network(path, outdir):
    """convert a legacy network into the project (outdir/<name>.knet, the
    file the project is then generated from and recorded with)"""
    import convert_legacy
    path = Path(path)
    name = path.name[len("react_"):] if path.name.startswith("react_") \
        else path.stem
    out = Path(outdir) / (name + ".knet")
    c = convert_legacy.Converter()
    c.convert(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(c.render(path.name))
    print("NOTE: legacy network %s converted to %s (kromex/convert_legacy.py;"
          " %d reactions)" % (path, out, c.n_reactions), file=sys.stderr)
    for w in c.warnings:
        print("WARNING: %s: %s" % (path.name, w), file=sys.stderr)
    return out


def build_parser():
    """the command-line parser: every generation option, its help and
    default (also the source of the options template)"""
    ap = argparse.ArgumentParser(
        description="The KROME-X code generator, driven by kromex."
                    "Preprocessor (and python3 -m kromex): these arguments are its "
                    "internal interface, docs/config-design.md")
    ap.add_argument("--version", action="version",
                    version="KROME-X %s" % krome_version(),
                    help="print the KROME-X version and exit")
    ap.add_argument("network_pos", nargs="?", type=Path, metavar="network",
                    help="the .knet network file (or -n / project.network)")
    ap.add_argument("-n", "--network", type=Path, default=None,
                    help="the .knet network file (legacy -n)")
    ap.add_argument("-o", "--outdir", type=Path, default=None,
                    help="the project directory (default generated/, or "
                         "build_NAME with --test)")
    ap.add_argument("--test", metavar="NAME", default=None,
                    help="also write the files of the shipped test NAME "
                         "(examples/tests/NAME: its driver, plot script and "
                         "data; no example), whose options the caller gives "
                         "(kromex.Preprocessor.from_test)")
    ap.add_argument("--outside", choices=("asym", "clamp", "zero"),
                    default="asym",
                    help="rate outside its declared T range: zero below + "
                         "clamp above (asym, default), clamp both, or "
                         "legacy zero")
    ap.add_argument("--cooling", default="",
                    help="thermal cooling processes, FIXED in the "
                         "generated code (legacy -cooling= names): "
                         + ", ".join(sorted(Codegen.COOLING_NAMES))
                         + ". Implies the tables they need (CO/OH/H2O/HCN "
                         "molecular, Z metal-line, Z_CIE/Z_CIENOUV Shen+13, "
                         "Z_EXTENDED both, DH enthalpies)")
    ap.add_argument("--heating", default="",
                    help="thermal heating processes, FIXED in the "
                         "generated code (legacy -heating= names): "
                         + ", ".join(sorted(Codegen.HEATING_NAMES)))
    ap.add_argument("--cooling-floor", nargs="?", const="DEFAULT",
                    default="", metavar="TERMS",
                    help="cooling floor at the per-cell parameter Tfloor: "
                         "each floored cooler is evaluated as cool(Tgas) - "
                         "cool(Tfloor). Without TERMS: every configured "
                         "radiative line cooler; with TERMS (legacy "
                         "-useIndividualFloor): only those, from "
                         + ", ".join(sorted(Codegen.FLOOR_NAMES)))
    ap.add_argument("--shielding", choices=("WG11", "DB96", "R14", "NONE"),
                    default="WG11", type=str.upper,
                    help="H2 self-shielding model behind the legacy "
                         "krome_fshield() (legacy -shielding): Wolcott-"
                         "Green+ 2011 (default), Draine & Bertoldi 1996, "
                         "Richings+ 2014, or NONE (fshield = 1, legacy's "
                         "default)")
    ap.add_argument("--evolve-rad", action="store_true",
                    help="co-evolve the photon field with the chemistry "
                         "(as ln Phi, updated in place); without it the "
                         "field stays fixed at the caller's values")
    ap.add_argument("--evolve-Td", action="store_true",
                    help="per-bin dust temperatures as solver state, "
                         "following each grain's energy equation (heat "
                         "capacity of Draine & Li 2001; legacy -usedTdust "
                         "tracked the balance instead; needs "
                         "thermal.cooling = ['DUSTBINS']); without it they are the "
                         "radiative-collisional balance of the state")
    ap.add_argument("--evolve-dust", default="", metavar="PROCESSES",
                    help="the grain radius of every dust bin as solver "
                         "state (per cell, CellsOpts::adust / the adust "
                         "argument), changed by GROWTH (accretion of the "
                         "key element C / Si), SPUTTER (thermal sputtering, "
                         "Nozawa et al. 2006 yields) and EVAP (sublimation, "
                         "Grassi et al. 2017 sec. 2.3); the atoms go to / "
                         "from the gas-phase key element (legacy "
                         "-dustOptions=GROWTH,SPUTTER,EVAP); needs dust "
                         "bins at runtime")
    ap.add_argument("--secondary-ionization", default="",
                    help="comma-separated species receiving secondary "
                         "ionization by fast photoelectrons (SvS85 split), "
                         "e.g. H,He")
    ap.add_argument("--num2col", choices=("default", "jeans", "jeans40",
                                          "length"),
                    default="default",
                    help="column density estimator N(n) of the molecular "
                         "coolants, the H2 opacity, the line escape "
                         "probability and the shielding (legacy "
                         "-columnDensityMethod): default 1.87e21 "
                         "(n/1e3)^(2/3), jeans / jeans40 half a Jeans "
                         "length, length n Lcol with the per-cell "
                         "parameter Lcol [cm]")
    ap.add_argument("--metal-cooling", default="",
                    help="with thermal.cooling = ['Z'] or Z_EXTENDED: restrict the "
                         "metal-line coolants to these coolZ.dat metals "
                         "(C,Si,Fe,O,S,C+,O+,Si+,Fe+; default: all present)")
    ap.add_argument("--cool-file", default="", metavar="FILE[,FILE]",
                    help="coolant files in the legacy coolZ.dat format "
                         "(LAMDA / CHIANTI conversions, e.g. "
                         "data/cooling/coolO2.dat;"
                         " legacy -coolFile), added to the embedded "
                         "coolZ.dat; select their coolants with thermal.cooling="
                         "<name> (CII, FeII for atoms and ions, the species "
                         "for molecules: O2) or thermal.cooling = ['FILE'] (all), or "
                         "with thermal.metals; each is its own cooling "
                         "term")
    ap.add_argument("--cool-escape", choices=("", "LVG"), default="",
                    help="optical depth of the line coolants (the coolZ.dat "
                         "metals and the thermal.coolant_files coolants): LVG, the "
                         "Sobolev escape probability (1 - exp(-tau)) / tau "
                         "per line, tau from the species column (radiation.num2col)"
                         " and the line width sqrt(2 kB T / m + vturb^2) "
                         "(per-cell parameter vturb [cm s^-1], default 0);"
                         " default: optically thin")
    ap.add_argument("--cool-levels", type=int, default=None, metavar="L",
                    help="keep only levels 0..L of every metal / file "
                         "coolant, with their lines (legacy -coolLevels)")
    ap.add_argument("--dust-tables",
                    choices=("", "HM2012", "MC", "3D", "PHOTO3D"),
                    default="",
                    help="embed the legacy dust tables (Tdust, dust "
                         "cooling, and the rate of krate_H2dust_table()): "
                         "2D for the HM2012 or MC radiation field, 3D with "
                         "the per-cell Av, or PHOTO3D with Av from the "
                         "dust absorption of the photon-bin field")
    ap.add_argument("--gas-total", action="store_true",
                    help="gas+total representation: each gas<->surface "
                         "pair's surface slot stores TOTAL = gas + "
                         "surface (stability option, removes the stiff "
                         "freezeout/evaporation exchange from it)")
    ap.add_argument("--single", action="store_true",
                    help="write one self-contained net_<name>.hpp instead "
                         "of the split core/rates/thermal/jacobian files")
    ap.add_argument("--fortran", action="store_true",
                    help="also emit the optional Fortran binding: module "
                         "krome (runtime/interfaces/fortran/krome.f90), the generated "
                         "constants module krome_net.f90, example.f90 and "
                         "the Makefile targets fortran / example_f90")
    ap.add_argument("--python", action="store_true",
                    help="also emit the optional Python binding: krome.py "
                         "(ctypes + numpy over the C interface), "
                         "example.py and the Makefile target shared "
                         "(libkrome.so/.dylib)")
    ap.add_argument("--gamma", default="CONST", type=str.upper,
                    metavar="{%s}" % ",".join(GAMMA_MODES),
                    help="adiabatic index of the thermal ODE (legacy "
                         "-gamma): CONST = the per-cell parameter gamma_ad "
                         "(default 5/3); EXACT, POPOVAS compute it from "
                         "the composition and Tgas")
    ap.add_argument("--h2-opacity", default="", type=str.upper,
                    choices=("", "RIPAMONTI", "OMUKAI"),
                    help="optically thick suppression of the H2 line "
                         "cooling (legacy -H2opacity): "
                         "Ripamonti & Abel 2004 or the Omukai escape table")
    ap.add_argument("--photoelectric", default="BT94", type=str.upper,
                    choices=("BT94", "WD01"),
                    help="the dust photoelectric heating of PHOTODUST / "
                         "PHOTODUSTNET: Bakes & Tielens 1994 (default, "
                         "legacy) or Weingartner & Draine 2001 eqs. 44-45")
    ap.add_argument("--wd01", default="", metavar="RV,bC,case,field",
                    help="the WD01 Table 2/3 row for thermal.photoelectric = 'WD01' "
                         "(default %s: R_V 3.1, b_C 6e-5, case A, the "
                         "ISRF; field B0 = a 3e4 K blackbody)"
                    % WD01_DEFAULT)
    ap.add_argument("--photo-induced", action="store_true",
                    help="radiative pumping of the metal-line levels by the "
                         "photon field (legacy -usePhotoInduced; stimulated "
                         "emission and absorption)")
    ap.add_argument("--cmb-pumping", action="store_true",
                    help="the CMB (blackbody at the per-cell parameter "
                         "Tcmb) in the metal-line level balance: stimulated "
                         "emission and absorption; the cooling floor is "
                         "dropped for the pumped lines (Z, Z_EXTENDED)")
    ap.add_argument("--habing-dust-shielding", action="store_true",
                    help="attenuate the photon-field Habing flux by dust, "
                         "exp(-2e-21 (d2g/0.01) 0.665 N_H) with num2col "
                         "columns (legacy -shieldHabingDust)")
    ap.add_argument("--no-runtime", action="store_true",
                    help="bare generated files only: do not copy the "
                         "runtime headers or write Makefile/example.cpp "
                         "into the output (project) directory")
    ap.add_argument("--rhs", choices=("auto", "unrolled", "loop"),
                    default="auto",
                    help="the mass-action rhs, fluxes and Jacobian as "
                         "generated statements (unrolled) or as loops over "
                         "index tables (loop: smaller code, much faster to "
                         "compile for large networks; the same operations "
                         "in the same order; legacy -iRHS); auto: loops "
                         "above %d reactions" % Codegen.RHS_LOOP_MIN)
    ap.add_argument("--reactant-multiplicity", action="store_true",
                    help="flux k prod(n)/c! for a species appearing c "
                         "times among the kinetic reactants (legacy "
                         "-nuclearMult; opt-in: use it when the rates are "
                         "per ordered pair)")
    ap.add_argument("--screening", action="store_true",
                    help="plasma screening of the reactions between charged "
                         "nuclei (legacy krome_stars: Itoh / Graboske / "
                         "DeWitt), from the state composition and density")
    ap.add_argument("--binding-energies", choices=("KROME", "SEMENOV"),
                    default="KROME",
                    help="surface binding energies: data/Ebare_ice.dat "
                         "(E_diff = 2/3 E_bind) or the OSU set of Semenov+10 "
                         "(Ebare_ice_Semenov.dat, E_diff = 0.77 E_bind)")
    ap.add_argument("--surface-nu", choices=("HASEGAWA", "LEGACY"),
                    default="HASEGAWA",
                    help="surface vibration frequency: Hasegawa+92 "
                         "sqrt(2 n_s E_b/(pi^2 m)) or legacy 1e12 s^-1")
    ap.add_argument("--no-surface-competition", action="store_true",
                    help="two-body surface rates without reaction-diffusion "
                         "competition (legacy)")
    ap.add_argument("--no-surface-tunnelling", action="store_true",
                    help="two-body surface reactions with a barrier Ea by "
                         "thermal hopping only, exp(-Ea/Tdust) (the legacy "
                         "Semenov framework); default: the larger of that "
                         "and tunnelling through a 1 A rectangular barrier")
    ap.add_argument("--reactive-desorption",
                    choices=("NONE", "FIXED", "MINISSALE"), default="NONE",
                    help="gas branch of the two-body surface reactions: as "
                         "written (NONE), the legacy database fractions "
                         "(FIXED), or Minissale+16 (ice: /10)")
    ap.add_argument("--surface-bins", type=int, default=0, metavar="NB",
                    help="per-bin surface species (legacy @surface): every "
                         "surface species and reaction cloned for NB dust "
                         "bins; the runtime DustBins must have NB bins "
                         "(nd x grain types). Default: one ice abundance "
                         "per species, rates area-weighted over the bins")
    ap.add_argument("--cr-desorption", default="", metavar="PARAM",
                    help="cosmic-ray spot-heating desorption (Hasegawa & "
                         "Herbst 1993) driven by the CR ionization-rate "
                         "parameter PARAM (e.g. crate)")
    ap.add_argument("--rtol", type=float, default=None,
                    help="the project's default relative tolerance "
                         "(SolverOpts::rtol starts from it; default 1e-8; "
                         "legacy -RTOL); still changeable at runtime")
    ap.add_argument("--atol", type=float, default=None,
                    help="the project's default absolute tolerance of the "
                         "species [cm^-3] (SolverOpts::atol_n; default "
                         "1e-20; legacy -ATOL); still changeable at "
                         "runtime")
    ap.add_argument("--custom-atol", default="", metavar="FILE",
                    help="per-species default absolute tolerances: a file "
                         "of 'SPECIES ATOL' lines (legacy -customATOL); "
                         "the other species keep solver.atol. The starting "
                         "SolverOpts::atol_species, changeable at runtime")
    ap.add_argument("--custom-rtol", default="", metavar="FILE",
                    help="per-species default relative tolerances: a file "
                         "of 'SPECIES RTOL' lines (legacy -customRTOL); "
                         "the other species keep solver.rtol. The starting "
                         "SolverOpts::rtol_species, changeable at runtime")
    ap.add_argument("--check-reverse", action="store_true",
                    help="write reverse_check.txt: each reaction's reverse "
                         "in the network (or its absence) and dH(300 K) "
                         "from the thermochemical data (legacy "
                         "-checkReverse)")
    ap.add_argument("--check", action="store_true",
                    help="check the network instead of generating a project "
                         "(with the other options given): the generated "
                         "code against an independent evaluation of the "
                         "network, rates that are not finite, negative or "
                         "implausibly large, species never formed or never "
                         "destroyed, reactions written twice; exit status 1 "
                         "on a failure (kromex/netcheck.py)")
    return ap


# the example value of each option that takes one and has no choices
OPTION_EXAMPLES = {
    "network": "networks/primordial.knet", "outdir": "myproject",
    "cooling": "ATOMIC,H2", "heating": "CHEM,COMPRESS",
    "metal_cooling": "C,C+,O,Si,Fe", "cool_file": "data/cooling/coolO2.dat",
    "cool_levels": "10", "secondary_ionization": "H,He",
    "evolve_dust": "GROWTH,SPUTTER,EVAP", "gamma": "EXACT",
    "wd01": WD01_DEFAULT, "surface_bins": "10", "cr_desorption": "crate",
    "rtol": "1e-6", "atol": "1e-30", "custom_atol": "atol_species.txt",
    "custom_rtol": "rtol_species.txt",
}
def option_lines(args):
    """the options of args that differ from the defaults, as
    "section.key = value" lines (NETWORK.md); the network first"""
    import config
    import tomlio
    vals, d0 = config.from_namespace(args), config.defaults()
    lines = ["project.network = %s"
             % tomlio.value(str(Path(args.network_file).resolve()))]
    for o in config.OPTIONS:
        v = vals[o.section][o.key]
        if o.name in ("project.network", "project.outdir") or v is None \
                or v == d0[o.section][o.key]:
            continue
        lines.append("%s = %s" % (o.name, tomlio.value(v)))
    return lines


def project_record(args):
    """krome.toml: every option this project was generated with, the
    network as an absolute path, and [meta] (version, the network's
    SHA-256, the test); the project directory is the file's own"""
    import hashlib
    import config
    vals = config.from_namespace(args)
    net = Path(args.network_file).resolve()
    vals["project"]["network"] = str(net)
    vals["project"]["outdir"] = None
    meta = {"version": krome_version(),
            "network_sha256": hashlib.sha256(net.read_bytes()).hexdigest()}
    if args.test:
        meta["test"] = args.test
    return config.to_toml(vals, meta, header=(
        "# The configuration this project was generated with (written by "
        "KROME-X,\n# overwritten on regeneration). Regenerate in place "
        "with\n#     python3 -m kromex krome.toml\n# or, in Python, "
        "kromex.Preprocessor.from_file(\"krome.toml\").run()."))


class GenerationError(Exception):
    """a network or a configuration the generator rejects (the message
    says why); main() prints it, kromex.Preprocessor raises it"""


def main(argv=None):
    ap = build_parser()
    try:
        args = ap.parse_args(sys.argv[1:] if argv is None else list(argv))
    except (ValueError, OSError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        sys.exit(1)
    try:
        res = generate(args, ap)
    except GenerationError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        sys.exit(1)
    if args.check:
        sys.exit(0 if res["check"]["ok"] else 1)


def generate(args, ap=None):
    """generate the project of the parsed arguments args (or, with
    args.check, check the network); GenerationError on a rejected network
    or configuration. Returns {"outdir", "network", "warnings", "check"}"""
    ap = ap or build_parser()
    if args.outdir is None:
        args.outdir = Path("build_%s" % args.test if args.test
                           else "generated")
    args.network_file = args.network_pos or args.network
    if args.network_file is None:
        raise GenerationError("no network given (positional, -n FILE, or "
                              "in the -options file)")

    if args.check:  # nothing is written but the temporary probe
        import shutil
        import tempfile
        check_tmp = Path(tempfile.mkdtemp(prefix="krome_check_"))
    if is_legacy_network(args.network_file):
        args.network_file = convert_legacy_network(
            args.network_file, check_tmp if args.check else args.outdir)
    net = knet.parse_network(args.network_file)
    if net.errors:
        raise GenerationError("\nERROR: ".join(net.errors))
    try:
        gen = Codegen(net, outside=args.outside,
                      secion=args.secondary_ionization.split(","),
                      num2col=args.num2col,
                      metal_cooling=args.metal_cooling.split(","),
                      dust_tables=args.dust_tables,
                      gas_total=args.gas_total,
                      cooling=args.cooling.split(","),
                      heating=args.heating.split(","),
                      cooling_floor=args.cooling_floor,
                      evolve_rad=args.evolve_rad,
                      evolve_Td=args.evolve_Td,
                      shielding=args.shielding,
                      gamma=args.gamma,
                      h2_opacity=args.h2_opacity,
                      photo_induced=args.photo_induced,
                      photoelectric=args.photoelectric, wd01=args.wd01,
                      cmb_pumping=args.cmb_pumping,
                      habing_shield=args.habing_dust_shielding,
                      reactant_multiplicity=args.reactant_multiplicity,
                      screening=args.screening,
                      binding_energies=args.binding_energies,
                      surface_nu=args.surface_nu,
                      surface_competition=not args.no_surface_competition,
                      surface_tunnelling=not args.no_surface_tunnelling,
                      reactive_desorption=args.reactive_desorption,
                      cr_desorption=args.cr_desorption,
                      surface_bins=args.surface_bins,
                      rtol=args.rtol, atol=args.atol,
                      cool_files=[f for f in args.cool_file.split(",")
                                  if f],
                      cool_levels=args.cool_levels,
                      evolve_dust=args.evolve_dust,
                      custom_atol=args.custom_atol,
                      custom_rtol=args.custom_rtol,
                      cool_escape=args.cool_escape,
                      rhs=args.rhs)
    except (ValueError, FileNotFoundError) as exc:
        raise GenerationError(str(exc)) from exc
    info = {"outdir": Path(args.outdir), "network": net.name,
            "warnings": sorted(set(gen.warnings)), "check": None}
    if args.check:
        import netcheck
        try:
            res = netcheck.check(net, gen, check_tmp)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            raise GenerationError(str(exc)) from exc
        finally:
            shutil.rmtree(check_tmp, ignore_errors=True)
        # including those raised while generating the probe's code
        info["warnings"] = sorted(set(gen.warnings))
        netcheck.report(net, res, info["warnings"])
        info["check"] = res
        return info
    args.outdir.mkdir(parents=True, exist_ok=True)
    try:
        if args.single:
            files = {"net_%s.hpp" % net.name: gen.generate(),
                     "krome_network.hpp": gen.interface_files()[
                         "krome_network.hpp"]}
        else:
            files = gen.generate_files(fortran=args.fortran)
    except (ValueError, FileNotFoundError) as exc:
        raise GenerationError(str(exc)) from exc
    if args.check_reverse:
        files["reverse_check.txt"] = gen.reverse_report()
    if not args.single and not args.no_runtime:
        import netcheck  # make check: reference rates of this network
        files["krome_check.cpp"] = netcheck.check_source(net, gen)
    files["NETWORK.md"] = gen.network_md(option_lines(args))
    files["krome.toml"] = project_record(args)
    for name, code in files.items():
        (args.outdir / name).write_text(code)
    extra = ""
    if not args.single and not args.no_runtime:
        extra = "; " + write_project(args.outdir, net.name,
                                     fortran=args.fortran,
                                     python=args.python,
                                     required=gen.param_required_by,
                                     test=args.test)
        if args.test:
            extra += "; " + write_test(args.outdir, args.test)
    print("wrote %s  (%d species, %d reactions%s%s)"
          % (args.outdir / ("net_%s.hpp" % net.name), len(net.species),
             len(net.reactions),
             "" if args.single else "; split into %d files"
             % (len(files) - 2 - ("krome_check.cpp" in files)),
             extra))
    # every warning, including those raised while writing the files
    info["warnings"] = sorted(set(gen.warnings))
    for wmsg in info["warnings"]:
        print("WARNING: " + wmsg)
    if getattr(gen, "burcat_used", False):
        print("NOTE: " + BURCAT_TERMS)
    return info


RUNTIME_DIR = paths.ROOT / "runtime"

MAKEFILE_HEAD = """\
# generated by KROME-X codegen for network '@NET@' -- written once,
# edit freely (regeneration does not overwrite an existing Makefile)
#
# KROME-X is a C/C++ library: the C++ headers can be used directly
# (header-only, plus net_@NET@_data.cpp), and libkrome.a carries the C
# interface (krome.h) for any caller.
#
#   make             libkrome.a + @DEFAULT@@TARGETS@
#   make lib         libkrome.a: the C interface (krome.h) + this network
NET      = @NET@
# make predefines CXX (and FC=f77): replace only the built-in defaults, so
# an environment or command-line choice still wins
ifeq ($(origin CXX),default)
CXX = c++
endif
# production default: -O3, fused multiply-add and the loop SIMD pragmas
# (#pragma omp simd, -fopenmp-simd: no OpenMP runtime). Results differ
# from a strict build at roundoff; for bitwise-reproducible builds use
# CXXFLAGS="-O2 -std=c++17 -ffp-contract=off" (the test suite does).
# Add -march=native (x86) or -mcpu=native (arm64) for the local CPU's
# widest vectors when the binary runs where it is built.
CXXFLAGS ?= -O3 -std=c++17 -ffp-contract=fast -fopenmp-simd
# position-independent objects, so the same objects also make the shared
# library (Python binding)
PIC       = -fPIC
# OpenMP over cells in krome_solve_cells: add -fopenmp to CXXFLAGS (and
# the matching runtime to LIBS; Apple clang needs libomp)
# no external library: the integrator (a port of CVODE's BDF) is part of
# the headers
INC       = -I.
LIBS      =
# LAPACK dense LU (SolverOpts::linsol; auto from 64 equations): a
# blocked, vectorized factorization, recommended for large networks
# (Apple Accelerate: 4x at 145 equations, 13x at 655). Set LAPACK to the
# vendor library, e.g. make LAPACK="-framework Accelerate" (macOS) or
# LAPACK=-lopenblas (single-threaded with OpenMP over cells:
# OPENBLAS_NUM_THREADS=1)
LAPACK   ?=
ifneq ($(strip $(LAPACK)),)
INC      += -DKROME_USE_LAPACK
LIBS     += $(LAPACK)
endif
# device backend (Kokkos: path "device" of krome_solve_cells): KOKKOS = the
# Kokkos prefix, CXX = the compiler Kokkos was built with (nvcc_wrapper for
# CUDA, hipcc for HIP, a host compiler with -fopenmp for OpenMP), e.g.
#   make lib KOKKOS=$$HOME/kokkos CXX=$$HOME/kokkos/bin/nvcc_wrapper
# The device build inlines the big tables (KROME_TABLES_INLINE: device-safe
# accessors from net_$(NET)_tables.hpp); net_$(NET)_data.cpp is then empty.
KOKKOS   ?=
DATAFLAGS =
ifneq ($(strip $(KOKKOS)),)
INC      += -DKROME_USE_KOKKOS -DKROME_TABLES_INLINE -I$(KOKKOS)/include
DATAFLAGS = -DKROME_TABLES_INLINE
KOKKOS_LIB ?= $(firstword $(wildcard $(KOKKOS)/lib64) $(KOKKOS)/lib)
LIBS     += -L$(KOKKOS_LIB) -lkokkoscore -Wl,-rpath,$(KOKKOS_LIB)
endif
# C++ runtime, for linking libkrome.a from C or Fortran
ifeq ($(shell uname),Darwin)
CXXLIB   ?= -lc++
SOEXT     = dylib
else
CXXLIB   ?= -lstdc++
SOEXT     = so
endif

all: @ALL@

# the big embedded tables: compiled once, relinked as needed
net_$(NET)_data.o: net_$(NET)_data.cpp net_$(NET)_thermal.hpp
\t$(CXX) $(CXXFLAGS) $(PIC) $(DATAFLAGS) -I. -c $< -o $@

HDRS      = krome.hpp krome_network.hpp net_$(NET)_core.hpp \
            $(wildcard core/*.hpp) $(wildcard core/cvode/*.hpp)

@EXAMPLE@
# the C interface library
krome_c.o: krome.cpp krome.h $(HDRS)
\t$(CXX) $(CXXFLAGS) $(PIC) $(INC) -c krome.cpp -o $@

libkrome.a: krome_c.o net_$(NET)_data.o
\tar rcs $@ $^

lib: libkrome.a
"""

MAKEFILE_FORTRAN = """
# Fortran binding (codegen project.bindings = ['fortran']): modules krome + krome_net
ifeq ($(origin FC),default)
FC = gfortran
endif
FFLAGS   ?= -O2

krome_net.mod: krome_net.f90
\t$(FC) $(FFLAGS) -c krome_net.f90

krome.mod: krome.f90 krome_net.mod
\t$(FC) $(FFLAGS) -c krome.f90

fortran: libkrome.a krome.mod
@EXAMPLE_F90@"""

MAKEFILE_EXAMPLE = """example: example.cpp net_$(NET)_data.o $(HDRS)
\t$(CXX) $(CXXFLAGS) $(INC) example.cpp net_$(NET)_data.o $(LIBS) -o $@
"""

MAKEFILE_EXAMPLE_F90 = """
example_f90: example.f90 libkrome.a krome.mod
\t$(FC) $(FFLAGS) example.f90 krome.o krome_net.o libkrome.a $(LIBS) \\
\t  $(CXXLIB) -o $@
"""

MAKEFILE_PYTHON = """
# the shared library: the Python binding krome.py (codegen project.bindings = ['python'], or
# kromex.Preprocessor.build) loads it
libkrome.$(SOEXT): krome_c.o net_$(NET)_data.o
\t$(CXX) $(CXXFLAGS) -shared krome_c.o net_$(NET)_data.o $(LIBS) -o $@

shared: libkrome.$(SOEXT)
"""

MAKEFILE_TEST = """
# the test driver (python3 -m kromex --test @TEST@): make && ./test
test: test.cpp krome_test_utils.hpp net_$(NET)_data.o $(HDRS)
\t$(CXX) $(CXXFLAGS) $(INC) test.cpp net_$(NET)_data.o $(LIBS) -o $@
"""

MAKEFILE_TAIL = """
# make check: build every target, then compare the rate coefficients of
# this build (compiler and flags) with reference values computed when the
# project was generated (krome_check.cpp, an independent evaluation)
krome_check: krome_check.cpp net_$(NET)_data.o $(HDRS)
\t$(CXX) $(CXXFLAGS) $(INC) krome_check.cpp net_$(NET)_data.o $(LIBS) -o $@

check: @CHECKDEPS@ krome_check
\t./krome_check

clean:
\trm -f krome_check @CLEAN@*.o *.mod libkrome.a \\
\t  libkrome.so libkrome.dylib
"""


def makefile(netname, fortran=False, python=False, test=None):
    """the project's Makefile: with --test the default target builds the
    test driver (and there is no example), else the example"""
    targets = ""
    alls = ["lib", "test" if test else "example"]
    body = MAKEFILE_HEAD.replace("@EXAMPLE@", "" if test else
                                 MAKEFILE_EXAMPLE)
    if test:
        default = "the test driver (test.cpp)"
    else:
        default = "the C++ example (example.cpp)"
        targets += "\n#   make example     the C++ example"
    if fortran:
        targets += ("\n#   make fortran     + the Fortran modules krome, "
                    "krome_net")
        if not test:
            targets += "\n#   make example_f90 the Fortran example"
        body += MAKEFILE_FORTRAN.replace(
            "@EXAMPLE_F90@", "" if test else MAKEFILE_EXAMPLE_F90)
    if test:
        body += MAKEFILE_TEST.replace("@TEST@", test)
    # the shared library is always a target (a Python binding added later,
    # e.g. by kromex.Preprocessor.build, needs no new Makefile); with
    # --python make builds it too
    targets += ("\n#   make shared      libkrome.so/.dylib for the "
                "Python module krome.py")
    body += MAKEFILE_PYTHON
    if python:
        alls.append("shared")
    body += MAKEFILE_TAIL.replace("@CLEAN@", "test " if test
                                  else "example example_f90 ")
    targets += ("\n#   make check       build everything and check the "
                "rates of this build\n#                    against "
                "reference values (krome_check.cpp)")
    body = body.replace("@CHECKDEPS@", " ".join(
        alls + (["fortran"] + ([] if test else ["example_f90"])
                if fortran else [])))
    return (body.replace("@TARGETS@", targets).replace("@DEFAULT@", default)
            .replace("@ALL@", " ".join(alls)).replace("@NET@", netname))


EXAMPLE = """\
// KROME-X example driver for network '@NET@' -- written once, edit
// freely (regeneration does not overwrite an existing example.cpp).
#include <cstdio>
#include "krome.hpp"  // the C++ API: solve(), solve_cells(), queries, ...
using namespace krome;

int main() {
  double n[NSPEC];
  for (int i = 0; i < NSPEC; ++i) n[i] = 1e-20;
  // TODO: set the initial abundances [cm^-3], e.g.:
  // n[IDX_H] = 1e4; n[IDX_E] = n[IDX_Hj] = 1e0;
  double Tgas = 1e2;

  Params p;  // runtime parameters (listed in NETWORK.md)
@REQ_CPP@
  SolverOpts o;         // solver options (see krome.hpp)
  o.evolve_T = true;    // false: fixed-Tgas chemistry

  const int flag = solve(n, Tgas, p, 1e10, o);
  printf("flag %d (%s), Tgas = %g K\\n", flag, flag_message(flag), Tgas);
  for (int i = 0; i < NSPEC; ++i)
    printf("%-10s %.6e\\n", species_names[i], n[i]);
  return 0;
}
"""


EXAMPLE_F90 = """\
! KROME-X Fortran example for network '@NET@' -- written once, edit
! freely (regeneration does not overwrite an existing example.f90).
! Build: make example_f90
program example
  use krome
  use krome_net
  implicit none
  integer, parameter :: ncell = 4
  real(8) :: n(krome_nspec, ncell), Tgas(ncell), dt(ncell)
  real(8) :: par(krome_npar, ncell)
  integer(c_int) :: flags(ncell)
  type(c_ptr) :: h
  integer :: c, ierr

  n = 1d-20
  ! TODO: set the initial abundances [cm^-3], e.g.
  ! n(krome_idx_H, :) = 1d4
  Tgas = 1d2
  dt = 1d10

  call krome_params_default(par)   ! then per-cell edits
@REQ_F90@
  h = krome_opts_new()
  call krome_set(h, "evolve_T", 1d0)

  ierr = krome_solve(n, Tgas, dt, h, par=par, flags=flags)
  print *, "ierr", ierr, " flags", flags
  do c = 1, ncell
    print *, "cell", c, " Tgas", Tgas(c)
  end do
  call krome_opts_free(h)
end program example
"""


EXAMPLE_PY = """\
# KROME-X Python example for network '@NET@' -- written once, edit
# freely (regeneration does not overwrite an existing example.py).
# Build the library first: make shared   (needs numpy at runtime)
import krome

net = krome.load()        # libkrome next to krome.py (or $KROME_LIB)
print(net)

# TODO: set the initial abundances [cm^-3], e.g.
# n = net.densities({"H": 1e4, "H+": 1e-4, "e-": 1e-4})
n = net.densities({}, fill=1e-20)

opts = net.options(evolve_T=True)

# per-cell parameters by name (see net.params)
par = {}
@REQ_PY@
res = net.solve(n, Tgas=1e2, dt=1e10, par=par, opts=opts)
print("Tgas = %g K" % res.Tgas)
for s, v in zip(net.species, res.n):
    print("%-10s %.6e" % (s, v))
"""


def _required_lines(required, fmt, head):
    if not required:
        return ""
    return "\n".join([head] + [fmt % (nm, by) for nm, by
                               in sorted(required.items())])


def write_project(outdir, netname, fortran=False, python=False,
                  required=None, test=None):
    """copy the pre-cooked runtime files into the project directory and
    write the Makefile and examples once (user files are never
    overwritten). The C/C++ library is always there; the Fortran and
    Python bindings only on request."""
    import shutil
    ncopy = 0
    # the C++ API, and the interfaces (runtime/interfaces/: the C one
    # always, the bindings on request), all at the project's top
    ifc = RUNTIME_DIR / "interfaces"
    srcs = sorted(RUNTIME_DIR.glob("*.hpp"))
    srcs += [f for pat in ("*.h", "*.cpp")
             for f in sorted((ifc / "C").glob(pat))]
    if fortran:
        srcs += sorted((ifc / "fortran").glob("*.f90"))
    if python:
        srcs += sorted((ifc / "python").glob("*.py"))
    for f in srcs:
        shutil.copy2(f, outdir / f.name)
        ncopy += 1
    # the implementation headers (krome.hpp includes them) and, in
    # core/cvode/, the integrator (CVODE's BDF, ported) with its SUNDIALS
    # license. core/ is KROME's: headers the runtime no longer has
    # (renamed or removed) go, and so do the earlier layouts' detail/ and
    # cvode/ (only KROME's files; a directory left empty is removed)
    for sub in ("core", "core/cvode"):
        (outdir / sub).mkdir(exist_ok=True)
        keep = {f.name for f in (RUNTIME_DIR / sub).glob("*.hpp")}
        for f in (outdir / sub).glob("*.hpp"):
            if f.name not in keep:
                f.unlink()
        for f in sorted((RUNTIME_DIR / sub).glob("*")):
            if f.is_file():
                shutil.copy2(f, outdir / sub / f.name)
                ncopy += 1
    for old in ("detail", "cvode"):
        d = outdir / old
        if d.is_dir():
            for f in list(d.glob("*.hpp")) + list(d.glob("LICENSE.sundials")):
                f.unlink()
            if not any(d.iterdir()):
                d.rmdir()
    made, notes = [], []
    req = required or {}
    # with --test the test driver is the project's program: no examples
    todo = [("Makefile", makefile(netname, fortran, python, test))]
    if not test:
        todo.append(("example.cpp", EXAMPLE.replace("@NET@", netname).replace(
            "@REQ_CPP@", _required_lines(
                req, "  // p.%s = ...;  // %s",
                "  // TODO: parameters without a default (the solver "
                "rejects them unset):"))))
    if fortran and not test:
        todo.append(("example.f90", EXAMPLE_F90.replace(
            "@NET@", netname).replace("@REQ_F90@", _required_lines(
                req, "  ! par(krome_par_%s, :) = ...  ! %s",
                "  ! TODO: parameters without a default (the solver "
                "rejects them unset):"))))
    if python and not test:
        todo.append(("example.py", EXAMPLE_PY.replace(
            "@NET@", netname).replace("@REQ_PY@", _required_lines(
                req, "# par[\"%s\"] = ...  # %s",
                "# TODO: parameters without a default (net.required; the "
                "solver rejects them unset):"))))
    for name, text in todo:
        path = outdir / name
        if not path.exists():
            path.write_text(text)
            made.append(name)
    mk = (outdir / "Makefile").read_text()
    for flag, target in ((fortran, "fortran:"), (python, "shared:"),
                         (test, "test:")):
        if flag and target not in mk:
            notes.append("the existing Makefile has no '%s' target (remove "
                         "it to regenerate)" % target[:-1])
    return "%d runtime files copied%s%s" % (
        ncopy, (", wrote " + " + ".join(made)) if made else "",
        "".join("; NOTE: " + x for x in notes))


def write_test(outdir, name):
    """the test's files (examples/tests/common + examples/tests/NAME: the
    driver, the plot script, data; never overwriting)"""
    import shutil
    made = []
    for src in (TESTS_DIR / "common", TESTS_DIR / name):
        for f in sorted(src.iterdir()) if src.is_dir() else []:
            if f.name == "options.toml" or not f.is_file():
                continue
            if not (outdir / f.name).exists():
                shutil.copy2(f, outdir / f.name)
                made.append(f.name)
    return ("test %s: %s; build and run: make && ./test, plot: "
            "python3 plot.py" % (name, " + ".join(made) if made else
                                 "files kept"))

if __name__ == "__main__":
    main()
