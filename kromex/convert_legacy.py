#!/usr/bin/env python3
"""Convert a legacy KROME network file (react_*) to the draft KROME-X format (.knet).

Usage:
    python3 convert_legacy.py <legacy_file> -o <output.knet> [--check] [--quiet]

Commented-out legacy reaction blocks are converted too (kept as comments, in
new syntax). --check numerically compares each plain analytic rate
before/after conversion over its declared temperature range (photo/xray/
surface/auto/extern/density-dependent rates are skipped).
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402
import constants  # noqa: E402

# ----------------------------------------------------------------------------
# lexical helpers
# ----------------------------------------------------------------------------

ELEMENTS = [
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar", "Ca", "Sc",
    "Ti", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "H", "D", "C", "N", "O",
    "F", "P", "S", "K", "V",
]
ELEMENTS_BY_LEN = sorted(ELEMENTS, key=len, reverse=True)

BUILTIN_FUNCS = {"exp", "log", "log10", "sqrt", "abs", "min", "max", "tanh",
                 "where"}


def f90_ops(e: str) -> str:
    """Fortran comparison operators -> symbolic (.gt. -> >, ...)."""
    for f, p in ((r"\.ge\.", ">="), (r"\.le\.", "<="), (r"\.gt\.", ">"),
                 (r"\.lt\.", "<"), (r"\.eq\.", "=="), (r"\.ne\.", "!=")):
        e = re.sub(f, p, e, flags=re.I)
    return e


def if_to_where(e: str) -> str:
    """Legacy conditional rate 'if(cond): expr' -> 'where(cond, expr, 0.0)'."""
    m = re.match(r"^\s*if\s*\((.+)\)\s*:\s*(.+)$", e, flags=re.S)
    return "where(%s, %s, 0.0)" % (m.group(1), m.group(2)) if m else e


def merge_to_where(e: str) -> str:
    """Fortran merge(a, b, mask) -> where(mask, a, b)."""
    while True:
        m = re.search(r"\bmerge\s*\(", e)
        if not m:
            return e
        j = e.index("(", m.start())
        depth, k = 0, j
        while k < len(e):
            if e[k] == "(":
                depth += 1
            elif e[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        a, b, mask = (x.strip() for x in smart_split(e[j + 1:k])[:3])
        e = e[:m.start()] + "where(%s, %s, %s)" % (mask, a, b) + e[k + 1:]


def normalize_expr(e: str) -> str:
    e = d2e(unicode_ascii(e.strip()))
    # Fortran kind suffixes on literals: 0.0_8 -> 0.0 (never species labels)
    e = re.sub(r"(?<![A-Za-z_\d])(\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?)_(?:4|8|16)\b",
               r"\1", e)
    return if_to_where(merge_to_where(f90_ops(e)))


def lazy_where(e: str) -> str:
    """where(c, a, b) -> Python conditional expression (lazy, checker only)."""
    while True:
        m = re.search(r"\bwhere\s*\(", e)
        if not m:
            return e
        j = e.index("(", m.start())
        depth, k = 0, j
        while k < len(e):
            if e[k] == "(":
                depth += 1
            elif e[k] == ")":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        c, a, b = (x.strip() for x in smart_split(e[j + 1:k])[:3])
        e = e[:m.start()] + "((%s) if (%s) else (%s))" % (a, c, b) + e[k + 1:]
CONSTRUCTORS = {"photo", "surface"}

# legacy shortcut symbols -> KROME-X definitions (emitted as `var` preamble if
# used; full list from legacy kromelib.get_Tshortcut, matched case-insensitively
# like the legacy preprocessor does)
SHORTCUTS = {
    "Te":         "Tgas_eV",
    "lnTe":       "log(Tgas_eV)",
    "invTe":      "1.0/Tgas_eV",
    "logT":       "log10(Tgas)",
    "lnT":        "log(Tgas)",
    "invT":       "1.0/Tgas",
    "invTgas":    "1.0/Tgas",
    "sqrTgas":    "sqrt(Tgas)",
    "invsqrT":    "1.0/sqrt(Tgas)",
    "T32":        "Tgas/300.0",
    "t3":         "Tgas/300.0",
    "t4":         "Tgas*1.0e-4",
    "invT32":     "300.0/Tgas",
    "sqrT32":     "sqrt(Tgas/300.0)",
    "invsqrT32":  "sqrt(300.0/Tgas)",
    "Tgas2":      "Tgas*Tgas",
    "Tgas3":      "Tgas*Tgas*Tgas",
    "Tgas4":      "Tgas*Tgas*Tgas*Tgas",
    "T0":         "288.0",
    "T02":        "288.0**2",
    "T03":        "288.0**3",
    "T04":        "288.0**4",
    "T0inv":      "1.0/288.0",
}
SHORTCUT_ORDER = ("Te", "lnTe", "invTe", "logT", "lnT", "invT", "invTgas",
                  "sqrTgas", "invsqrT", "T32", "t3", "t4", "invT32", "sqrT32",
                  "invsqrT32", "Tgas2", "Tgas3", "Tgas4", "T0", "T02", "T03",
                  "T04", "T0inv")

# legacy derived constants (kromeobj krome_constants) that rate expressions
# may name: as an expression of the network-language constants, or as the
# value from kromex/constants.py (CODATA 2018) when it needs m_p, h or c
_C = constants
LEGACY_CONSTANTS = {
    "pre_kvgas_sqrt": "sqrt(8.0*boltzmann_erg/pi)",
    "kvgas_erg": repr(8.0 * _C.boltzmann_erg / _C.pi / _C.p_mass),
    "ip_mass": repr(1.0 / _C.p_mass),
    "pre_planck": repr(2.0 * _C.planck_erg / _C.clight ** 2),
    "exp_planck": repr(_C.planck_erg / _C.boltzmann_erg),
}

# implicit-state arguments stripped from extension-function calls
STATE_ARGS = {"n", "n(:)", "tgas", "tdust", "tabtdust"}

KB_EV = constants.boltzmann_eV  # [eV/K] (kromex/constants.py)


def d2e(s: str) -> str:
    """Fortran d-exponents -> e-exponents (1.3d-32 -> 1.3e-32)."""
    return re.sub(r"(?<=[0-9.])[dD](?=[+\-]?[0-9])", "e", s)


def smart_split(s: str, sep=","):
    """Split on `sep` at zero parenthesis/brace depth."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch in "({":
            depth += 1
        elif ch in ")}":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def canon_species(name: str):
    """Canonicalize a legacy species name; None means 'drop' (photon)."""
    s = name.strip()
    if not s:
        return None
    # legacy isotopes in brackets: [60Fe] -> 60Fe, [13C]O -> 13CO (KROME-X
    # writes the mass number first)
    s = re.sub(r"\[(\d+)([A-Z][a-z]?)\]", r"\1\2", s)
    # excited-state labels: legacy O(1D) -> O_1D (as legacy f90 names did)
    s = s.replace("(", "_").replace(")", "").replace("[", "_").replace("]", "")
    # interior charge signs are invalid in KROME-X: OCS+H2 -> OCSH2+
    body = s.rstrip("+-")
    if "+" in body or "-" in body:
        trail = s[len(body):]
        s = (body.replace("+", "").replace("-", "") + trail
             + "+" * body.count("+") + "-" * body.count("-"))
    if s in ("g", "G", "photon", "CR"):
        return None
    if s.upper() in ("E", "E-"):
        return "e-"
    if s.startswith("{") and s.endswith("}"):
        inner = canon_species(s[1:-1])
        return "{%s}" % inner if inner else None
    m = re.match(r"^(\d*)(.*?)(\++|-+)?$", s)
    iso, core, charge = m.group(1), m.group(2), m.group(3) or ""
    # phase suffix (_ice/_dust and legacy variants) preserved verbatim
    phase = ""
    pm = re.match(r"^(.*?)(_\w+)$", core)
    if pm:
        core, phase = pm.group(1), pm.group(2)
    # the elements as written (CO carbon monoxide, Co cobalt, SC, Sc), else
    # case-insensitively as legacy (Semenov networks write FE, MG, SI)
    out = _segment(core, True)
    if out is None:
        out = _segment(core, False)
    return iso + out + phase + charge


def _segment(core, exact):
    """a species core as element symbols (with counts); exact: the case as
    written, None if it does not segment; else case-insensitive, an
    unknown remainder kept verbatim"""
    out, i = "", 0
    while i < len(core):
        for el in ELEMENTS_BY_LEN:
            seg = core[i:i + len(el)]
            if seg == el if exact else seg.lower() == el.lower():
                out += el
                i += len(el)
                while i < len(core) and core[i].isdigit():
                    out += core[i]
                    i += 1
                break
        else:
            if exact:
                return None
            out += core[i:]  # unknown token: keep verbatim
            break
    return out


IDX_MARK = re.compile(r"\x01(\w+)\x02")


def decode_idx(name: str) -> str:
    """A legacy index name (idx_<name>) as a placeholder, resolved after the
    whole file is read (Converter.resolve_idx): against the network's own
    species by legacy's rule for the f90 names, case-insensitively as
    Fortran, so labelled ions (idx_H3j_para -> H3_para+) and any spelling
    case (idx_H2_ORTHO -> H2_ortho) resolve; guess_idx otherwise."""
    return "\x01%s\x02" % name


def legacy_f90idx(raw: str) -> str:
    """legacy KROME's index name of a species as written in the network
    (kromelib: idx_ + name, + -> j, - -> k, ( -> _, ) and [ dropped,
    ] -> _), lower-cased (Fortran is case-insensitive)"""
    return (raw.replace("+", "j").replace("-", "k").replace("(", "_")
            .replace(")", "").replace("[", "").replace("]", "_").lower())


def guess_idx(name: str) -> str:
    """Legacy Fortran-safe index name -> species (idx_Hj -> H+, idx_Dk -> D-,
    idx_E -> e-); charges >1 are encoded as j2/k2 etc."""
    if name.upper() == "E":
        return "e-"
    m = re.match(r"^(.+?)j(\d*)$", name)
    if m:
        return canon_species(m.group(1) + "+" * int(m.group(2) or 1))
    m = re.match(r"^(.+?)k(\d*)$", name)
    if m:
        return canon_species(m.group(1) + "-" * int(m.group(2) or 1))
    return canon_species(name)


def unicode_ascii(e: str) -> str:
    """typographic characters pasted into legacy rate expressions (minus
    and dashes, which gfortran rejects) -> ASCII"""
    for u, a in (("\u2212", "-"), ("\u2013", "-"), ("\u2014", "-"),
                 ("\u00d7", "*"), ("\u00a0", " ")):
        e = e.replace(u, a)
    return e


def parse_limit(tok: str) -> str:
    """'.LE.5.5e3' / '>1160' / 'NONE' -> '5.5e3' / '1160' / ''."""
    t = tok.strip()
    if not t or t.upper() == "NONE":
        return ""
    t = re.sub(r"^(\.GE\.|\.LE\.|\.GT\.|\.LT\.|>=|<=|>|<)", "", t, flags=re.I)
    return d2e(t.strip())


DEFAULT_FMT = ["idx", "R", "R", "R", "P", "P", "P", "P", "Tmin", "Tmax", "rate"]
FMT_MAP = {"shield": "shield",
           "idx": "idx", "r": "R", "p": "P", "tmin": "Tmin", "tmax": "Tmax",
           "rate": "rate", "k": "rate", "qeff": "qeff", "qpp": "qeff"}


class BlockState:
    """Parser state for one stream (active lines, or commented-out lines)."""

    def __init__(self):
        self.fmt = list(DEFAULT_FMT)
        self.photo = False
        self.xray = False
        self.surface = False
        self.cooling = False   # inside @cooling/@heating block
        self.cool_kind = "cooling"


# ----------------------------------------------------------------------------
# converter
# ----------------------------------------------------------------------------

class Converter:
    def __init__(self):
        self.used_shortcuts = set()
        self.externs = set()
        self.semenov_mapped = 0  # legacy grain functions -> surface(...)
        self.grains = None     # (n, a2) of the legacy single-size grains
        self.params = {}       # new name -> legacy name
        self.tabvars = set()   # names defined via @tabvar (not params)
        self.items = []        # ordered output items (dicts)
        self.checks = []       # (lineno, reaction, old_expr, new_expr, lo, hi,
        #                          legacy k() index or None)
        self.kold = {}         # legacy k() index -> (old_expr, lo, hi)
        self.warnings = []
        self.n_reactions = 0
        self.n_commented = 0
        self.n_cool = 0
        self.uses_bare_T = False
        self.fake_species = set()  # legacy hardcoded fakes (FK1, FK2, ...)
        self.photo_auto_eqs = set()  # dedup photo vs legacy XRAY blocks
        self.hooks = []        # legacy modifier blocks {kind, body, ...}
        self.kseq = 0          # legacy k() numbering of active reactions
        self.idx_map = {}      # legacy f90 index name (lower) -> species

    # -- expression handling --------------------------------------------------

    def clean_extern_args(self, e: str) -> str:
        """Strip implicit-state args and idx_X -> X in non-builtin calls."""
        out, i = "", 0
        while i < len(e):
            m = re.match(r"[A-Za-z_]\w*", e[i:])
            if m:
                name = m.group(0)
                j = i + len(name)
                if (j < len(e) and e[j] == "(" and name not in BUILTIN_FUNCS
                        and name != "n"):
                    depth, k = 0, j
                    while k < len(e):
                        if e[k] == "(":
                            depth += 1
                        elif e[k] == ")":
                            depth -= 1
                            if depth == 0:
                                break
                        k += 1
                    inner = self.clean_extern_args(e[j + 1:k])
                    kept = []
                    for a in (x.strip() for x in smart_split(inner)):
                        if a.replace(" ", "").lower() in STATE_ARGS:
                            continue
                        kept.append(re.sub(r"\bidx_(\w+)\b",
                                           lambda m: decode_idx(m.group(1)),
                                           a))
                    out += name + "(" + ", ".join(kept) + ")"
                    i = k + 1
                    continue
                out += name
                i = j
                continue
            out += e[i]
            i += 1
        return out

    def convert_expr(self, e: str, register=True) -> str:
        e = normalize_expr(e)  # d->e literals, .gt.-> >, if(): -> where()
        e = re.sub(r"\bsum\s*\(\s*n\s*\(\s*1\s*:\s*nmols\s*\)\s*\)", "ntot", e)
        e = re.sub(r"\bdexp\b", "exp", e)
        for c, v in LEGACY_CONSTANTS.items():
            e = re.sub(r"\b%s\b" % c, "(%s)" % v, e)
        e = re.sub(r"\benergy_eV\b", "E_eV", e)
        e = re.sub(r"\buser_Tdust\b", "Tdust", e)        # builtin, not a param
        e = re.sub(r"\bphys_(zredshift|Tcmb)\b", r"\1", e)  # KROME-X params
        e = re.sub(r"\bget_Hnuclei\s*\([^()]*(?:\([^()]*\)[^()]*)*\)",
                   "Hnuclei", e)  # builtin state symbol in KROME-X
        e = re.sub(r"\bn\s*\(\s*idx_(\w+)\s*\)",
                   lambda m: "n(%s)" % decode_idx(m.group(1)), e)
        if register:
            for legacy in re.findall(r"\buser_(\w+)\b", e):
                if legacy not in self.tabvars:
                    self.params.setdefault(legacy, "user_" + legacy)
        e = re.sub(r"\buser_(\w+)\b", r"\1", e)
        e = self.clean_extern_args(e)
        # legacy krate_evap_semenov_nonthermal(idx) (the network's
        # user_crflux): the cosmic-ray term of krate_nonthermal_evaporation,
        # i.e. that function with no photodesorption (G0 = yield = 0)
        if re.search(r"\bkrate_evap_semenov_nonthermal\s*\(\s*[^,()]+\)",
                     e):
            e = re.sub(r"\bkrate_evap_semenov_nonthermal\s*\(\s*"
                       r"([^,()]+?)\s*\)",
                       r"krate_nonthermal_evaporation(\1, 0.0, 0.0, crflux, "
                       r"0.0)", e)
            if register:
                self.params.setdefault("crflux", "user_crflux")
                self.externs.discard("krate_evap_semenov_nonthermal")
                self.externs.add("krate_nonthermal_evaporation")
        if re.search(r"\bT\b", e):
            self.uses_bare_T = True
        for name in SHORTCUTS:
            pat = r"\b%s\b" % re.escape(name)
            if re.search(pat, e, flags=re.I):
                e = re.sub(pat, name, e, flags=re.I)  # normalize case
                self.used_shortcuts.add(name)
        if register:
            for fn in re.findall(r"\b([A-Za-z_]\w*)\s*\(", e):
                if fn not in BUILTIN_FUNCS and fn != "n" \
                        and fn not in CONSTRUCTORS:
                    self.externs.add(fn)
        return e

    # -- directive handling ---------------------------------------------------

    def handle_directive(self, s, st, commented, lineno):
        low = s.lower()
        emit = lambda **kw: self.items.append(dict(commented=commented, **kw))

        if low.startswith("@format:"):
            toks = [t.strip().lower() for t in s.split(":", 1)[1].split(",")]
            st.fmt = [FMT_MAP.get(t, t) for t in toks]
        elif low.startswith("@var:") or low.startswith("@var "):
            body = s.split(":", 1)[1] if ":" in s else s[4:]
            name, expr = body.split("=", 1)
            name = name.strip()
            if "(" in name:
                # legacy per-species array var, e.g. imsqrt(nspec)=...;
                # its uses imsqrt(idx_X) become extension-function calls
                emit(t="raw", text="# [legacy array var %s dropped: its "
                                   "uses convert to an extension function]"
                                   % name)
                return
            if name.startswith("user_"):  # match the rename inside exprs
                name = name[5:]
            expr = self.convert_expr(expr, register=not commented)
            if name == expr:  # e.g. legacy 'ntot = sum(n(1:nmols))' -> builtin
                emit(t="raw", text="# [legacy 'var %s' dropped: builtin in KROME-X]"
                                   % name)
            else:
                emit(t="var", name=name, expr=expr)
        elif low.startswith("@common") or low.startswith("@commons"):
            for tok in s.split(":", 1)[1].split(","):
                tok = tok.strip()
                if not tok or tok == "user_Tdust":
                    continue
                new = tok[5:] if tok.startswith("user_") else tok
                self.params.setdefault(new, tok)
        elif low.startswith("@ref:") or low.startswith("@ref "):
            emit(t="ref", text=s.split(":", 1)[1].strip() if ":" in s
                 else s[4:].strip())
        elif low.startswith("@ghost:"):
            names = [canon_species(x) for x in s.split(":", 1)[1].split(",")]
            emit(t="species", names=[x for x in names if x])
        elif low.startswith("@ice:"):
            toks = smart_split(s.split(":", 1)[1])
            if len(toks) < 3:
                self.warnings.append("line %d: malformed @ice" % lineno)
                emit(t="raw", text="# [UNCONVERTED] " + s)
                return
            sp = canon_species(toks[0])
            kind = toks[1].strip().lower()
            rate = self.convert_expr(",".join(toks[2:]),
                                     register=not commented)
            ice = sp + "_ice"
            rx = ("%s -> %s" % (sp, ice) if kind == "freezeout"
                  else "%s -> %s" % (ice, sp))
            self.n_reactions += not commented
            self.n_commented += commented
            emit(t="rx", rx=rx, rng="[]", rate=rate, ann="")
        elif low.startswith("@tabvar:"):
            body = s.split(":", 1)[1]
            name, rest = body.split("=", 1)
            name = name.strip()
            new = name[5:] if name.startswith("user_") else name
            self.tabvars.add(new)
            toks = smart_split(rest)
            fname = toks[0].strip()
            args = [self.convert_expr(a, register=not commented)
                    for a in toks[1:]]
            emit(t="table", name=new + "_tab", file=fname)
            emit(t="var", name=new,
                 expr="%s_tab(%s)" % (new, ", ".join(args)))
        elif re.match(r"@photo_?(start|begin)", low):
            st.photo = True
        elif re.match(r"@photo_?(stop|end)", low):
            st.photo = False
        elif re.match(r"@xray_?(start|begin)", low):
            st.xray = True
            emit(t="raw", text="# [legacy @XRAY block: converted to photo() "
                               "-- KROME-X does X-ray chemistry through the photo "
                               "machinery with bins extending to X-ray "
                               "energies; tabulated X-ray rates are retired]")
        elif re.match(r"@xray_?(stop|end)", low):
            st.xray = False
        elif re.match(r"@surface_?(start|begin)", low):
            st.surface = True
        elif re.match(r"@surface_?(stop|end)", low):
            st.surface = False
        elif re.match(r"@(cooling|heating)_?(start|begin)", low):
            st.cooling = True
            st.cool_kind = "cooling" if "cooling" in low else "heating"
        elif re.match(r"@(cooling|heating)_?(stop|end)", low):
            st.cooling = False
        elif low.startswith("@cooling:") or low.startswith("@heating:"):
            kind = "cooling" if low.startswith("@cooling") else "heating"
            self.n_cool += 1
            emit(t="cool", kind=kind, name="custom_%d" % self.n_cool,
                 expr=self.convert_expr(s.split(":", 1)[1],
                                        register=not commented))
        elif low.startswith("@next_h2photodissociation"):
            # legacy: the next reaction's rate is replaced by kpd_H2 (the
            # H2pdB.dat lines on the photon bins): KROME-X h2pd()
            st.next_h2pd = True
            emit(t="raw", text="# [legacy %s: the next reaction is h2pd()]"
                               % s)
        elif low == "@next_solomon":
            # legacy: the next reaction gets the Solomon pumping heat of
            # -heating=CHEM; KROME-X identifies the Solomon process itself (an
            # H2_solomonLW / J_LW rate, or a photo() H2 -> H + H in the
            # Lyman-Werner band) and gives it that heat
            emit(t="raw", text="# [legacy %s: KROME-X recognizes the Solomon "
                               "reaction itself (Solomon heating of CHEM)]"
                               % s)
        elif re.match(r"@(cr_|catalysis)", low):
            emit(t="raw", text="# [legacy %s marker removed: these are "
                               "ordinary reactions in KROME-X]" % s)
        elif re.match(r"@(notab|storeonce)", low):
            emit(t="raw", text="# [legacy %s dropped: handled outside the "
                               "network file in KROME-X]" % s)
        elif re.match(r"@(reactionmodifier|odemodifier)_?(start|begin)", low):
            st.hook = low.split("_")[0].lstrip("@")
            emit(t="raw", text="# [UNCONVERTED legacy %s block (raw "
                               "Fortran, no KROME-X equivalent):]" % st.hook,
                 hook_head=len(self.hooks))
            self.hooks.append({"kind": st.hook, "body": [],
                               "commented": commented})
            st.in_hook = True
        elif re.match(r"@(reactionmodifier|odemodifier)_?(stop|end)", low):
            st.in_hook = False
        else:
            self.warnings.append("line %d: unhandled directive %r "
                                 "passed through as comment" % (lineno, s))
            emit(t="raw", text="# [UNCONVERTED] " + s)

    # -- reaction handling ----------------------------------------------------

    def looks_like_limit(self, tok):
        t = tok.strip()
        if t.upper() == "NONE" or re.match(r"^(\.[GL][ET]\.|>=?|<=?)", t):
            return True
        try:
            float(d2e(t))
            return True
        except ValueError:
            return False

    def handle_reaction(self, s, st, commented, lineno):
        fields = smart_split(s)
        # some legacy files declare a @format without Tmin,Tmax but write
        # reactions that carry them anyway; the legacy parser tolerated it
        if (len(fields) == len(st.fmt) + 2 and "Tmin" not in st.fmt
                and self.looks_like_limit(fields[-3])
                and self.looks_like_limit(fields[-2])):
            st.fmt = st.fmt[:-1] + ["Tmin", "Tmax", "rate"]
            if not commented:
                self.warnings.append(
                    "line %d: legacy @format lacks Tmin,Tmax but reactions "
                    "carry them — columns inferred" % lineno)
        if len(fields) > len(st.fmt):
            fields = fields[:len(st.fmt) - 1] + \
                     [",".join(fields[len(st.fmt) - 1:])]
        if len(fields) != len(st.fmt):
            if not commented:
                self.warnings.append(
                    "line %d: %d fields but format has %d — kept as comment"
                    % (lineno, len(fields), len(st.fmt)))
                self.items.append({"t": "raw", "commented": False,
                                   "text": "# [UNCONVERTED] " + s})
            else:
                self.items.append({"t": "raw", "commented": False,
                                   "text": "#" + s})
            return

        rea, pro, lo, hi, rate, qeff = [], [], "", "", "", ""
        fshield = ""  # @format 'shield' column (KROME branches that read it)
        for key, val in zip(st.fmt, fields):
            if key == "R":
                sp = canon_species(val)
                if sp:
                    rea.append(sp)
                    self.idx_map.setdefault(legacy_f90idx(val.strip()), sp)
            elif key == "P":
                sp = canon_species(val)
                if sp:
                    pro.append(sp)
                    self.idx_map.setdefault(legacy_f90idx(val.strip()), sp)
            elif key == "Tmin":
                lo = parse_limit(val)
            elif key == "Tmax":
                hi = parse_limit(val)
            elif key == "rate":
                rate = val.strip()
            elif key == "qeff":
                q = val.strip()
                if q:
                    qeff = self.convert_expr(q, register=not commented)
            elif key == "shield":
                f = val.split("!")[0].split("#")[0].strip()
                if f.upper() not in ("", "N", "NONE", "N/A", "NO"):
                    fshield = self.convert_expr(f, register=not commented)

        if not rea and not pro:
            # legacy dummy reactions (e.g. '1,,,,,0.d0') kept the old
            # preprocessor happy in cooling-only/custom-ODE projects;
            # KROME-X allows networks with no reactions, so they are dropped
            self.items.append({"t": "raw", "commented": False,
                               "text": "# [legacy dummy reaction dropped: "
                                       "KROME-X allows empty networks]"})
            return
        if not rea or not pro or not rate:
            self.items.append({"t": "raw", "commented": False,
                               "text": ("#" if commented else
                                        "# [UNCONVERTED] ") + s})
            if not commented:
                self.warnings.append("line %d: could not assemble reaction — "
                                     "kept as comment" % lineno)
            return

        reaction = "%s -> %s" % (" + ".join(rea), " + ".join(pro))
        for sp in rea + pro:
            if re.match(r"^FK\d+$", sp):  # legacy hardcoded fake species
                self.fake_species.add(sp)
        register = not commented
        checkable = False

        # photo reactions marked by rate token instead of a @photo block:
        # 'krome_kph_auto', 'krome_kph_auto=<xsec(E)>', or a bare expression
        # in energy_eV; T-limit columns then hold photon-energy bounds in eV
        is_kph = rate.lower().startswith("krome_kph_auto")
        if not st.photo and (is_kph or re.search(r"\benergy_eV\b", rate)):
            def energy_bound(v, sentinel):
                try:
                    x = float(d2e(v))
                    return "" if (sentinel == "lo" and x == 0.0) or \
                                 (sentinel == "hi" and x >= 1e19) else v
                except ValueError:
                    return v
            lo, hi = energy_bound(lo, "lo"), energy_bound(hi, "hi")
            if is_kph:
                rate = "auto" if "=" not in rate else rate.split("=", 1)[1]
            st.photo, photo_restore = True, True
        else:
            photo_restore = False

        if st.photo:
            if rate.lower() == "auto":
                xsec = "auto"
            elif rate.lower().startswith("@xsecfile="):
                # legacy replaces the token by the photo-bin rate and keeps
                # the rest: "@xsecFile=SWRI * f" is the SWRI rate times f,
                # i.e. a shielding factor
                m = re.match(r"@xsecfile=\s*(\w+)\s*(.*)$", rate, re.I)
                xsec = '"%s"' % m.group(1).lower()
                # end-of-line comments are stripped as legacy does
                rest = m.group(2).split("!")[0].split("#")[0].strip()
                if rest.startswith("*"):
                    shield = self.convert_expr(rest[1:].strip(), register)
                elif rest:
                    shield = None
                    self.warnings.append(
                        "line %d: '%s' after the cross-section file is not "
                        "a factor: ignored" % (lineno, rest))
                else:
                    shield = None
            else:
                xsec = self.convert_expr(rate, register)
            args = ["xsec=%s" % xsec]
            if not rate.lower().startswith("@xsecfile="):
                shield = None
            # the 'shield' column and a '* f' after the file multiply
            sh = [x for x in (shield, fshield) if x]
            if sh:
                args.append("shield=%s" % "*".join(
                    "(%s)" % x if len(sh) > 1 else x for x in sh))
            if lo:
                args.append("Emin=%s" % lo)
            if hi:
                args.append("Emax=%s" % hi)
            new_rate, rng = "photo(%s)" % ", ".join(args), "[]"
            if xsec == "auto":
                self.photo_auto_eqs.add(reaction)
            if photo_restore:
                st.photo = False
        elif st.xray and rate.lower() == "auto":
            if reaction in self.photo_auto_eqs:
                # in KROME-X the XRAY block IS photochemistry: a reaction already
                # present as photo(xsec=auto) would be double-counted
                self.items.append({"t": "raw", "commented": False,
                                   "text": "# [legacy XRAY duplicate of "
                                           "photo reaction dropped: %s]"
                                           % reaction})
                return
            new_rate, rng = "photo(xsec=auto)", "[]"
            self.photo_auto_eqs.add(reaction)
        elif st.surface and rate.lower() == "auto":
            new_rate, rng = "surface(auto)", "[]"
        elif re.match(r"^\s*dust_adsorption_rate\s*\(", rate):
            # dust_adsorption_rate(ndust, imsqrt, stick, a^2, sqrt(Tgas)) =
            # ndust pi a^2 v_th stick: the KROME-X adsorption on the network's
            # single grain population (surface_grains, written by render)
            # with its sticking coefficient
            a = smart_split(rate[rate.index("(") + 1:rate.rindex(")")])
            grains = (self.convert_expr(a[0], register),
                      self.convert_expr(a[3], register))
            if self.grains is None:
                self.grains = grains
            elif self.grains != grains:
                self.warnings.append("dust_adsorption_rate with different "
                                     "grain arguments %s vs %s: the first "
                                     "kept" % (grains, self.grains))
            new_rate = "surface(stick=%s)" % self.convert_expr(a[2],
                                                               register)
            rng = "[]"
            self.semenov_mapped += 1
        elif re.match(r"^\s*krate_evaporation_total\s*\(", rate):
            # legacy Semenov grain functions: the KROME-X surface model
            # (codegen --binding-energies=SEMENOV, --cr-desorption=<CR param>)
            new_rate, rng = "surface(auto)", "[]"
            self.semenov_mapped += 1
        elif re.match(r"^\s*krate_2bodySi\s*\(", rate):
            # krate_2bodySi(alpha, Ea, ...): Langmuir-Hinshelwood with the
            # barrier Ea [K] and the branching alpha; for identical reactants
            # Semenov's alpha holds the 1/2 that KROME-X applies itself
            a = rate[rate.index("(") + 1:].split(",")
            alpha, Ea = float(d2e(a[0])), float(d2e(a[1]))
            if len(rea) == 2 and rea[0] == rea[1]:
                alpha *= 2.0
            new_rate = "surface(Ea=%r, branch=%r)" % (Ea, alpha)
            rng = "[]"
            self.semenov_mapped += 1
        elif rate.lower() == "auto":
            new_rate = "auto"
            rng = "[]" if (not lo and not hi) else \
                  "[%s, %s]" % (lo or "0", hi or "inf")
        else:
            check_old = normalize_expr(rate)
            new_rate = self.convert_expr(rate, register)
            checkable = True
            rng = "[]" if (not lo and not hi) else \
                  "[%s, %s]" % (lo or "0", hi or "inf")

        if commented:
            self.n_commented += 1
        else:
            self.n_reactions += 1
            if checkable:
                self.checks.append((lineno, reaction, check_old,
                                    new_rate, lo, hi, self.kseq + 1))
                self.kold[self.kseq + 1] = (check_old, lo, hi)
        if fshield and not new_rate.startswith("photo(") and not commented:
            # the KROME branches that read the column apply it to the
            # photo-processes (rate and photoheating) only
            self.warnings.append("line %d: 'shield' column on a non-photo "
                                 "reaction ignored (%s)" % (lineno, reaction))
        item = {"t": "rx", "commented": commented, "rx": reaction,
                "rng": rng, "rate": new_rate,
                "ann": ("  ; qeff=%s" % qeff) if qeff else ""}
        if getattr(st, "next_h2pd", False) and not commented:
            # legacy removed no photons: an effective pdiss (dissociations
            # per absorbed line photon) for the KROME-X photon bookkeeping
            item.update(rate="h2pd()", rng="[]",
                        ann="  ; pdiss=0.15  # effective (legacy: none)")
            st.next_h2pd = False
        if not commented:
            # legacy k() numbering: sequential order of active reactions
            self.kseq += 1
            item["kidx"] = self.kseq
        self.items.append(item)

    # -- legacy k() references and rate modifiers --------------------------------

    def mark_reverse(self):
        """legacy -reverse (every reaction reversed) as per-reaction <->.
        Legacy also reversed reactions whose reverse was already in the
        network, duplicating them; those are skipped here."""
        def split(txt):
            lhs, _, rhs = txt.partition(" -> ")
            return ([x.strip() for x in lhs.split(" + ")],
                    [x.strip() for x in rhs.split(" + ")])
        rxs = [i for i in self.items if i["t"] == "rx" and
               not i["commented"]]
        keys = {tuple(map(tuple, map(sorted, split(i["rx"])))) for i in rxs}
        n = 0
        for it in rxs:
            R, P = split(it["rx"])
            if (it["rate"].startswith(("photo(", "surface(")) or len(P) > 3
                    or any(x.startswith("{") for x in R)
                    or (tuple(sorted(P)), tuple(sorted(R))) in keys):
                continue
            it["rx"] = it["rx"].replace(" -> ", " <-> ", 1)
            n += 1
        self.warnings.append("--reverse: %d reactions written as <->" % n)

    def _rx_items(self):
        return {i["kidx"]: i for i in self.items
                if i["t"] == "rx" and "kidx" in i}

    def _after_hook(self, h, lines, head_text):
        """retitle hook h's head comment and insert `lines` (KROME-X directives)
        right after the hook's quoted legacy body."""
        k = next(j for j, i in enumerate(self.items)
                 if i.get("hook_head") == h)
        self.items[k]["text"] = head_text
        k += 1
        while k < len(self.items) and "hook_body" in self.items[k]:
            k += 1
        for ln in reversed(lines):
            self.items.insert(k, {"t": "raw", "commented": False,
                                  "text": ln})

    def apply_rate_modifiers(self):
        """translate the legacy @reactionModifier patterns found in the
        corpus (Jacobian-consistent, spec 6.2):
          k(:) = k(:) / C          -> network line `rate_scale = 1 / (C)`
          if (cond) k(a:b) = 0     -> rates a..b wrapped in
                                      where(cond, 0.0, rate)
        A block with any other statement stays an UNCONVERTED comment."""
        for h, hook in enumerate(self.hooks):
            if hook["kind"] != "reactionmodifier" or hook["commented"]:
                continue
            stmts = [re.sub(r"!.*$", "", x).strip() for x in hook["body"]]
            stmts = [x for x in stmts if x]
            actions = []
            for s in stmts:
                m = re.match(r"^k\(:\)\s*=\s*k\(:\)\s*/\s*(.+)$", s,
                             re.I)
                if m:
                    actions.append(("scale", self.convert_expr(m.group(1))))
                    continue
                m = re.match(r"^if\s*\((.+)\)\s*k\((\d+)\s*:\s*(\d+)\)"
                             r"\s*=\s*0(?:\.0*)?(?:d0)?$", s, re.I)
                if m:
                    actions.append(("zero", self.convert_expr(m.group(1)),
                                    int(m.group(2)), int(m.group(3))))
                    continue
                actions = None
                break
            if not actions:
                self.warnings.append("reactionModifier block not "
                                     "convertible: kept as comment")
                continue
            rxs = self._rx_items()
            directives = []
            for a in actions:
                if a[0] == "scale":
                    directives.append("rate_scale = 1 / (%s)" % a[1]
                                      if not re.match(r"^\w+$", a[1])
                                      else "rate_scale = 1 / %s" % a[1])
                else:
                    _, cond, lo, hi = a
                    for kk in range(lo, hi + 1):
                        i = rxs.get(kk)
                        if i is None or i["rate"] == "auto" or \
                                i["rate"].startswith(("photo(", "surface(")):
                            self.warnings.append(
                                "rate modifier: k(%d) not a plain rate, "
                                "left unmodified" % kk)
                            continue
                        i["rate"] = "where(%s, 0.0, %s)" % (cond, i["rate"])
            # the block is now expressed in KROME-X terms: mark it converted
            what = " and ".join(
                x for x, used in (
                    ("rate_scale", directives),
                    ("where() in the rates below",
                     any(a[0] == "zero" for a in actions))) if used)
            self._after_hook(h, directives,
                             "# [legacy reactionModifier converted into %s:]"
                             % what)

    def apply_ode_modifiers(self):
        """translate legacy @ODEmodifier blocks that only zero species
        derivatives into `freeze` directives (spec 6.2):
          if (cond) then; dn(idx_X) = 0; ...; end if -> freeze X, ... where cond
          if (cond) dn(idx_X) = 0                   -> freeze X where cond
          dn(idx_X) = 0                             -> freeze X
        A block with any other statement stays an UNCONVERTED comment."""
        zero = r"dn\(\s*idx_(\w+)\s*\)\s*=\s*0(?:\.0*)?(?:d0)?"
        for h, hook in enumerate(self.hooks):
            if hook["kind"] != "odemodifier" or hook["commented"]:
                continue
            stmts = [re.sub(r"!.*$", "", x).strip() for x in hook["body"]]
            stmts = [x for x in stmts if x]
            out, cond, ok = [], None, True
            for s in stmts:
                m = re.match(r"^if\s*\((.+)\)\s*then$", s, re.I)
                if m and cond is None:
                    cond = (self.convert_expr(m.group(1)), [])
                    continue
                if re.match(r"^end\s*if$", s, re.I) and cond is not None:
                    if cond[1]:
                        out.append(cond)
                    cond = None
                    continue
                m = re.match(r"^if\s*\((.+)\)\s*%s$" % zero, s, re.I)
                if m and cond is None:
                    out.append((self.convert_expr(m.group(1)),
                                [decode_idx(m.group(2))]))
                    continue
                m = re.match(r"^%s$" % zero, s, re.I)
                if m:
                    sp = decode_idx(m.group(1))
                    if cond is not None:
                        cond[1].append(sp)
                    else:
                        out.append((None, [sp]))
                    continue
                ok = False
                break
            if not ok or cond is not None or not out:
                self.warnings.append("ODEmodifier block not convertible: "
                                     "kept as comment")
                continue
            lines = ["freeze %s%s" % (", ".join(sps),
                                      " where %s" % c if c else "")
                     for c, sps in out]
            self._after_hook(h, lines, "# [legacy ODEmodifier converted "
                                       "into freeze directives:]")

    def inline_k_refs(self):
        """legacy rate references k(N) (the rate of reaction N) -> the
        referenced rate expression inlined (chains resolved); only for
        referenced reactions without a T range, else left + warned."""
        rxs = self._rx_items()
        pat = re.compile(r"(?<![A-Za-z0-9_])k\((\d+)\)")
        for _ in range(16):  # chains k(49) -> k(9)
            changed = False
            for i in rxs.values():
                def sub(m):
                    j = rxs.get(int(m.group(1)))
                    if j is None or j["rng"] != "[]" or \
                            j["rate"].startswith(("photo(", "surface(")) or \
                            j["rate"] == "auto":
                        self.warnings.append(
                            "k(%s) in %s: target not an unranged plain "
                            "rate, left as is" % (m.group(1), i["rx"]))
                        return m.group(0)
                    return "(%s)" % j["rate"]
                new = pat.sub(sub, i["rate"])
                if new != i["rate"]:
                    i["rate"] = new
                    changed = True
            if not changed:
                break
        if not any(pat.search(i["rate"]) for i in rxs.values()):
            self.externs.discard("k")

    # -- main loop ------------------------------------------------------------

    def convert(self, path: Path):
        st = BlockState()          # active stream
        st_c = BlockState()        # commented-out stream
        st.in_hook = st_c.in_hook = False

        for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
            line = raw.rstrip()
            s = line.strip()

            if st.in_hook:
                if re.match(r"@(reactionmodifier|odemodifier)_?(stop|end)",
                            s.lower()):
                    st.in_hook = False
                else:
                    self.items.append({"t": "raw", "commented": False,
                                       "text": "#   " + s,
                                       "hook_body": st.hook})
                    self.hooks[-1]["body"].append(s)
                continue

            if not s:
                self.items.append({"t": "raw", "commented": False,
                                   "text": ""})
                continue

            if s.startswith("#"):
                inner = s.lstrip("#").strip()
                if inner.lower().startswith("@ref"):
                    # refs to active reactions are conventionally commented
                    self.handle_directive(inner, st, False, lineno)
                elif inner.startswith("@"):
                    self.handle_directive(inner, st_c, True, lineno)
                elif re.match(r"^\d+\s*,", inner):
                    # a commented block may rely on the surrounding *active*
                    # @format rather than carrying its own
                    n = len(smart_split(inner))
                    use = st_c
                    if n != len(st_c.fmt) and n == len(st.fmt):
                        use = st
                    self.handle_reaction(inner, use, True, lineno)
                else:
                    self.items.append({"t": "raw", "commented": False,
                                       "text": line})
                continue

            if s.startswith("@"):
                self.handle_directive(s, st, False, lineno)
            else:
                self.handle_reaction(s, st, False, lineno)

        self.apply_rate_modifiers()
        self.inline_k_refs()
        self.apply_ode_modifiers()
        self.resolve_idx()

        defined_vars = {i["name"] for i in self.items
                        if i["t"] == "var" and not i.get("commented")}
        if self.uses_bare_T and "T" not in defined_vars:
            self.items.insert(0, {
                "t": "var", "name": "T",
                "expr": "Tgas  # legacy file used T without defining it "
                        "(legacy bug); assumed Tgas",
                "commented": False})
            self.warnings.append("bare symbol 'T' used but never defined in "
                                 "the legacy file — emitted 'var T = Tgas'")

    def resolve_idx(self):
        """the idx_<name> placeholders (decode_idx) -> species names: by the
        network's own species (legacy f90 index names, case-insensitive),
        else guess_idx (a species the reactions never name)"""
        def name(m):
            key = m.group(1).lower()
            if key == "e":
                return "e-"
            return self.idx_map.get(key) or guess_idx(m.group(1))

        def walk(v):
            if isinstance(v, str):
                return IDX_MARK.sub(name, v)
            if isinstance(v, list):
                return [walk(x) for x in v]
            if isinstance(v, tuple):
                return tuple(walk(x) for x in v)
            if isinstance(v, dict):
                return {k: walk(x) for k, x in v.items()}
            return v
        self.items = walk(self.items)
        self.checks = walk(self.checks)
        self.hooks = walk(self.hooks)

    # -- rendering ------------------------------------------------------------

    def render(self, src_name: str) -> str:
        out = ["krome_format 1",
               "# converted from legacy network '%s' by convert_legacy.py"
               % src_name, ""]
        if self.semenov_mapped:
            # legacy Semenov grain functions became surface(...) reactions
            out += ["# %d legacy Semenov grain-surface rates are KROME-X surface "
                    "reactions:" % self.semenov_mapped,
                    "# generate with surface.binding_energies = \"SEMENOV\" and "
                    "surface.cr_desorption = the CR-rate parameter"]
            # legacy -useSemenov: one grain size, user_xdust grains of
            # radius^2 user_gsize2 (sites, adsorption, ice coverage)
            g = self.grains
            if g is None and {"xdust", "gsize2"} <= set(self.params):
                g = ("xdust", "gsize2")
            if g is not None:
                out += ["# the legacy Semenov framework's grains: one size",
                        "surface_grains single(n = %s, a2 = %s)" % g]
            eb = set()
            f = paths.ROOT / "data" / \
                "Ebare_ice_Semenov.dat"
            if f.exists():  # names compared case-insensitively (legacy)
                eb = {ln.split()[0].upper() for ln in f.read_text()
                      .splitlines() if ln.strip() and not ln.startswith("#")}
            surf = set()
            for i in self.items:
                if i["t"] == "rx" and not i["commented"] and \
                        i["rate"].startswith("surface("):
                    for tok in re.split(r"\s*(?:->|\+)\s*", i["rx"]):
                        t = tok.strip()
                        if t.lower().endswith("_dust"):
                            surf.add(t[:-5])
            miss = sorted(x for x in surf if x.upper() not in eb)
            if miss:
                out += ["# no binding energy in Ebare_ice_Semenov.dat for "
                        "these ices: inert as in legacy", "# (1e99 K: they "
                        "never desorb or diffuse)"]
                out += ["ebind %s = 1e99, 1e99" % x for x in miss]
            out.append("")
        var_names = {i["name"] for i in self.items
                     if i["t"] == "var" and not i.get("commented")}
        params = {k: v for k, v in self.params.items() if k not in var_names}
        if params:
            for new, old in sorted(params.items()):
                out.append("param %-8s = 0.0    # runtime parameter "
                           "(legacy @common:%s); set default!" % (new, old))
            out.append("")
        if self.fake_species:
            out.append("species %s  ; opaque   # legacy fake species"
                       % ", ".join(sorted(self.fake_species)))
            out.append("")
        if self.used_shortcuts:
            out.append("# legacy shortcut symbols, now ordinary vars:")
            for name in SHORTCUT_ORDER:
                if name in self.used_shortcuts:
                    out.append("var %-8s = %s" % (name, SHORTCUTS[name]))
            out.append("")
        if self.externs:
            out.append("# functions provided by runtime physics modules:")
            for name in sorted(self.externs):
                out.append("extern %s" % name)
            out.append("")

        rxs = [i for i in self.items if i["t"] == "rx"]
        w_rx = min(max((len(i["rx"]) for i in rxs), default=0), 46)
        w_rng = min(max((len(i["rng"]) for i in rxs), default=0), 16)
        for it in self.items:
            pre = "# " if it.get("commented") else ""
            if it["t"] == "raw":
                out.append(it["text"])
            elif it["t"] == "var":
                out.append("%svar %s = %s" % (pre, it["name"], it["expr"]))
            elif it["t"] == "ref":
                out.append("%sref %s" % (pre, it["text"]))
            elif it["t"] == "species":
                out.append("%sspecies %s" % (pre, ", ".join(it["names"])))
            elif it["t"] == "table":
                out.append('%stable %s = "%s"' % (pre, it["name"], it["file"]))
            elif it["t"] == "cool":
                out.append("%s%s %s = %s"
                           % (pre, it["kind"], it["name"], it["expr"]))
            else:
                out.append("%s%-*s  %-*s  %s%s"
                           % (pre, w_rx, it["rx"], w_rng, it["rng"],
                              it["rate"], it["ann"]))
        return "\n".join(out) + "\n"

    # -- numerical cross-check ------------------------------------------------

    def run_check(self, quiet=False):
        var_items = [i for i in self.items
                     if i["t"] == "var" and not i.get("commented")]
        env_funcs = {"exp": math.exp, "dexp": math.exp,
                     "log": math.log, "log10": math.log10,
                     "sqrt": math.sqrt, "abs": abs, "min": min, "max": max,
                     "tanh": math.tanh}

        def base_env(T):
            Te = T * KB_EV
            e = dict(env_funcs)
            e.update(Tgas=T, Tgas_eV=Te, Te=Te, lnTe=math.log(Te),
                     invTe=1.0 / Te, invT=1.0 / T, invTgas=1.0 / T,
                     logT=math.log10(T), lnT=math.log(T),
                     sqrTgas=math.sqrt(T), invsqrT=1.0 / math.sqrt(T),
                     T32=T / 300.0, t3=T / 300.0, t4=T * 1e-4,
                     invT32=300.0 / T, sqrT32=math.sqrt(T / 300.0),
                     invsqrT32=math.sqrt(300.0 / T),
                     Tgas2=T * T, Tgas3=T ** 3, Tgas4=T ** 4,
                     T0=288.0, T02=288.0 ** 2, T03=288.0 ** 3,
                     T04=288.0 ** 4, T0inv=1.0 / 288.0,
                     Tdust=max(T, 10.0), ntot=1.0e4,
                     user_Tdust=max(T, 10.0))
            for p in self.params:
                e[p] = 1.3e-17
                e["user_" + p] = 1.3e-17
            return e

        def with_vars(env):
            failed = set()
            for it in var_items:
                try:
                    env[it["name"]] = eval(lazy_where(it["expr"]),
                                           {"__builtins__": {}}, env)
                except Exception:
                    failed.add(it["name"])
            return failed

        def uses(expr, names):
            return any(re.search(r"\b%s\b" % re.escape(x), expr)
                       for x in names)

        def legacy_k(T):
            """k(N) of the legacy expressions: reaction N's legacy rate
            at the same T (zero outside its range, as legacy)."""
            def k(N):
                old, lo, hi = self.kold[int(N)]
                if (lo and T < float(lo)) or (hi and T > float(hi)):
                    return 0.0
                e = base_env(T)
                with_vars(e)
                e["k"] = k
                return eval(lazy_where(old), {"__builtins__": {}}, e)
            return k

        final = {kk: i["rate"] for kk, i in self._rx_items().items()}
        n_ok = n_fail = n_skip = 0
        for lineno, rx, old, new, lo, hi, kidx in self.checks:
            # the final KROME-X rate (after modifier folding and k() inlining)
            new = final.get(kidx, new)
            if "n(" in old or "n(" in new or uses(new, self.externs):
                n_skip += 1
                continue
            old, new = lazy_where(old), lazy_where(new)
            try:
                tmin = float(lo) if lo else 10.0
                tmax = float(hi) if hi else 1.0e8
            except ValueError:
                n_skip += 1
                continue
            tmin, tmax = max(tmin, 2.7), min(tmax, 1.0e9)
            if tmax <= tmin:
                n_skip += 1
                continue
            pts = [10 ** (math.log10(tmin * 1.0001)
                          + k * (math.log10(tmax * 0.9999)
                                 - math.log10(tmin * 1.0001)) / 15.0)
                   for k in range(16)]
            bad = None
            for T in pts:
                e_old = base_env(T)
                e_new = base_env(T)
                failed = with_vars(e_old)
                with_vars(e_new)
                e_old["k"] = legacy_k(T)
                if uses(old, failed) or uses(new, failed):
                    bad = "SKIP"
                    break
                try:
                    v_old = eval(old, {"__builtins__": {}}, e_old)
                    v_new = eval(new, {"__builtins__": {}}, e_new)
                except Exception as exc:
                    bad = "eval failed at T=%.3e: %s" % (T, exc)
                    break
                denom = max(abs(v_old), abs(v_new), 1e-300)
                if abs(v_old - v_new) / denom > 1e-10:
                    bad = ("mismatch at T=%.3e: old=%.6e new=%.6e"
                           % (T, v_old, v_new))
                    break
            if bad == "SKIP":
                n_skip += 1
            elif bad:
                n_fail += 1
                print("CHECK FAIL line %d  %s\n   %s" % (lineno, rx, bad))
            else:
                n_ok += 1
        if not quiet or n_fail:
            print("check: %d rates verified, %d failed, %d skipped "
                  "(photo/xray/surface/auto/extern/density-dep)"
                  % (n_ok, n_fail, n_skip))
        return n_fail == 0, n_ok, n_fail, n_skip


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("legacy", type=Path)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--reverse", action="store_true",
                    help="legacy -reverse: write every eligible reaction as "
                         "A <-> B (reverse from detailed balance); photo "
                         "and curly-reactant reactions, reactions with more "
                         "than 3 products, and those whose reverse is "
                         "already in the network are left as they are")
    args = ap.parse_args()

    c = Converter()
    c.convert(args.legacy)
    if args.reverse:
        c.mark_reverse()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(c.render(args.legacy.name))

    msg = "converted %d reactions (+%d commented-out) -> %s" \
          % (c.n_reactions, c.n_commented, args.output)
    if not args.quiet:
        print(msg)
        if c.params:
            print("params:  " + ", ".join(sorted(c.params)))
        if c.externs:
            print("externs: " + ", ".join(sorted(c.externs)))
    for w in c.warnings:
        print("WARNING: %s: %s" % (args.legacy.name, w))

    ok = True
    if args.check:
        ok, n_ok, n_fail, n_skip = c.run_check(quiet=args.quiet)
        if args.quiet:
            print("%-28s %4d rx (+%3d commented)  check: %d ok / %d fail / "
                  "%d skip" % (args.legacy.name, c.n_reactions, c.n_commented,
                               n_ok, n_fail, n_skip))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
