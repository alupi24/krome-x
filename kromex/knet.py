#!/usr/bin/env python3
"""Parser for krome_format 1 network files (.knet) -> validated NetworkIR.

The IR is a JSON-serializable structure holding the species table (with
composition, charge, mass), params, vars (as expression ASTs), tables,
externs, custom cooling/heating terms, and the reaction list with typed
rates (expression / auto / photo / xray / surface), temperature ranges,
piecewise-group and additive-channel structure.

CLI:
    python3 knet.py file1.knet [file2.knet ...] [-o ir_dir] [-v]

AST node forms (JSON dicts; codegen.Codegen.diff and _cx walk them):
    {"n": 1.3e-32}                      number
    {"s": "Tgas", "kind": KIND}         symbol; KIND: "state" (Tgas, E_eV,
                                        ...), "const", "param", "var";
                                        codegen adds "ncol" (the column
                                        argument of a column-form shield)
                                        and "unclamped": True (the actual
                                        gas T, not the rate's clamped one)
    {"sp": "H2"}                        species literal (n() and extern args)
    {"sps": ["H", "H2"]}                species list (extern args)
    {"str": "swri"}                     string literal
    {"u": "-", "x": AST}                unary
    {"b": "+", "l": AST, "r": AST}      binary (+ - * / **)
    {"c": ">", "l": AST, "r": AST}      comparison (only as where() arg 1)
    {"f": "exp", "a": [AST, ...], "kind": FKIND}
                                        call; FKIND: "builtin", "table",
                                        "extern", "density" (n(species));
                                        codegen adds "tfunc" (a generated
                                        function of Tgas and constants)
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import constants as _constants  # noqa: E402

STATE_SYMBOLS = {"Tgas", "Tgas_eV", "Tdust", "ntot", "Hnuclei", "J_LW"}
# J_LW: mean intensity of the photon bins in the Lyman-Werner band
# (11.2-13.6 eV) [eV s^-1 cm^-2 Hz^-1 sr^-1]; a rate linear in it removes
# rate/pdiss photons from those bins (network-format.md 5.4)
LW_BAND = (11.2, 13.6)
# externs linear in their first argument (a J_LW rate may pass J_LW
# through them)
LINEAR_EXTERNS = {"H2_solomonLW"}
# runtime parameters every network has without declaring them
# (network-format.md 4.3)
IMPLICIT_PARAMS = {"G0": 0.0, "d2g": 1.0e-2, "Av": 0.0,
                   "gamma_ad": 5.0 / 3.0, "Tcmb": 2.725, "tff": 0.0,
                   "zredshift": 0.0, "hubble0": 0.0, "omega0": 0.3153,
                   "visc_nu": 0.0, "visc_omega": 0.0, "opratio": 3.0,
                   "total_Z": 1.0, "Tfloor": 2.725}
PHOTO_STATE = {"E_eV"}  # valid only inside photo(xsec=<expr>)
BUILTIN_FUNCS = {"exp": 1, "log": 1, "log10": 1, "sqrt": 1, "abs": 1,
                 "min": 2, "max": 2, "tanh": 1, "where": 3}
CONSTRUCTORS = {"photo", "surface"}

MASS_AMU = {
    "H": 1.008, "D": 2.014, "He": 4.0026, "Li": 6.94, "Be": 9.012,
    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098,
    "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38,
}
# the physical constants network expressions can name (kromex/constants.py)
CONSTANTS = {k: _constants.CONSTANTS[k][0] for k in (
    "Rgas_kJ", "boltzmann_eV", "boltzmann_erg", "N_avogadro", "pi",
    "planck_eV", "eV_to_erg")}
ELECTRON_MASS_AMU = 5.48579909e-4
ELEMENTS_BY_LEN = sorted(MASS_AMU, key=len, reverse=True)


class ParseError(Exception):
    def __init__(self, msg, lineno=None):
        super().__init__("line %s: %s" % (lineno or "?", msg))
        self.lineno = lineno


# ----------------------------------------------------------------------------
# species
# ----------------------------------------------------------------------------

SPECIES_RE = re.compile(
    r"^(?P<curly>\{)?(?P<iso>\d+)?(?P<core>[A-Za-z][A-Za-z0-9]*)"
    r"(?P<labels>(?:_[A-Za-z0-9]+)*)(?P<charge>\++|-+)?(?P<curly2>\})?$")


def parse_species(tok, lineno=None):
    """-> dict(name, composition, charge, mass_amu, phase, label, iso,
    kinetic, opaque). Curly braces mark a non-kinetic reactant. A core that
    is not an element composition yields an opaque species (composition
    None, zero mass) rather than an error."""
    tok = tok.strip()
    if tok == "e-":
        return {"name": "e-", "composition": {}, "charge": -1,
                "mass_amu": ELECTRON_MASS_AMU, "phase": "gas", "label": "",
                "iso": 0, "kinetic": True, "opaque": False}
    m = SPECIES_RE.match(tok)
    if not m or bool(m.group("curly")) != bool(m.group("curly2")):
        raise ParseError("malformed species %r" % tok, lineno)
    labels = [x for x in (m.group("labels") or "").split("_") if x]
    # surface phase: a trailing _ice / _dust label (any case: legacy
    # Semenov networks write X_DUST); _c_dust marks a chemisorbed species
    phase = (labels.pop().lower() if labels and labels[-1].lower()
             in ("ice", "dust") else "gas")
    core, comp, i, opaque = m.group("core"), {}, 0, False
    while i < len(core):
        for el in ELEMENTS_BY_LEN:
            if core.startswith(el, i):
                i += len(el)
                cnt = ""
                while i < len(core) and core[i].isdigit():
                    cnt += core[i]
                    i += 1
                comp[el] = comp.get(el, 0) + int(cnt or 1)
                break
        else:
            opaque, comp = True, None
            break
    iso = int(m.group("iso") or 0)
    charge = m.group("charge") or ""
    chg = len(charge) * (1 if charge.startswith("+") else -1) if charge else 0
    if opaque:
        mass = 0.0
    else:
        mass = iso * 1.0 if iso else sum(MASS_AMU[e] * n
                                         for e, n in comp.items())
    name = ((m.group("iso") or "") + core + (m.group("labels") or "")
            + charge)
    return {"name": name, "composition": comp, "charge": chg,
            "mass_amu": mass - chg * ELECTRON_MASS_AMU,
            "phase": phase, "label": "_".join(labels), "iso": iso,
            "kinetic": not m.group("curly"), "opaque": opaque}


# ----------------------------------------------------------------------------
# expression lexer + parser
# ----------------------------------------------------------------------------

TOKEN_RE = re.compile(r"""
    (?P<isotope>(?<=\bn\()\d+[A-Z][A-Za-z0-9_]*[+\-]*(?=\s*\)))  # n(60Fe)
  | (?P<num>(\d+\.?\d*|\.\d+)([eE][+\-]?\d+)?)
  | (?P<species>[A-Za-z_]\w*[+\-]+(?=\s*[,)\]]))  # charged species literal arg
  | (?P<name>[A-Za-z_]\w*)
  | (?P<str>"[^"]*")
  | (?P<op>\*\*|<=|>=|==|!=|[+\-*/(),=<>\[\]])
  | (?P<ws>\s+)
""", re.X)


def lex(text, lineno):
    toks, pos = [], 0
    while pos < len(text):
        m = TOKEN_RE.match(text, pos)
        if not m:
            raise ParseError("unexpected character %r in expression %r"
                             % (text[pos], text), lineno)
        pos = m.end()
        if m.lastgroup == "ws":
            continue
        # an isotope in n() (mass number first: n(60Fe)) is a species
        kind = "species" if m.lastgroup == "isotope" else m.lastgroup
        toks.append((kind, m.group(0)))
    toks.append(("end", ""))
    return toks


class ExprParser:
    """Recursive descent over the §4 grammar; emits JSON AST nodes."""

    def __init__(self, toks, lineno):
        self.t, self.i, self.lineno = toks, 0, lineno

    def peek(self):
        return self.t[self.i]

    def next(self):
        tok = self.t[self.i]
        self.i += 1
        return tok

    def expect(self, val):
        k, v = self.next()
        if v != val:
            raise ParseError("expected %r, got %r" % (val, v), self.lineno)

    def parse(self):
        node = self.cmp_()
        if self.peek()[0] != "end":
            raise ParseError("trailing input at %r" % self.peek()[1],
                             self.lineno)
        return node

    def cmp_(self):
        node = self.sum_()
        if self.peek()[1] in ("<", "<=", ">", ">=", "==", "!="):
            op = self.next()[1]
            node = {"c": op, "l": node, "r": self.sum_()}
        return node

    def sum_(self):
        node = self.term()
        while self.peek()[1] in ("+", "-"):
            op = self.next()[1]
            node = {"b": op, "l": node, "r": self.term()}
        return node

    def term(self):
        node = self.factor()
        while self.peek()[1] in ("*", "/"):
            op = self.next()[1]
            node = {"b": op, "l": node, "r": self.factor()}
        return node

    def factor(self):
        if self.peek()[1] in ("+", "-"):
            op = self.next()[1]
            x = self.factor()
            return x if op == "+" else {"u": "-", "x": x}
        return self.power()

    def power(self):
        node = self.atom()
        if self.peek()[1] == "**":
            self.next()
            node = {"b": "**", "l": node, "r": self.factor()}
        return node

    def atom(self):
        kind, val = self.next()
        if kind == "num":
            return {"n": float(val)}
        if kind == "str":
            return {"str": val[1:-1]}
        if kind == "species":
            return {"sp": val}
        if kind == "name":
            if self.peek()[1] == "(":
                self.next()
                args, kwargs = [], {}
                if self.peek()[1] != ")":
                    while True:
                        if (self.peek()[0] == "name"
                                and self.t[self.i + 1][1] == "="):
                            key = self.next()[1]
                            self.next()
                            kwargs[key] = self.cmp_()
                        else:
                            args.append(self.cmp_())
                        if self.peek()[1] != ",":
                            break
                        self.next()
                self.expect(")")
                node = {"f": val, "a": args}
                if kwargs:
                    node["kw"] = kwargs
                return node
            return {"s": val}
        if val == "(":
            node = self.cmp_()
            self.expect(")")
            return node
        if val == "[":  # species list (extern arguments only)
            sps = []
            while True:
                k, v = self.next()
                if k not in ("name", "species"):
                    raise ParseError("expected species in list, got %r" % v,
                                     self.lineno)
                sps.append(v)
                k, v = self.next()
                if v == "]":
                    break
                if v != ",":
                    raise ParseError("expected ',' or ']' in species list",
                                     self.lineno)
            return {"sps": sps}
        raise ParseError("unexpected token %r" % val, self.lineno)


def parse_expr(text, lineno):
    return ExprParser(lex(text, lineno), lineno).parse()


def walk(node):
    yield node
    for key in ("l", "r", "x"):
        if key in node:
            yield from walk(node[key])
    for a in node.get("a", []):
        yield from walk(a)
    for a in node.get("kw", {}).values():
        yield from walk(a)


# ----------------------------------------------------------------------------
# file-level parser
# ----------------------------------------------------------------------------

def strip_comment(line):
    out, inq = [], False
    for ch in line:
        if ch == '"':
            inq = not inq
        if ch == "#" and not inq:
            break
        out.append(ch)
    return "".join(out).rstrip()


def logical_lines(path):
    """(lineno, text) with '\\' continuations joined, comments stripped."""
    buf, start = "", 0
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        s = strip_comment(raw)
        if buf:
            s = buf + " " + s.strip()
        if s.rstrip().endswith("\\"):
            buf, start = s.rstrip()[:-1], (start or lineno)
            continue
        if s.strip():
            yield (start or lineno), s.strip()
        buf, start = "", 0
    if buf.strip():
        yield start, buf.strip()


class Network:
    def __init__(self, name):
        self.name = name
        self.species = {}          # name -> species dict
        self.params = {}           # name -> default
        self.vars = {}             # name -> AST
        self.var_order = []
        self.tables = {}           # name -> file
        self.table_log = set()     # tables interpolated in log10 axes
        self.externs = set()
        self.cooling = []
        self.heating = []
        self.reactions = []
        # AST multiplying every rate (network-format.md 6.2)
        self.rate_scale = None
        self.freeze = []           # [{species, cond (AST or None), line}]
        self.ebind = {}            # gas name -> (E_bare, E_ice) [K] (5.3)
        # grain population of the surface chemistry (5.3): None = the
        # effective MRN population; {"kind": "single", "n": ast, "a2": ast,
        # "line": k} = one grain size (legacy Semenov framework)
        self.surface_grains = None
        self.errors = []
        self.notices = []

    def add_species(self, tok, lineno):
        sp = parse_species(tok, lineno)
        key = sp["name"]
        if key not in self.species:
            self.species[key] = {k: sp[k] for k in
                                 ("name", "composition", "charge", "mass_amu",
                                  "phase", "label", "iso", "opaque")}
            if sp["opaque"]:
                self.notices.append("opaque species %r (unknown composition):"
                                    " conservation checks skipped" % key)
        return sp


RANGE_RE = re.compile(r"^\[\s*(?:([^\s,\]]+)\s*,\s*([^\s,\]]*)\s*)?\]$")


def parse_range(tok, lineno):
    m = RANGE_RE.match(tok)
    if not m:
        raise ParseError("malformed range %r" % tok, lineno)
    if m.group(1) is None:
        return None
    lo = float(m.group(1)) if m.group(1) not in ("", "0") else 0.0
    hi = float(m.group(2)) if m.group(2) not in ("", "inf") else math.inf
    if hi < lo:
        raise ParseError("empty range %r" % tok, lineno)
    return [lo, hi]


def split_reaction_line(text, lineno):
    """-> (lhs, rhs, range_tok, rate_text, annotations dict, reversible);
    `A <-> B` declares the reverse reaction from detailed balance (spec
    3.2)"""
    reversible = "<->" in text.split("[", 1)[0]
    lhs, _, rest = text.partition("<->" if reversible else "->")
    if not _:
        raise ParseError("missing '->'", lineno)
    b0 = rest.find("[")
    if b0 < 0:
        raise ParseError("missing [Tmin, Tmax] range field", lineno)
    b1 = rest.find("]", b0)
    if b1 < 0:
        raise ParseError("unterminated range field", lineno)
    rhs, rng = rest[:b0].strip(), rest[b0:b1 + 1]
    tail = rest[b1 + 1:].strip()
    parts, depth, cur = [], 0, ""
    for ch in tail:
        depth += ch in "("
        depth -= ch in ")"
        if ch == ";" and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    rate = parts[0].strip()
    ann = {}
    for p in parts[1:]:
        k, _, v = p.partition("=")
        if not _:
            raise ParseError("malformed annotation %r" % p.strip(), lineno)
        ann[k.strip()] = v.strip()
    if not rate:
        raise ParseError("missing rate", lineno)
    ann["__reversible__"] = reversible
    return lhs.strip(), rhs, rng, rate, ann


def parse_stoich(text, net, lineno, allow_curly):
    out = []
    for tok in re.split(r"\s\+\s", " " + text + " "):
        tok = tok.strip()
        if not tok:
            raise ParseError("empty species in %r" % text, lineno)
        sp = net.add_species(tok, lineno)
        if not sp["kinetic"] and not allow_curly:
            raise ParseError("curly (non-kinetic) species %r only allowed "
                             "among reactants" % tok, lineno)
        out.append({"sp": sp["name"], "kinetic": sp["kinetic"]})
    return out


def parse_rate(rate, ann, net, lineno):
    if rate == "auto":
        return {"kind": "auto"}
    if re.match(r"^h2pd\s*\(\s*\)$", rate):
        # H2 photodissociation from the H2pdB.dat line data on the photon
        # bins (legacy @next_H2photodissociation, kpd_H2;
        # network-format.md 5.4)
        return {"kind": "photo", "xsec": {"kind": "h2pd"}}
    call = re.match(r"^(photo|surface)\s*\((.*)\)$", rate, flags=re.S)
    if call:
        kind, inner = call.group(1), call.group(2).strip()
        if kind == "surface":
            # surface(auto): the surface model resolves the process;
            # surface(Ea=K | prob=p [, branch=f]): a two-body reaction with
            # its barrier [K] or reaction probability given
            # (network-format.md 5.3)
            if inner == "auto":
                return {"kind": kind, "xsec": {"kind": "auto"}}
            node = parse_expr(rate, lineno)
            kw = node.get("kw", {})
            # surface(stick=expr): an adsorption with its sticking
            # coefficient given (instead of the HM79 one)
            two = set(kw) - {"stick"}
            if node["a"] or set(kw) - {"Ea", "prob", "branch", "stick"} or \
                    ("stick" in kw and two) or \
                    ("stick" not in kw and ("Ea" in kw) == ("prob" in kw)):
                raise ParseError("surface(auto), surface(stick=...) or "
                                 "surface(Ea=... | prob=...[, branch=...])",
                                 lineno)
            out = {"kind": kind, "xsec": {"kind": "auto"}}
            for k in kw:
                if k == "stick":
                    out["stick"] = kw[k]
                    continue
                if "n" not in kw[k]:
                    raise ParseError("surface() %s must be a number" % k,
                                     lineno)
                out[k] = kw[k]["n"]
            return out
        node = parse_expr(rate, lineno)  # reuse call parsing for kwargs
        kw = node.get("kw", {})
        if "xsec" not in kw:
            raise ParseError("photo() requires xsec=", lineno)
        x = kw["xsec"]
        if x == {"s": "auto"}:
            xsec = {"kind": "auto"}
        elif "str" in x:
            xsec = {"kind": "ref", "ref": x["str"]}
        else:
            xsec = {"kind": "expr", "ast": x}
        out = {"kind": "photo", "xsec": xsec}
        # shielding factor of the rate (network-format.md 5.4)
        if "shield" in kw:
            out["shield"] = kw["shield"]
        for k in ("Emin", "Emax"):
            if k in kw:
                if "n" not in kw[k]:
                    raise ParseError("photo() %s must be a number" % k, lineno)
                out[k] = kw[k]["n"]
        if set(kw) - {"xsec", "Emin", "Emax", "shield"} or node["a"]:
            raise ParseError("unexpected photo() arguments", lineno)
        return out
    return {"kind": "expr", "ast": parse_expr(rate, lineno)}


def resolve_symbols(net):
    """Classify names; check unknowns, arities, var cycles, E_eV scope,
    comparison placement."""
    def classify(node, ctx):
        for nd in walk(node):
            if "s" in nd and "kind" not in nd:
                s = nd["s"]
                if s in STATE_SYMBOLS:
                    nd["kind"] = "state"
                elif s in PHOTO_STATE:
                    if ctx != "xsec":
                        net.errors.append("E_eV used outside a photo cross "
                                          "section (in %s)" % ctx)
                    nd["kind"] = "state"
                elif s in CONSTANTS:
                    nd["kind"] = "const"
                elif s in net.params or s in IMPLICIT_PARAMS:
                    nd["kind"] = "param"
                elif s in net.vars:
                    nd["kind"] = "var"
                elif s in net.species:
                    nd.clear()
                    nd["sp"] = s
                else:
                    net.errors.append("unknown symbol %r in %s" % (s, ctx))
            elif "sp" in nd and len(nd) == 1:
                if nd["sp"] not in net.species:
                    net.errors.append("unknown species %r in %s"
                                      % (nd["sp"], ctx))
            elif "f" in nd and "kind" not in nd:
                f = nd["f"]
                if f == "n":
                    a = nd["a"]
                    if len(a) != 1 or not ("sp" in a[0] or "s" in a[0]):
                        net.errors.append("n() takes one species (in %s)"
                                          % ctx)
                    elif "s" in a[0]:
                        name = a[0]["s"]
                        a[0].clear()
                        a[0]["sp"] = name
                        if name not in net.species:
                            net.errors.append("unknown species %r in n() "
                                              "(in %s)" % (name, ctx))
                    nd["kind"] = "density"
                elif f in BUILTIN_FUNCS:
                    if len(nd["a"]) != BUILTIN_FUNCS[f]:
                        net.errors.append("%s() takes %d args (in %s)"
                                          % (f, BUILTIN_FUNCS[f], ctx))
                    nd["kind"] = "builtin"
                elif f in net.tables:
                    nd["kind"] = "table"
                elif f in net.externs:
                    nd["kind"] = "extern"
                    for a in nd["a"]:
                        if "s" in a and a["s"] in net.species:
                            name = a["s"]
                            a.clear()
                            a["sp"] = name
                else:
                    net.errors.append("unknown function %r in %s" % (f, ctx))
            if "c" in nd and ctx_allows_cmp.get(id(nd)) is not True:
                net.errors.append("comparison outside where() in %s" % ctx)

    # mark comparisons that sit directly under where() as allowed
    ctx_allows_cmp = {}

    def mark_where(node):
        for nd in walk(node):
            if nd.get("f") == "where" and nd.get("a"):
                if "c" in nd["a"][0]:
                    ctx_allows_cmp[id(nd["a"][0])] = True

    sg = net.surface_grains
    everything = ([("var %s" % k, v) for k, v in net.vars.items()]
                  + ([("rate_scale", net.rate_scale)]
                     if net.rate_scale is not None else [])
                  + ([("surface_grains %s" % k, sg[k]) for k in ("n", "a2")]
                     if sg is not None else [])
                  + [("freeze condition (line %d)" % f["line"], f["cond"])
                     for f in net.freeze if f["cond"] is not None]
                  + [("cooling %s" % c["name"], c["expr"]) for c in net.cooling]
                  + [("heating %s" % c["name"], c["expr"]) for c in net.heating])
    for rx in net.reactions:
        r = rx["rate"]
        if r["kind"] == "expr":
            everything.append(("reaction %s" % rx["eq"], r["ast"]))
        elif r["kind"] == "photo" and r["xsec"]["kind"] == "expr":
            everything.append(("xsec", r["xsec"]["ast"]))
        if r["kind"] == "photo" and r.get("shield") is not None:
            everything.append(("shield of %s" % rx["eq"], r["shield"]))
        if r["kind"] == "surface" and r.get("stick") is not None:
            everything.append(("stick of %s" % rx["eq"], r["stick"]))
        if rx.get("qeff") is not None:
            everything.append(("qeff of %s" % rx["eq"], rx["qeff"]))

    def check_sps_lists(node, ctx):
        allowed = set()
        for nd in walk(node):
            if nd.get("f") in net.externs:
                for a in nd.get("a", []):
                    if "sps" in a:
                        allowed.add(id(a))
        for nd in walk(node):
            if "sps" in nd:
                if id(nd) not in allowed:
                    net.errors.append("species list outside an extern call "
                                      "in %s" % ctx)
                for s0 in nd["sps"]:
                    if s0 not in net.species:
                        net.errors.append("unknown species %r in species "
                                          "list (in %s)" % (s0, ctx))

    for ctx, node in everything:
        mark_where(node)
        if ctx.startswith("freeze") and "c" in node:
            ctx_allows_cmp[id(node)] = True
    for ctx, node in everything:
        classify(node, "xsec" if ctx == "xsec" else ctx)
        check_sps_lists(node, ctx)

    # var dependency cycles + evaluation order
    deps = {k: {nd["s"] for nd in walk(v) if nd.get("kind") == "var"}
            for k, v in net.vars.items()}
    order, seen, stack = [], set(), set()

    def visit(k):
        if k in seen:
            return
        if k in stack:
            net.errors.append("cyclic var definition involving %r" % k)
            return
        stack.add(k)
        for d in sorted(deps.get(k, ())):
            visit(d)
        stack.discard(k)
        seen.add(k)
        order.append(k)

    for k in net.var_order:
        visit(k)
    net.var_order = [k for k in order if k in net.vars]


def check_conservation(net):
    for rx in net.reactions:
        sps = [net.species[r["sp"]] for r in rx["reactants"]] \
            + [net.species[p["sp"]] for p in rx["products"]]
        nuclear = any(s["iso"] for s in sps)
        rx["nuclear"] = nuclear
        if any(s.get("opaque") for s in sps):
            continue  # opaque species: conservation not checkable
        R = [net.species[r["sp"]] for r in rx["reactants"]]
        P = [net.species[p["sp"]] for p in rx["products"]]
        if nuclear:
            a = sum(s["iso"] or round(sum(MASS_AMU[e] * n for e, n in
                    s["composition"].items())) for s in R)
            b = sum(s["iso"] or round(sum(MASS_AMU[e] * n for e, n in
                    s["composition"].items())) for s in P)
            if a != b:
                net.errors.append("nucleon imbalance in %s (%d vs %d)"
                                  % (rx["eq"], a, b))
            continue
        for side, other in ((R, P),):
            ea = {}
            eb = {}
            for s in R:
                for e, n in s["composition"].items():
                    ea[e] = ea.get(e, 0) + n
            for s in P:
                for e, n in s["composition"].items():
                    eb[e] = eb.get(e, 0) + n
            if ea != eb:
                net.errors.append("element imbalance in %s (%s vs %s)"
                                  % (rx["eq"], ea, eb))
        if sum(s["charge"] for s in R) != sum(s["charge"] for s in P):
            net.errors.append("charge imbalance in %s" % rx["eq"])


def group_reactions(net):
    """Piecewise groups vs additive channels (§3.1); photo E-window overlap."""
    by_eq = {}
    for i, rx in enumerate(net.reactions):
        by_eq.setdefault(rx["eq"], []).append(i)
    gid = 0
    for eq, idxs in by_eq.items():
        limited = sorted((i for i in idxs if net.reactions[i]["trange"]),
                         key=lambda i: net.reactions[i]["trange"])
        unlimited = [i for i in idxs if not net.reactions[i]["trange"]]
        chains = []  # each chain: consecutive tiling pieces of one rate
        for i in limited:
            lo, hi = net.reactions[i]["trange"]
            if lo == hi:
                net.notices.append("degenerate T range [%g, %g] for %s"
                                   % (lo, hi, eq))
            for ch in chains:
                if net.reactions[ch[-1]]["trange"][1] == lo:
                    ch.append(i)
                    break
            else:
                # not a continuation: a new chain overlapping existing
                # coverage is a parallel channel; beyond all coverage it is
                # an interior gap in a piecewise rate
                if chains and all(
                        lo >= net.reactions[ch[-1]]["trange"][1]
                        for ch in chains):
                    net.errors.append("gap in piecewise T ranges for %s" % eq)
                chains.append([i])
        for ch in chains:
            gid += 1
            for i in ch:
                net.reactions[i]["group"] = gid
        nch = len(chains) + len(unlimited)
        if nch > 1:
            net.notices.append("%d additive channels for %s" % (nch, eq))
        photo = [i for i in idxs
                 if net.reactions[i]["rate"]["kind"] == "photo"]
        wins = sorted((net.reactions[i]["rate"].get("Emin", 0.0),
                       net.reactions[i]["rate"].get("Emax", math.inf))
                      for i in photo)
        for (a1, b1), (a2, b2) in zip(wins, wins[1:]):
            if a2 < b1:
                net.errors.append("overlapping photo energy windows for %s"
                                  % eq)


def _split_where(body):
    """'A, B where cond' -> ('A, B', 'cond'); no where -> (body, '')."""
    m = re.match(r"^(.*?)\s+where\s+(.+)$", body)
    return (m.group(1), m.group(2).strip()) if m else (body, "")


def jlw_degree(node, net):
    """degree of an expression in J_LW: 0 (independent), 1 (linear:
    J_LW times a J_LW-free factor), or None (anything else)"""
    if not any(nd.get("s") == "J_LW" for nd in walk(node)):
        return 0
    if node.get("s") == "J_LW":
        return 1
    if "u" in node:
        return jlw_degree(node["x"], net)
    if "b" in node:
        dl, dr = jlw_degree(node["l"], net), jlw_degree(node["r"], net)
        if None in (dl, dr):
            return None
        if node["b"] == "*":
            return dl + dr if dl + dr <= 1 else None
        if node["b"] == "/":
            return dl if dr == 0 else None
        if node["b"] in ("+", "-"):
            return dl if dl == dr else None
        return None
    if "f" in node:
        f, a = node["f"], node.get("a", [])
        if f in LINEAR_EXTERNS and a and all(
                jlw_degree(x, net) == 0 for x in a[1:]):
            return jlw_degree(a[0], net)
        if f == "where" and len(a) == 3 and jlw_degree(a[0], net) == 0:
            d = [jlw_degree(x, net) for x in a[1:]]
            zero = [x.get("n") == 0.0 for x in a[1:]]
            if d == [1, 1] or (d == [1, 0] and zero[1]) \
                    or (d == [0, 1] and zero[0]):
                return 1
        return None
    return None


def check_reverse(net):
    """`<->` reactions (network-format.md 3.2): an expression or `auto` rate
    between kinetic species only; the reverse must not also be written out."""
    def key(rx):
        return (tuple(sorted(r["sp"] for r in rx["reactants"])),
                tuple(sorted(p["sp"] for p in rx["products"])))
    keys = {key(rx) for rx in net.reactions}
    for rx in net.reactions:
        if not rx.get("reverse"):
            continue
        r = rx["rate"]
        bad = None
        if r["kind"] not in ("expr", "auto"):
            bad = "only expression and auto rates can be reversed"
        elif any(not q["kinetic"] for q in rx["reactants"]):
            bad = "curly (non-kinetic) reactants cannot be reversed"
        elif rx.get("lw") or rx.get("pdiss") is not None:
            bad = "radiative channels cannot be reversed"
        elif rx.get("qeff") is not None:
            bad = "a qeff reaction cannot be reversed"
        elif (key(rx)[1], key(rx)[0]) in keys:
            bad = ("its reverse is also written in the network (drop one, "
                   "or use ->)")
        if bad:
            net.errors.append("%s (<->): %s" % (rx["eq"], bad))


def check_radiative(net):
    """J_LW rates (network-format.md 5.4): linear in J_LW, one kinetic
    reactant (the absorber), no temperature range, pdiss given. pdiss is
    allowed on photo() and J_LW rates only; J_LW appears in reaction rates
    only."""
    for k, v in net.vars.items():
        if any(nd.get("s") == "J_LW" for nd in walk(v)):
            net.errors.append("J_LW used in var %s: write it in the rate "
                              "expression itself" % k)
    for kind in ("cooling", "heating"):
        for c in getattr(net, kind):
            if any(nd.get("s") == "J_LW" for nd in walk(c["expr"])):
                net.errors.append("J_LW used in %s %s" % (kind, c["name"]))
    for rx in net.reactions:
        r = rx["rate"]
        if r["kind"] == "photo" and r.get("shield") is not None:
            if any(nd.get("s") in ("J_LW", "E_eV")
                   for nd in walk(r["shield"])):
                net.errors.append("shield of %s: J_LW / E_eV not allowed"
                                  % rx["eq"])
        lw = r["kind"] == "expr" and any(nd.get("s") == "J_LW"
                                         for nd in walk(r["ast"]))
        rx["lw"] = lw
        if lw:
            if jlw_degree(r["ast"], net) != 1:
                net.errors.append("%s: a J_LW rate must be linear in J_LW "
                                  "(J_LW times a J_LW-free factor)"
                                  % rx["eq"])
            if sum(1 for q in rx["reactants"] if q["kinetic"]) != 1:
                net.errors.append("%s: a J_LW rate needs exactly one kinetic "
                                  "reactant (the absorber)" % rx["eq"])
            if rx["trange"]:
                net.errors.append("%s: a J_LW rate takes no temperature "
                                  "range ([])" % rx["eq"])
            if "pdiss" not in rx:
                net.errors.append("%s: a J_LW rate needs ; pdiss=<events "
                                  "per absorbed photon>" % rx["eq"])
        elif r["kind"] == "photo" and r["xsec"]["kind"] == "h2pd" and (
                sorted(q["sp"] for q in rx["reactants"]) != ["H2"]
                or sorted(q["sp"] for q in rx["products"]) != ["H", "H"]
                or "pdiss" not in rx):
            net.errors.append("%s: h2pd() is H2 -> H + H with ; pdiss= "
                              "(dissociations per absorbed photon)"
                              % rx["eq"])
        elif "pdiss" in rx and r["kind"] != "photo":
            net.errors.append("%s: pdiss is only for photo() and J_LW rates"
                              % rx["eq"])


def resolve_modifiers(net):
    """rate_scale / freeze checks (network-format.md 6.2). Both are evaluated
    with the rates, so they may use what a rate may use, except that a freeze
    condition must be constant during a solver call (parameters and
    constants only): a switch on the evolving state would make the RHS
    discontinuous."""
    for fr in net.freeze:
        for s in fr["species"]:
            if s not in net.species:
                net.errors.append("freeze (line %d): unknown species %r"
                                  % (fr["line"], s))
        c = fr["cond"]
        if c is None:
            continue
        for nd in walk(c):
            k = nd.get("kind")
            if k not in (None, "param", "const", "builtin") or "sp" in nd:
                net.errors.append(
                    "freeze condition (line %d) may use only parameters and "
                    "constants (it must stay fixed during a solver call)"
                    % fr["line"])
                break
    sg = net.surface_grains
    if sg is not None:
        for k in ("n", "a2"):
            for nd in walk(sg[k]):
                kk = nd.get("kind")
                if kk not in (None, "param", "const", "builtin") \
                        or "sp" in nd:
                    net.errors.append(
                        "surface_grains (line %d): %s may use only "
                        "parameters and constants (a cell property, fixed "
                        "during a solver call)" % (sg["line"], k))
                    break
    if net.rate_scale is not None and any(
            rx["rate"]["kind"] == "photo" for rx in net.reactions):
        net.errors.append("rate_scale with photo() reactions is not "
                          "supported (a photo rate is sigma x flux; scale "
                          "the cross section instead)")


def parse_network(path):
    net = Network(Path(path).stem)
    lines = list(logical_lines(Path(path)))
    if not lines or not re.match(r"^krome_format\s+1$", lines[0][1]):
        net.errors.append("missing 'krome_format 1' header")
        return net
    pending_ref = None
    ncool = {"cooling": 0, "heating": 0}

    for lineno, text in lines[1:]:
        try:
            m = re.match(r"^rate_scale\s*=\s*(.+)$", text)
            if m:
                if net.rate_scale is not None:
                    raise ParseError("rate_scale given twice", lineno)
                net.rate_scale = parse_expr(m.group(1).strip(), lineno)
                continue
            # ebind X = E_bare, E_ice [K]: a surface binding energy that
            # completes or overrides the --binding-energies data set
            m = re.match(r"^ebind\s+(\S+)\s*=\s*([^,]+),\s*(.+)$", text)
            if m:
                try:
                    net.ebind[m.group(1)] = (float(m.group(2)),
                                             float(m.group(3).split("#")[0]))
                except ValueError:
                    raise ParseError("ebind X = E_bare, E_ice [K]", lineno)
                continue
            # surface_grains single(n = expr, a2 = expr | a = expr): the
            # grain population of the surface chemistry (default: MRN)
            m = re.match(r"^surface_grains\s+(.+)$", text)
            if m and "->" not in text:
                if net.surface_grains is not None:
                    raise ParseError("surface_grains given twice", lineno)
                body = m.group(1).split("#")[0].strip()
                if body.lower() == "mrn":
                    continue
                node = parse_expr(body, lineno)
                kw = node.get("kw", {})
                if node.get("f") != "single" or node.get("a") or \
                        "n" not in kw or ("a2" in kw) == ("a" in kw) or \
                        set(kw) - {"n", "a", "a2"}:
                    raise ParseError("surface_grains mrn, or surface_grains "
                                     "single(n = <grains cm^-3>, a2 = "
                                     "<radius^2 cm^2> | a = <radius cm>)",
                                     lineno)
                a2 = kw["a2"] if "a2" in kw else \
                    {"b": "*", "l": kw["a"], "r": kw["a"]}
                net.surface_grains = {"kind": "single", "n": kw["n"],
                                      "a2": a2, "line": lineno}
                continue
            m = re.match(r"^freeze\s+(.+)$", text)
            if m and "->" not in text:
                sps, cond = _split_where(m.group(1))
                names = [s.strip() for s in sps.split(",") if s.strip()]
                if not names:
                    raise ParseError("freeze needs at least one species",
                                     lineno)
                net.freeze.append({
                    "species": names, "line": lineno,
                    "cond": parse_expr(cond, lineno) if cond else None})
                continue
            m = re.match(r"^(param|var|table|extern|species|ref|cooling|"
                         r"heating)\s+(.*)$", text)
            if m and "->" not in text.split("=")[0]:
                kw, body = m.group(1), m.group(2).strip()
                if kw == "param":
                    name, _, val = body.partition("=")
                    net.params[name.strip()] = float(val.strip())
                elif kw == "var":
                    name, _, val = body.partition("=")
                    name = name.strip()
                    net.vars[name] = parse_expr(val.strip(), lineno)
                    net.var_order.append(name)
                elif kw == "table":
                    # table NAME = "file" [log]
                    name, _, val = body.partition("=")
                    tm = re.match(r'^"?([^"\s]+)"?\s*(\w+)?\s*$',
                                  val.strip())
                    if not tm or tm.group(2) not in (None, "log"):
                        raise ParseError("table needs NAME = \"file\" "
                                         "[log]", lineno)
                    net.tables[name.strip()] = tm.group(1)
                    if tm.group(2) == "log":
                        net.table_log.add(name.strip())
                elif kw == "extern":
                    net.externs.add(body.split()[0])
                elif kw == "species":
                    names, _, s_ann = body.partition(";")
                    force_opaque = s_ann.split("#")[0].strip() == "opaque"
                    for tok in names.split(","):
                        sp = net.add_species(tok.strip(), lineno)
                        if force_opaque:
                            rec = net.species[sp["name"]]
                            rec.update(opaque=True, composition=None,
                                       mass_amu=0.0)
                elif kw == "ref":
                    pending_ref = body
                else:  # cooling / heating
                    name, _, val = body.partition("=")
                    ncool[kw] += 1
                    getattr(net, kw).append(
                        {"name": name.strip(),
                         "expr": parse_expr(val.strip(), lineno)})
                continue

            lhs, rhs, rng_tok, rate, ann = split_reaction_line(text, lineno)
            rx = {
                "id": len(net.reactions) + 1,
                "reactants": parse_stoich(lhs, net, lineno, allow_curly=True),
                "products": parse_stoich(rhs, net, lineno, allow_curly=False),
                "trange": parse_range(rng_tok, lineno),
                "rate": parse_rate(rate, ann, net, lineno),
                "line": lineno,
            }
            rx["eq"] = "%s -> %s" % (
                " + ".join(r["sp"] + ("" if r["kinetic"] else "{}")
                           for r in rx["reactants"]),
                " + ".join(p["sp"] for p in rx["products"]))
            rx["reverse"] = ann.pop("__reversible__")
            rx["qeff"] = (parse_expr(ann.pop("qeff"), lineno)
                          if "qeff" in ann else None)
            if "pdiss" in ann:
                try:
                    rx["pdiss"] = float(ann.pop("pdiss"))
                except ValueError:
                    raise ParseError("pdiss must be a number", lineno)
                if not 0.0 < rx["pdiss"] <= 1.0:
                    raise ParseError("pdiss must be in (0, 1]", lineno)
            if ann:
                raise ParseError("unknown annotation(s) %s" % list(ann),
                                 lineno)
            if pending_ref:
                rx["ref"] = pending_ref
                pending_ref = None
            net.reactions.append(rx)
        except ParseError as e:
            net.errors.append(str(e))

    resolve_symbols(net)
    resolve_modifiers(net)
    check_radiative(net)
    check_reverse(net)
    check_conservation(net)
    group_reactions(net)
    return net


def to_json(net):
    return {
        "format": 1,
        "name": net.name,
        "species": [net.species[k] for k in sorted(net.species)],
        "params": [{"name": k, "default": v}
                   for k, v in sorted(net.params.items())],
        "vars": [{"name": k, "expr": net.vars[k]} for k in net.var_order],
        "tables": [{"name": k, "file": v, "log": k in net.table_log}
                   for k, v in sorted(net.tables.items())],
        "externs": sorted(net.externs),
        "cooling": net.cooling,
        "heating": net.heating,
        "reactions": net.reactions,
        "rate_scale": net.rate_scale,
        "surface_grains": net.surface_grains,
        "freeze": net.freeze,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("-o", "--out", type=Path, help="directory for IR JSON")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    n_err = 0
    for f in args.files:
        net = parse_network(f)
        status = "OK " if not net.errors else "ERR"
        print("%s %-28s %4d rx  %3d species  %2d vars  %2d params  "
              "%d externs%s"
              % (status, f.name, len(net.reactions), len(net.species),
                 len(net.vars), len(net.params), len(net.externs),
                 ("  [%d notices]" % len(net.notices)) if net.notices else ""))
        for e in net.errors:
            n_err += 1
            print("   ERROR: %s" % e)
        if args.verbose:
            for n in net.notices:
                print("   note: %s" % n)
        if args.out and not net.errors:
            args.out.mkdir(parents=True, exist_ok=True)
            (args.out / (f.stem + ".json")).write_text(
                json.dumps(to_json(net), indent=1))
    sys.exit(1 if n_err else 0)


if __name__ == "__main__":
    main()
