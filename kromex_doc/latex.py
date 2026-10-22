"""LaTeX from a parsed network: the rate expressions (from their ASTs, no
sympy), the species and reactions (mhchem for the HTML pages, plain math
for the table: mhchem is too slow for thousands of reactions), and the
reaction table of a paper (longtable, amsmath)."""

import math
import re

import knet

# precedence: sum 1, product 2, unary minus 3, power 4, atom 5
STATE_TEX = {"Tgas": "T", "Tgas_eV": r"T_\mathrm{eV}", "Tdust": r"T_\mathrm{d}",
             "ntot": r"n_\mathrm{tot}", "Hnuclei": r"n_\mathrm{H}",
             "J_LW": r"J_\mathrm{LW}", "E_eV": "E"}
CONST_TEX = {"pi": r"\pi", "boltzmann_eV": r"k_\mathrm{B}",
             "boltzmann_erg": r"k_\mathrm{B}", "planck_eV": "h",
             "Rgas_kJ": "R", "N_avogadro": r"N_\mathrm{A}",
             "eV_to_erg": r"\mathrm{eV}"}
FUNC_TEX = {"exp": r"\exp", "log": r"\ln", "log10": r"\log_{10}",
            "tanh": r"\tanh", "min": r"\min", "max": r"\max"}
CMP_TEX = {">": ">", "<": "<", ">=": r"\geq", "<=": r"\leq", "==": "=",
           "!=": r"\neq"}


def num(v):
    """a number as LaTeX: plain when short, a x 10^b otherwise"""
    if v == 0:
        return "0"
    if float(v).is_integer() and abs(v) < 1e5:
        return "%d" % v
    s = "%.6g" % v
    if "e" not in s and 1e-3 <= abs(v) < 1e5:
        return s
    m, e = ("%.6e" % v).split("e")
    m = m.rstrip("0").rstrip(".")
    e = int(e)
    return r"10^{%d}" % e if m == "1" else r"%s\times10^{%d}" % (m, e)


def name_tex(s):
    """a param or var name: upright, underscores as subscripts"""
    head, _, rest = s.partition("_")
    head = r"\mathrm{%s}" % head
    return head + (r"_\mathrm{%s}" % rest.replace("_", r"\_") if rest
                   else "")


def species_parts(sp):
    """a KROME-X species name in parts: (prefix, isotope, core, label,
    charge, ice). The prefix: o-/p-/m- (ortho/para/meta), c-/l-/t-/i-/n-
    (structure); the label: what follows an underscore (an excited state
    O_1D, a site); ice: X_DUST, X_ice; the grains as gr"""
    name, ice, pre = sp, False, ""
    m = re.match(r"^(.*?)_(dust|ice)$", name, re.I)
    if m:
        name, ice = m.group(1), True
    m = re.match(r"^(.*?)([+\-]*)$", name)
    name, charge = m.group(1), m.group(2)
    m = re.match(r"^(.*?)_(ortho|para|meta)$", name)
    if m:
        name = m.group(1)
        pre = {"ortho": "o-", "para": "p-", "meta": "m-"}[m.group(2)]
    m = re.match(r"^([cltin])_(.+)$", name)
    if m:
        pre, name = pre + m.group(1) + "-", m.group(2)
    m = re.match(r"^GRAIN(0?)$", name)
    if m:
        name = "gr"
    m = re.match(r"^(\d*)([^_]*)_?(.*)$", name)
    iso, core, label = m.groups()
    if label in ("c", "l", "t", "i", "n"):   # an isomer label: C3H_l
        pre, label = pre + label + "-", ""
    return pre, iso, core, label, charge, ice


def species_tex(sp):
    """a species in mhchem: X_DUST and X_ice as X(s), ortho/para as o-/p-,
    O_1D as O(^1D), the grains as gr"""
    pre, iso, core, label, charge, ice = species_parts(sp)
    out = pre + ("^{%s}" % iso if iso else "") + core
    if label:
        m = re.match(r"^(\d+)([A-Za-z]+)$", label)
        out += ("(^{%s}%s)" % m.groups() if m else
                "(%s)" % label.replace("_", ","))
    return out + charge + ("(s)" if ice else "")


def reaction_tex(rx, net=None):
    """the reaction in mhchem (no $...$): A + B -> C + D (the HTML pages:
    MathJax)"""
    lhs = " + ".join(species_tex(r["sp"]) for r in rx["reactants"])
    rhs = " + ".join(species_tex(p["sp"]) for p in rx["products"])
    if rx["rate"]["kind"] == "photo":
        lhs += r" + $\gamma$"
    return r"\ce{%s -> %s}" % (lhs, rhs)


def species_math(sp):
    """a species in plain math mode (the LaTeX table: mhchem's \\ce is
    too slow for thousands of reactions): H2+ -> \\mathrm{H}_{2}^{+}"""
    pre, iso, core, label, charge, ice = species_parts(sp)
    out = r"\mathrm{%s}" % pre if pre else ""
    if iso:
        out += "{}^{%s}" % iso
    for tok in re.findall(r"[A-Za-z]+|\d+|.", core):
        out += "_{%s}" % tok if tok.isdigit() else r"\mathrm{%s}" % tok
    if label:
        m = re.match(r"^(\d+)([A-Za-z]+)$", label)
        out += (r"({}^{%s}\mathrm{%s})" % m.groups() if m else
                r"(\mathrm{%s})" % label.replace("_", ","))
    if charge:
        n = len(charge)
        out += "^{%s%s}" % (n if n > 1 else "", charge[0])
    if ice:
        out += r"(\mathrm{s})"
    return out


def reaction_math(rx):
    """the reaction in plain math mode, $A + B \\rightarrow C$"""
    lhs = " + ".join(species_math(r["sp"]) for r in rx["reactants"])
    rhs = " + ".join(species_math(p["sp"]) for p in rx["products"])
    if rx["rate"]["kind"] == "photo":
        lhs += r" + \gamma"
    return r"$%s \rightarrow %s$" % (lhs, rhs)


class Tex:
    """an expression as LaTeX; the vars of inline (a short definition, a
    few nodes) are written in place, the others by name"""

    def __init__(self, net, inline=None):
        self.net = net
        self.inline = (set(inline) if inline is not None else
                       {v for v, a in net.vars.items() if size(a) <= 3})

    def __call__(self, nd):
        return self.go(nd)[0]

    def go(self, nd):
        """(LaTeX, precedence)"""
        if "n" in nd:
            t = num(nd["n"])
            return t, (5 if nd["n"] >= 0 and r"\times" not in t else 2)
        if "s" in nd:
            k, s = nd.get("kind"), nd["s"]
            if k == "var" and s in self.inline:
                return self.go(self.net.vars[s])
            if k == "state" or s in STATE_TEX:
                return STATE_TEX.get(s, name_tex(s)), 5
            if k == "const":
                return CONST_TEX.get(s, name_tex(s)), 5
            return name_tex(s), 5
        if "sp" in nd:
            return r"\ce{%s}" % species_tex(nd["sp"]), 5
        if "str" in nd:
            return r"\text{%s}" % nd["str"], 5
        if "u" in nd:
            if "n" in nd["x"]:      # a negative number
                return "-" + num(nd["x"]["n"]), 3
            t, p = self.go(nd["x"])
            return "-" + (t if p >= 3 else paren(t)), 3
        if "c" in nd:
            a, b = self(nd["l"]), self(nd["r"])
            return "%s %s %s" % (a, CMP_TEX[nd["c"]], b), 0
        if "b" in nd:
            return self.binary(nd)
        if "f" in nd:
            return self.func(nd)
        return r"\text{?}", 5

    def binary(self, nd):
        op = nd["b"]
        (a, pa), (b, pb) = self.go(nd["l"]), self.go(nd["r"])

        def par(t, p, need):
            return t if p >= need else paren(t)
        if op == "+":
            return "%s + %s" % (a, b), 1
        if op == "-":
            return "%s - %s" % (a, par(b, pb, 2)), 1
        if op in ("*", "/"):
            # a chain of products and quotients as one fraction
            nums, dens = self.factors(nd)
            top = self.product(nums) if nums else ("1", 5)
            if not dens:
                return top
            return r"\frac{%s}{%s}" % (top[0], self.product(dens)[0]), 5
        if op == "**":
            # (1/x)^p as x^{-p}
            l = self.deref(nd["l"])
            if l.get("b") == "/" and self.deref(l["l"]).get("n") == 1:
                x, px = self.go(l["r"])
                if px < 5 or "^{" in x or x.startswith((r"\frac", r"\sqrt")):
                    x = paren(x)
                return "{%s}^{-%s}" % (x, b if pb >= 5 else
                                       paren(b)), 4
            # a fraction, a power, an exponential or a root as the base
            # takes parentheses
            if pa < 5 or "^{" in a or a.startswith((r"\frac", r"\sqrt")):
                a = paren(a)
            return "{%s}^{%s}" % (a, b), 4
        return r"\text{?}", 5

    def deref(self, nd):
        """an inlined var as its definition"""
        while nd.get("kind") == "var" and nd["s"] in self.inline:
            nd = self.net.vars[nd["s"]]
        return nd

    def factors(self, nd):
        """the numerator and denominator factors of a * / chain (through
        the inlined vars)"""
        nd = self.deref(nd)
        op = nd.get("b")
        if op in ("*", "/"):
            n1, d1 = self.factors(nd["l"])
            n2, d2 = self.factors(nd["r"])
            return (n1 + n2, d1 + d2) if op == "*" else (n1 + d2, d1 + n2)
        if nd.get("n") == 1:
            return [], []
        return [nd], []

    def product(self, nodes):
        """factors as one product: numbers first, x between numbers"""
        nodes = sorted(nodes, key=lambda x: "n" not in x)
        out, prev_num = "", False
        for i, x in enumerate(nodes):
            t, p = self.go(x)
            if p < 2:
                t = paren(t)
            isnum = "n" in x
            if i:
                out += r" \times " if isnum and prev_num else r"\,"
            out += t
            prev_num = isnum
        return out, (2 if len(nodes) > 1 else self.go(nodes[0])[1])

    def func(self, nd):
        f, k = nd["f"], nd.get("kind")
        args = [self(a) for a in nd["a"]]
        if k == "density":
            return r"n_{\ce{%s}}" % species_tex(nd["a"][0]["sp"]), 5
        if f == "sqrt":
            return r"\sqrt{%s}" % args[0], 5
        if f == "abs":
            return r"\left|%s\right|" % args[0], 5
        if f == "exp":
            a, p = self.go(nd["a"][0])
            if len(a) <= 40 and r"\frac" not in a:
                return r"e^{%s}" % a, 5
            return r"\exp" + paren(a), 5
        if f == "where":
            return (r"\begin{cases} %s & %s\\ %s & \text{otherwise}"
                    r"\end{cases}" % (args[1], args[0], args[2])), 5
        head = FUNC_TEX.get(f, name_tex(f))
        return head + paren(", ".join(args)), 5


def paren(t):
    """parentheses: sized (\\left, \\right) around a fraction, plain
    otherwise, so that TeX can break a long formula inside them"""
    if r"\frac" in t or r"\sqrt" in t or r"\begin" in t:
        return r"\left(%s\right)" % t
    return "(%s)" % t


def size(nd):
    return sum(1 for _ in knet.walk(nd))


def num_text(v):
    """a number as plain text: plain when short, 1.5e+08 otherwise"""
    if float(v).is_integer() and abs(v) < 1e5:
        return "%d" % v
    return "%.6g" % v


def _trange(tr, fmt, le, ge, unit):
    """T < hi, T >= lo, or lo <= T < hi, the unit after each number;
    empty for the whole T axis"""
    if not tr:
        return ""
    lo, hi = tr
    has_lo, has_hi = lo > 0, not math.isinf(hi)
    if has_lo and has_hi:
        return "%s%s %s T < %s%s" % (fmt(lo), unit, le, fmt(hi), unit)
    if has_hi:
        return "T < %s%s" % (fmt(hi), unit)
    if has_lo:
        return "T %s %s%s" % (ge, fmt(lo), unit)
    return ""


def trange_tex(tr):
    """a T range as LaTeX math, e.g. $10\,\mathrm{K} \leq T <
    6700\,\mathrm{K}$ (all inside the math, for MathJax too)"""
    s = _trange(tr, num, r"\leq", r"\geq", r"\,\mathrm{K}")
    return "$%s$" % s if s else ""


def trange_text(tr):
    """the same as plain text, e.g. 10 K \u2264 T < 6700 K (plot labels)"""
    return _trange(tr, num_text, "\u2264", "\u2265", " K")


def rate_tex(rx, tex):
    """the rate coefficient column: the expression, or what computes it"""
    r = rx["rate"]
    k = r["kind"]
    if k == "expr":
        return "$%s$" % tex(r["ast"])
    if k == "photo":
        x = r["xsec"]
        if x["kind"] == "expr":
            s = r"photo, $\sigma(E) = %s$\,cm$^2$" % tex(x["ast"])
        elif x["kind"] == "ref":
            s = r"photo, $\sigma$ from %s" % x["ref"]
        elif x["kind"] == "h2pd":
            s = r"photo, H$_2$ line data"
        else:
            s = r"photo, $\sigma$ from the database"
        if "Emin" in r or "Emax" in r:
            s += r", $%s \leq E < %s$\,eV" % (num(r.get("Emin", 0.0)),
                                             num(r.get("Emax", math.inf))
                                             if "Emax" in r else r"\infty")
        return s
    if k == "surface":
        if "Ea" in r:
            return r"surface, $E_\mathrm{a} = %s$\,K" % num(r["Ea"])
        if "prob" in r:
            return r"surface, probability $%s$" % num(r["prob"])
        if "stick" in r:
            return r"adsorption, $S = %s$" % tex(r["stick"])
        return "surface model"
    if k == "auto":
        return "from the database"
    return k


def escape(s):
    return (s.replace("\\", r"\textbackslash{}").replace("_", r"\_")
            .replace("%", r"\%").replace("&", r"\&").replace("#", r"\#"))


def table(net, path_body, path_doc, title):
    """the reaction table: path_body (the longtable, for \\input) and
    path_doc (a standalone document around it, with the definitions of the
    vars written by name)"""
    tex = Tex(net)
    rows = []
    prev = None
    for rx in net.reactions:
        same = prev is not None and rx["eq"] == prev["eq"] and \
            rx.get("group") is not None and \
            rx.get("group") == prev.get("group")
        ref = rx.get("ref", "")
        rows.append("%s & %s & %s & %s & %s \\\\" % (
            "" if same else str(rx["id"]),
            "" if same else reaction_math(rx),
            rate_tex(rx, tex), trange_tex(rx["trange"]),
            ref if ref.lstrip().startswith("\\") else escape(ref)))
        prev = rx
    body = [r"% generated by kromex_doc from " + escape(net.name),
            r"\begin{longtable}{rlp{0.45\textwidth}ll}",
            r"\caption{The reactions of %s. Rate coefficients in "
            r"cm$^{3(n-1)}$\,s$^{-1}$ for $n$ reactants; $T$ in K.}\\"
            % escape(net.name),
            r"\hline No. & Reaction & Rate coefficient & $T$ range & Ref. "
            r"\\ \hline", r"\endfirsthead",
            r"\hline No. & Reaction & Rate coefficient & $T$ range & Ref. "
            r"\\ \hline", r"\endhead", r"\hline", r"\endfoot"]
    body += rows + [r"\end{longtable}"]
    path_body.write_text("\n".join(body) + "\n")
    used = set()
    for rx in net.reactions:
        if rx["rate"]["kind"] == "expr":
            for nd in knet.walk(rx["rate"]["ast"]):
                if nd.get("kind") == "var":
                    used.add(nd["s"])
    # the vars named in the table, and those they name, in file order
    todo = list(used)
    while todo:
        v = todo.pop()
        for nd in knet.walk(net.vars[v]):
            if nd.get("kind") == "var" and nd["s"] not in used:
                used.add(nd["s"])
                todo.append(nd["s"])
    defs = [v for v in net.var_order if v in used and v not in tex.inline]
    doc = [r"\documentclass{article}",
           r"\usepackage[a4paper,margin=2cm,landscape]{geometry}",
           r"\usepackage{amsmath,longtable}",
           r"\begin{document}", r"\section*{%s}" % escape(title),
           r"\input{%s}" % path_body.name]
    if defs:
        doc += [r"\subsection*{Definitions}", r"\begin{longtable}{ll}"]
        doc += [r"$%s$ & $= %s$ \\" % (name_tex(v), tex(net.vars[v]))
                for v in defs]
        doc += [r"\end{longtable}"]
    doc += [r"$T$: the gas temperature [K]; $n_X$: the number density of "
            r"$X$ [cm$^{-3}$].", r"\end{document}"]
    path_doc.write_text("\n".join(doc) + "\n")
