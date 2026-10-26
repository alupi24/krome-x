"""The rate coefficients of a parsed network, evaluated in numpy from the
expressions of the .knet (no build): over a grid of Tgas or of another
quantity, the rest held at reference values.

What only the generated code computes is not evaluated here: photo and
surface rates, rates from the database (auto), tables and extern
functions. Their reactions are documented without plots."""

import math

import numpy as np

import knet

KB_EV = knet.CONSTANTS["boltzmann_eV"]
TMAX_RATE = 1e9   # the generated code evaluates rates at min(max(T, Tcmb), 1e9)


class NotEvaluable(Exception):
    pass


def sigma_v96(E, E0, s0, ya, P, yw, y0, y1):
    """Verner et al. 1996 photoionisation fit [cm^2], as the generated code
    (codegen.py emit_photo); networks call it in photo(xsec=...)"""
    x = E / E0 - y0
    y = np.sqrt(x * x + y1 * y1)
    return 1e-18 * s0 * ((x - 1.0) ** 2 + yw * yw) * y ** (0.5 * P - 5.5) \
        * (1.0 + np.sqrt(y / ya)) ** (-P)


# extern functions that are pure formulas, evaluated here too
PURE_EXTERNS = {"sigma_v96": sigma_v96}


def ast_deps(node, net, seen=None):
    """the names an expression depends on, through its vars: states and
    params by name, densities as "n(X)", tables and externs as "f()" """
    out = set()
    seen = set() if seen is None else seen
    for nd in knet.walk(node):
        if "s" in nd:
            k = nd.get("kind")
            if k == "var":
                if nd["s"] not in seen:
                    seen.add(nd["s"])
                    out |= ast_deps(net.vars[nd["s"]], net, seen)
            elif k in ("state", "param"):
                out.add(nd["s"])
        elif "f" in nd:
            k = nd.get("kind")
            if k == "density":
                out.add("n(%s)" % nd["a"][0]["sp"])
            elif k in ("table", "extern"):
                out.add("%s()" % nd["f"])
    return out


def rate_deps(rx, net):
    """the dependencies of a reaction's rate (None when it is not an
    expression this module evaluates)"""
    r = rx["rate"]
    if r["kind"] != "expr":
        return None
    d = ast_deps(r["ast"], net)
    if net.rate_scale is not None:
        d |= ast_deps(net.rate_scale, net)
    return d


def evaluable(rx, net):
    d = rate_deps(rx, net)
    return d is not None and not any(x.endswith("()") for x in d)


class Evaluator:
    """values: {name: scalar or array} for states, params and densities
    ("n(X)"); the missing params take their defaults"""

    def __init__(self, net, values):
        self.net = net
        self.v = dict(knet.IMPLICIT_PARAMS)
        self.v.update(net.params)
        self.v.update(values)
        T = np.asarray(self.v.get("Tgas", 100.0), dtype=float)
        Tcmb = float(self.v.get("Tcmb", 2.725))
        self.v["Tgas"] = np.minimum(np.maximum(T, Tcmb), TMAX_RATE)
        self.v.setdefault("Tdust", 10.0)
        self.v["Tgas_eV"] = self.v["Tgas"] * KB_EV
        self.v.setdefault("Hnuclei", 1.0)
        self.v.setdefault("ntot", self.v["Hnuclei"])
        self.v.setdefault("J_LW", 0.0)
        self.cache = {}

    def var(self, name):
        if name not in self.cache:
            self.cache[name] = self.ev(self.net.vars[name])
        return self.cache[name]

    def ev(self, nd):
        if "n" in nd:
            return np.float64(nd["n"])
        if "s" in nd:
            k, s = nd.get("kind"), nd["s"]
            if k == "var":
                return self.var(s)
            if k == "const":
                return np.float64(knet.CONSTANTS[s])
            if s in self.v:
                return np.asarray(self.v[s], dtype=float)
            raise NotEvaluable("no value for %s" % s)
        if "u" in nd:
            return -self.ev(nd["x"])
        if "b" in nd:
            a, b, op = self.ev(nd["l"]), self.ev(nd["r"]), nd["b"]
            if op == "+":
                return a + b
            if op == "-":
                return a - b
            if op == "*":
                return a * b
            if op == "/":
                return a / b
            if op == "**":
                return np.power(a, b)
        if "c" in nd:
            a, b, op = self.ev(nd["l"]), self.ev(nd["r"]), nd["c"]
            return {">": np.greater, "<": np.less, ">=": np.greater_equal,
                    "<=": np.less_equal, "==": np.equal,
                    "!=": np.not_equal}[op](a, b)
        if "f" in nd:
            k, f = nd.get("kind"), nd["f"]
            if k == "density":
                key = "n(%s)" % nd["a"][0]["sp"]
                if key not in self.v:
                    raise NotEvaluable("no value for %s" % key)
                return np.asarray(self.v[key], dtype=float)
            if k == "extern" and f in PURE_EXTERNS:
                return PURE_EXTERNS[f](*[self.ev(x) for x in nd["a"]])
            if k != "builtin":
                raise NotEvaluable("%s() is computed by the generated code"
                                   % f)
            a = [self.ev(x) for x in nd["a"]]
            if f == "where":
                return np.where(a[0], a[1], a[2])
            fn = {"exp": np.exp, "log": np.log, "log10": np.log10,
                  "sqrt": np.sqrt, "abs": np.abs, "tanh": np.tanh,
                  "min": np.minimum, "max": np.maximum}[f]
            return fn(*a)
        raise NotEvaluable("unsupported expression")

    def rate(self, rx):
        """the rate coefficient of a reaction as written (no T range)"""
        r = rx["rate"]
        if r["kind"] != "expr":
            raise NotEvaluable("a %s rate is computed by the generated code"
                               % r["kind"])
        with np.errstate(all="ignore"):
            k = self.ev(r["ast"])
            if self.net.rate_scale is not None:
                k = k * self.ev(self.net.rate_scale)
        return k


def in_range(T, trange):
    """the half-open [Tmin, Tmax) of a piece (all T without a range)"""
    if not trange:
        return np.ones_like(T, dtype=bool)
    return (T >= trange[0]) & (T < trange[1])


def joints(net, group, values):
    """the relative jumps of a piecewise rate at its interior bounds:
    [(T, k_left, k_right, |dk| / k_left)]"""
    out = []
    for a, b in zip(group, group[1:]):
        Tb = a["trange"][1]
        try:
            kl = float(Evaluator(net, dict(values, Tgas=Tb)).rate(a))
            kr = float(Evaluator(net, dict(values, Tgas=Tb)).rate(b))
        except NotEvaluable:
            continue
        err = abs(kl - kr) / abs(kl) if kl != 0.0 else (
            0.0 if kr == 0.0 else math.inf)
        out.append((Tb, kl, kr, err))
    return out
