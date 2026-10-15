"""The HTML documentation of a network: index, species, reactions (with
their rate plots), piecewise rates, missing chemistry, graphs. Static
pages; formulas by MathJax and graphs without Graphviz by viz.js, both
from a CDN."""

import datetime
import math
from html import escape

import numpy as np

from . import evaluate, graph, latex, missing, photo, svgplot

CSS = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 0;
       color: #222; background: #fff; }
header { background: #23395b; color: #fff; padding: 10px 20px; }
header a { color: #fff; margin-right: 16px; text-decoration: none; }
header b { margin-right: 24px; }
main { padding: 16px 20px; max-width: 1200px; }
table { border-collapse: collapse; margin: 8px 0 16px; }
th, td { border-bottom: 1px solid #ddd; padding: 3px 8px; text-align: left;
         vertical-align: top; font-size: 14px; }
th { background: #f0f0f0; }
tr:nth-child(even) td { background: #fafafa; }
code, pre { background: #f4f4f4; padding: 1px 4px; font-size: 13px; }
pre { padding: 8px; overflow-x: auto; }
.warn { color: #b00; }
.muted { color: #777; }
.small { font-size: 12px; }
"""

MATHJAX = """<script>window.MathJax = {tex: {packages: {'[+]': ['mhchem']},
inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}, loader: {load: ['[tex]/mhchem']},
svg: {fontCache: 'global'}};</script>
<script async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>"""


# interactive plots: plotly.js replaces each SVG (div.kplot) with the
# same data, zoomable; without it (offline, no JavaScript) the SVG stays
PLOTLY = """<script>function kplots() {
  document.querySelectorAll('.kplot').forEach(function (el) {
    var f = JSON.parse(el.dataset.fig);
    el.innerHTML = '';
    Plotly.newPlot(el, f.data, f.layout, {displaylogo: false,
                                         responsive: true});
  });
}</script>
<script defer onload="kplots()"
 src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>"""


def page(title, body, depth=0, net_name=""):
    up = "../" * depth
    nav = "".join('<a href="%s%s">%s</a>' % (up, href, lab) for href, lab in (
        ("index.html", "Home"), ("species.html", "Species"),
        ("reactions.html", "Reactions"), ("photo.html", "Photo"),
        ("piecewise.html", "Piecewise rates"),
        ("missing.html", "Missing reactions"),
        ("branches.html", "Missing branches"), ("graphs.html", "Graphs"),
        ("latex.html", "LaTeX")))
    return ("<!doctype html><html><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width, "
            "initial-scale=1'><title>%s</title><link rel='stylesheet' "
            "href='%sstyle.css'>%s%s</head><body><header><b>%s</b>%s"
            "</header><main>%s</main></body></html>"
            % (escape(title), up, MATHJAX,
               PLOTLY if 'class="kplot"' in body else "",
               escape(net_name), nav, body))


def ce(sp):
    return r"\(\ce{%s}\)" % escape(latex.species_tex(sp))


def fmt(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    return "%.3g" % v


def k_label(rx):
    """k with its units: cm^{3(n-1)} s^-1 for n reactants"""
    n = sum(1 for r in rx["reactants"] if r.get("kinetic", True))
    sup = {0: "", 1: "³", 2: "⁶", 3: "⁹"}
    return "k [%ss⁻¹]" % ("cm%s " % sup.get(n - 1, "^%d" % (3 * (n - 1)))
                          if n > 1 else "")


class Site:
    def __init__(self, net, source, outdir, values, ranges, tref,
                 thermo, max_graph=graph.MAX_NODES, plots=True,
                 outside="asym", interactive=True):
        self.net = net
        self.source = source
        self.lines = source.read_text().splitlines()
        self.out = outdir
        self.values = values           # reference values {name: value}
        self.ranges = ranges           # {name: (lo, hi)}, Tgas first
        self.tref = tref
        self.thermo = thermo
        self.max_graph = max_graph
        self.plots = plots
        self.interactive = interactive
        # outside a rate's T range (rates.outside): asym, zero below and
        # the bound's value above (the default); clamp; zero
        self.outside = outside
        self.tex = latex.Tex(net)
        self.sp = sorted(net.species)
        self.sidx = {s: i for i, s in enumerate(self.sp)}
        self.groups = {}
        for rx in net.reactions:
            if rx.get("group") is not None:
                self.groups.setdefault(rx["group"], []).append(rx)
        for g in self.groups.values():
            g.sort(key=lambda r: r["trange"][0])
        self.by_eq = {}
        for rx in net.reactions:
            self.by_eq.setdefault(rx["eq"], []).append(rx)
        self.xs, self.xs_warn = photo.xsecs(net)
        self.kref = {}
        self.reason = {}
        for rx in net.reactions:
            self.kref[rx["id"]] = self.k_at(rx, tref)

    # -- helpers --------------------------------------------------------------

    def w(self, rel, html):
        p = self.out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(html)

    def page(self, rel, title, body):
        self.w(rel, page(title, body, rel.count("/"), self.net.name))

    def slink(self, sp, depth):
        if sp not in self.sidx:
            return escape(sp)
        return '<a href="%sspecies/s%d.html">%s</a>' % (
            "../" * depth, self.sidx[sp], ce(sp))

    def rlink(self, rx, depth, text=None):
        return '<a href="%sreactions/r%d.html">%s</a>' % (
            "../" * depth, rx["id"],
            text if text is not None else r"\(%s\)" % escape(
                latex.reaction_tex(rx)))

    def k_at(self, rx, T):
        """k at the reference values and T as KROME-X evaluates it outside
        the group's range (rates.outside), or None"""
        if not evaluate.evaluable(rx, self.net):
            self.reason[rx["id"]] = self.why_not(rx)
            return None
        tr = rx["trange"]
        if tr and not (tr[0] <= T < tr[1]):
            # another piece of the group covers T, or outside the union
            g = self.groups.get(rx.get("group"), [rx])
            lo = min(r["trange"][0] for r in g)
            hi = max(r["trange"][1] for r in g)
            if T < lo or T >= hi:
                below = T < lo
                if self.outside == "zero" or (below and
                                              self.outside == "asym"):
                    return 0.0
                end = (min(g, key=lambda r: r["trange"][0]) if below else
                       max(g, key=lambda r: r["trange"][1]))
                if end is not rx:
                    return None
                T = lo if below else hi
            else:
                return None
        try:
            v = evaluate.Evaluator(self.net, dict(self.values, Tgas=T))
            return float(v.rate(rx))
        except evaluate.NotEvaluable as e:
            self.reason[rx["id"]] = str(e)
            return None

    def why_not(self, rx):
        k = rx["rate"]["kind"]
        if k != "expr":
            return "a %s rate: computed by the generated code" % k
        d = evaluate.rate_deps(rx, self.net)
        f = sorted(x for x in d if x.endswith("()"))
        return "uses %s: computed by the generated code" % ", ".join(f)

    def zero_params(self, rx):
        d = evaluate.rate_deps(rx, self.net) or set()
        v = dict(evaluate.knet.IMPLICIT_PARAMS)
        v.update(self.net.params)
        v.update(self.values)
        return sorted(x for x in d if x in v and v[x] == 0
                      and x not in self.ranges)

    # -- the pages ------------------------------------------------------------

    def build(self, miss, graphs, latex_files):
        self.w("style.css", CSS)
        self.index(miss, graphs)
        self.species_pages()
        self.reaction_pages()
        self.photo_page()
        self.piecewise_page()
        self.missing_pages(miss)
        self.graphs_page(graphs)
        self.latex_page(latex_files)

    def index(self, miss, graphs):
        net = self.net
        n_eval = sum(1 for rx in net.reactions
                     if evaluate.evaluable(rx, net))
        kinds = {}
        for rx in net.reactions:
            kinds[rx["rate"]["kind"]] = kinds.get(rx["rate"]["kind"], 0) + 1
        rows = [("file", "<code>%s</code>" % escape(str(self.source))),
                ("species", len(net.species)),
                ("reactions", "%d (%s)" % (len(net.reactions), ", ".join(
                    "%d %s" % (v, k) for k, v in sorted(kinds.items())))),
                ("rates evaluated here", "%d (the others: photo, surface, "
                 "database, tables, extern functions; computed by the "
                 "generated code)" % n_eval),
                ("parameters", ", ".join(sorted(net.params)) or "–"),
                ("vars", len(net.vars)), ("tables", ", ".join(
                    sorted(net.tables)) or "–"),
                ("extern functions", ", ".join(sorted(net.externs)) or "–")]
        body = ["<h1>%s</h1>" % escape(net.name),
                "<table>%s</table>" % "".join(
                    "<tr><th>%s</th><td>%s</td></tr>" % (k, v)
                    for k, v in rows)]
        body.append("<h2>Evaluation</h2><p>Rate coefficients evaluated from "
                    "the expressions of the network file at T = %g K (the "
                    "tables) and over %s (the plots), with:</p>"
                    % (self.tref, ", ".join("%s %g–%g" % (k, *v)
                                           for k, v in self.ranges.items())))
        vals = dict(evaluate.knet.IMPLICIT_PARAMS)
        vals.update(net.params)
        vals.update(self.values)
        used = set()
        for rx in net.reactions:
            used |= evaluate.rate_deps(rx, net) or set()
        body.append("<table><tr><th>quantity</th><th>value</th></tr>%s"
                    "</table>" % "".join(
                        "<tr><td>%s</td><td>%s</td></tr>" % (
                            escape(k), fmt(vals[k]) if k in vals else
                            '<span class="muted">not set</span>')
                        for k in sorted(used) if not k.endswith("()")
                        and k != "Tgas"))
        if net.notices:
            body.append("<h2>Parser notices</h2><ul>%s</ul>" % "".join(
                "<li>%s</li>" % escape(x) for x in net.notices))
        body.append('<p class="small muted">Generated by kromex_doc '
                    '(KROME-X) on %s.</p>'
                    % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.page("index.html", net.name, "\n".join(body))

    def species_pages(self):
        net, th = self.net, self.thermo
        form = {s: [] for s in self.sp}
        dest = {s: [] for s in self.sp}
        for rx in net.reactions:
            for r in {r["sp"] for r in rx["reactants"]}:
                dest.setdefault(r, []).append(rx)
            for p in {p["sp"] for p in rx["products"]}:
                form.setdefault(p, []).append(rx)
        rows = []
        for i, s in enumerate(self.sp):
            d = net.species[s]
            comp = ("opaque" if d["composition"] is None else " ".join(
                "%s%s" % (e, n if n > 1 else "")
                for e, n in sorted(d["composition"].items())) or "–")
            h = th.h298(s)
            rows.append("<tr><td>%s</td><td>%s</td><td>%.4g</td><td>%+d</td>"
                        "<td>%s</td><td>%s</td><td>%d</td><td>%d</td></tr>"
                        % (self.slink(s, 0), d["phase"], d["mass_amu"],
                           d["charge"], comp, fmt(h), len(form[s]),
                           len(dest[s])))
            body = ["<h1>%s</h1>" % ce(s), "<table>",
                    "<tr><th>name</th><td><code>%s</code></td></tr>" % escape(s),
                    "<tr><th>phase</th><td>%s</td></tr>" % d["phase"],
                    "<tr><th>mass</th><td>%.6g amu</td></tr>" % d["mass_amu"],
                    "<tr><th>charge</th><td>%+d</td></tr>" % d["charge"],
                    "<tr><th>composition</th><td>%s</td></tr>" % comp,
                    "<tr><th>&Delta;H<sub>f</sub>(298.15 K)</th><td>%s</td>"
                    "</tr>" % ("%.4g kJ/mol" % h if h is not None else
                               '<span class="muted">not in the tables</span>'),
                    "</table>"]
            for lab, lst in (("Formation", form[s]), ("Destruction", dest[s])):
                body.append("<h2>%s (%d)</h2>" % (lab, len(lst)))
                body.append(self.reaction_table(lst, 1, ref_species=s))
            self.page("species/s%d.html" % i, s, "\n".join(body))
        self.page("species.html", "Species", (
            "<h1>Species (%d)</h1><p class='muted'>&Delta;H<sub>f</sub>: the "
            "formation enthalpy at 298.15 K [kJ/mol] (data/thermo_burcat"
            ".toml: Burcat's database, Goos, Burcat &amp; Ruscic; "
            "data/thermoNIST.dat).</p><table><tr><th>species</th><th>phase"
            "</th><th>mass [amu]</th><th>charge</th><th>composition</th>"
            "<th>&Delta;H<sub>f</sub></th><th>formed by</th>"
            "<th>destroyed by</th></tr>%s</table>"
            % (len(self.sp), "".join(rows))))

    def reaction_table(self, lst, depth, ref_species=None):
        rows = []
        for rx in lst:
            k = self.kref[rx["id"]]
            rows.append("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td>"
                        "<td>%s</td><td class='small'>%s</td></tr>" % (
                            rx["id"], self.rlink(rx, depth),
                            escape(latex.trange_tex(rx["trange"])),
                            rx["rate"]["kind"], fmt(k) if k is not None
                            else '<span class="muted">–</span>',
                            escape(rx.get("ref", ""))))
        return ("<table><tr><th>#</th><th>reaction</th><th>T range</th>"
                "<th>rate</th><th>k(%g K)</th><th>ref</th></tr>%s</table>"
                % (self.tref, "".join(rows)))

    def reaction_pages(self):
        net = self.net
        for rx in net.reactions:
            self.page("reactions/r%d.html" % rx["id"], rx["eq"],
                      self.reaction_body(rx))
        self.page("reactions.html", "Reactions",
                  "<h1>Reactions (%d)</h1>%s" % (
                      len(net.reactions),
                      self.reaction_table(net.reactions, 0)))

    def reaction_body(self, rx):
        net = self.net
        b = ["<h1>%d. %s</h1>" % (rx["id"], r"\(%s\)" % escape(
            latex.reaction_tex(rx)))]
        src = self.lines[rx["line"] - 1] if rx["line"] <= len(self.lines) \
            else ""
        b.append("<p>Species: %s</p>" % ", ".join(
            self.slink(s, 1) for s in dict.fromkeys(
                [r["sp"] for r in rx["reactants"]] +
                [p["sp"] for p in rx["products"]])))
        b.append("<table>")
        b.append("<tr><th>rate</th><td>%s</td></tr>" % escape(
            latex.rate_tex(rx, self.tex)))
        b.append("<tr><th>T range</th><td>%s</td></tr>" % (escape(
            latex.trange_tex(rx["trange"])) or "all T"))
        if rx.get("ref"):
            b.append("<tr><th>reference</th><td>%s</td></tr>"
                     % escape(rx["ref"]))
        b.append("<tr><th>source</th><td><code>%s:%d</code><pre>%s</pre>"
                 "</td></tr>" % (escape(self.source.name), rx["line"],
                                 escape(src)))
        deps = evaluate.rate_deps(rx, net)
        if deps:
            b.append("<tr><th>depends on</th><td>%s</td></tr>"
                     % escape(", ".join(sorted(deps))))
        b.append("</table>")
        # the vars written by name in the formula
        names = []
        if rx["rate"]["kind"] == "expr":
            for nd in evaluate.knet.walk(rx["rate"]["ast"]):
                if nd.get("kind") == "var" and nd["s"] not in self.tex.inline \
                        and nd["s"] not in names:
                    names.append(nd["s"])
        if names:
            b.append("<p>with</p><ul>%s</ul>" % "".join(
                "<li>%s</li>" % escape("$%s = %s$" % (
                    latex.name_tex(v), self.tex(net.vars[v])))
                for v in names))
        others = [r for r in self.by_eq[rx["eq"]] if r is not rx]
        if others:
            b.append("<p>Same reaction, other lines: %s</p>" % ", ".join(
                self.rlink(r, 1, "#%d (%s, %s)" % (
                    r["id"], escape(latex.trange_tex(r["trange"])) or "all T",
                    "same piecewise rate" if r.get("group") is not None
                    and r.get("group") == rx.get("group") else
                    "additive channel")) for r in others))
        b.append(self.plots_html(rx))
        return "\n".join(b)

    def plots_html(self, rx):
        net = self.net
        if rx["id"] in self.xs:
            return self.photo_html(self.xs[rx["id"]])
        if not evaluate.evaluable(rx, net):
            return '<p class="muted">No plot: %s.</p>' % escape(
                self.why_not(rx))
        if not self.plots:
            return ""
        deps = evaluate.rate_deps(rx, net)
        out = []
        g = self.groups.get(rx.get("group"), [rx]) if rx["trange"] else [rx]
        try:
            svg = self.t_plot(g, rx)
        except evaluate.NotEvaluable as e:
            return '<p class="muted">No plot: %s (set it with --set).</p>' \
                % escape(str(e))
        if svg:
            out.append(svg)
        elif "Tgas" in deps or "Tgas_eV" in deps:
            z = self.zero_params(rx)
            out.append('<p class="muted">Zero over the range%s.</p>' % (
                " (%s = 0 at the reference values; set with --set)"
                % ", ".join(z) if z else ""))
        for var, (lo, hi) in self.ranges.items():
            if var == "Tgas" or var not in deps:
                continue
            x = np.logspace(math.log10(lo), math.log10(hi), 120)
            try:
                k = evaluate.Evaluator(net, dict(self.values, Tgas=self.tref,
                                                 **{var: x})).rate(rx)
            except evaluate.NotEvaluable:
                continue
            k = np.broadcast_to(k, x.shape)
            s = svgplot.figure([{"x": x, "y": k}], var, k_label(rx),
                               title="k at T = %g K" % self.tref,
                               interactive=self.interactive)
            if s:
                out.append(s)
        if not any("Tgas" in d for d in [deps]) and not ("Tgas_eV" in deps) \
                and not any(v in deps for v in self.ranges):
            k = self.kref[rx["id"]]
            z = self.zero_params(rx)
            out.append('<p>k = %s (independent of T%s)%s.</p>' % (
                fmt(k), " and of the plotted ranges" if len(self.ranges) > 1
                else "", ": %s = 0 at the reference values (set with "
                "--set)" % ", ".join(z) if k == 0 and z else ""))
        return "\n".join(out)

    def t_plot(self, group, this):
        """k(T) of a piecewise group: each piece in its range (solid), its
        formula outside (faint dashed), KROME-X's rate outside the union
        (zero below, the bound's value above: dotted)"""
        lo, hi = self.ranges["Tgas"]
        T = np.logspace(math.log10(lo), math.log10(hi), 240)
        series, pieces = [], []
        for i, rx in enumerate(group):
            tr = rx["trange"]
            Tp = T
            if tr:
                inner = [t for t in (tr[0], tr[1]) if lo < t < hi]
                Tp = np.unique(np.concatenate([T, inner]))
            k = np.broadcast_to(evaluate.Evaluator(
                self.net, dict(self.values, Tgas=Tp)).rate(rx), Tp.shape)
            m = evaluate.in_range(Tp, tr) if tr else np.ones_like(Tp, bool)
            if tr and not math.isinf(tr[1]):
                m = m | (Tp == tr[1])
            color = svgplot.COLORS[i % len(svgplot.COLORS)]
            lab = ("#%d %s" % (rx["id"], latex.trange_tex(tr).replace(
                "$", "").replace(r"\leq", "≤").replace(r"\,K", " K")
                .replace(r"\times10^", "e")) if len(group) > 1 else None)
            if tr:
                series.append({"x": Tp, "y": np.where(m, np.nan, k),
                               "style": "dashed", "color": color,
                               "opacity": 0.35, "width": 1.0})
            marks = []
            if tr:
                for t in tr:
                    if lo <= t <= hi and t > 0:
                        j = np.argmin(np.abs(Tp - t))
                        marks.append((Tp[j], k[j]))
            series.append({"x": Tp[m], "y": k[m], "color": color,
                           "width": 2.4 if rx is this else 1.6,
                           "label": lab, "markers": marks})
            pieces.append((tr, Tp, k))
        # KROME-X outside the union of the ranges (rates.outside)
        notes = []
        if all(p[0] for p in pieces):
            ends = (("below", min(pieces, key=lambda p: p[0][0]), 0),
                    ("above", max(pieces, key=lambda p: p[0][1]), 1))
            for side, piece, j in ends:
                Tb = piece[0][j]
                if not (lo < Tb < hi):
                    continue
                if self.outside == "zero" or (side == "below" and
                                              self.outside == "asym"):
                    notes.append("%s the range KROME-X uses k = 0" % side)
                    continue
                kb = float(np.broadcast_to(evaluate.Evaluator(
                    self.net, dict(self.values, Tgas=Tb)).rate(
                        group[pieces.index(piece)]), ()))
                Tc = T[T >= Tb] if side == "above" else T[T < Tb]
                series.append({"x": Tc, "y": np.full(Tc.shape, kb),
                               "style": "dotted", "color": "#000",
                               "label": "KROME-X %s the range" % side})
        note = "; ".join(notes) or None
        return svgplot.figure(series, "T [K]", k_label(this),
                              xrange=(lo, hi), note=note,
                              interactive=self.interactive)

    def photo_html(self, x):
        """a photo reaction: the cross section's source and window, the
        shielding, the rate in the Draine field, sigma(E)"""
        rows = [("cross section", escape(x.source)),
                ("energy window", "%s – %s eV" % (
                    fmt(x.Emin), "∞" if math.isinf(x.Emax) else fmt(x.Emax)))]
        if x.shield is not None:
            rows.append(("shielding factor", escape("$%s$" % self.tex(
                x.shield))))
        if x.ok:
            rows.append(("k, unshielded Draine field", "%s s<sup>-1</sup> "
                         "<span class='muted small'>(4&pi;&int;F&sigma;dE, "
                         "Draine 1978, 5–13.6 eV)</span>"
                         % fmt(x.draine_rate())))
        else:
            rows.append(("data", '<span class="warn">%s</span>'
                         % escape(x.reason)))
        out = ["<h2>Photo reaction</h2><table>%s</table>" % "".join(
            "<tr><th>%s</th><td>%s</td></tr>" % r for r in rows)]
        if x.ok and self.plots:
            lo = max(min(x.Emin if x.Emin > 0 else 5.0, 5.0) * 0.8, 0.5)
            hi = x.Emax * 1.2 if math.isfinite(x.Emax) else \
                max(1e3, 10 * lo)
            E = x.grid(lo, hi)
            svg = svgplot.figure([{"x": E, "y": x(E)}], "E [eV]",
                                 "σ [cm²]", title="cross section",
                                 decades=8.0, interactive=self.interactive)
            out.append(svg or '<p class="muted">The cross section is zero '
                       'over %g–%g eV.</p>' % (lo, hi))
        return "\n".join(out)

    def photo_page(self):
        rows = []
        for rid, x in sorted(self.xs.items()):
            rx = x.rx
            k = x.draine_rate() if x.ok else None
            rows.append("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td>"
                        "<td>%s</td><td>%s</td></tr>" % (
                            rid, self.rlink(rx, 0), escape(x.source),
                            fmt(x.Emin), fmt(k) if k is not None else
                            '<span class="warn">%s</span>'
                            % escape(x.reason),
                            "yes" if x.shield is not None else ""))
        nmiss = sum(1 for x in self.xs.values() if not x.ok)
        body = ["<h1>Photo reactions (%d)</h1>" % len(self.xs),
                "<p>Cross sections from the sources the code generator "
                "reads (data/database), and the unshielded rate in the "
                "Draine (1978) field, 5–13.6 eV. %s</p>" % (
                    '<span class="warn">%d without data: the generated code '
                    'gives them a zero rate.</span>' % nmiss if nmiss
                    else "Every cross section has its data."),
                "<table><tr><th>#</th><th>reaction</th><th>cross section"
                "</th><th>E<sub>min</sub> [eV]</th><th>k Draine [s<sup>-1"
                "</sup>]</th><th>shielded</th></tr>%s</table>"
                % "".join(rows)]
        if self.xs_warn:
            body.append("<h2>Notes on the data files</h2><ul>%s</ul>"
                        % "".join("<li>%s</li>" % escape(w)
                                  for w in self.xs_warn))
        self.page("photo.html", "Photo reactions", "\n".join(body))

    def piecewise_page(self):
        rows = []
        for gid, g in sorted(self.groups.items()):
            if len(g) < 2:
                continue
            js = evaluate.joints(self.net, g, self.values)
            worst = max((j[3] for j in js), default=None)
            flag = worst is not None and worst > 1e-3
            rows.append("<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td>"
                        "</tr>" % (
                            self.rlink(g[0], 0), len(g), " ".join(
                                "T=%g: %s/%s" % (j[0], fmt(j[1]), fmt(j[2]))
                                for j in js) or "–",
                            ('<span class="warn">%.2g%%</span>' if flag
                             else "%.2g%%") % (100 * worst)
                            if worst is not None else "–"))
        self.page("piecewise.html", "Piecewise rates", (
            "<h1>Piecewise rates</h1><p>Reactions written in pieces over T "
            "ranges (one rate). At each interior bound the two pieces should "
            "agree: the jump |k<sub>left</sub> &minus; k<sub>right</sub>| / "
            "k<sub>left</sub> is shown, in red above 0.1%%.</p><table><tr>"
            "<th>reaction</th><th>pieces</th><th>k left / right at the "
            "bounds</th><th>largest jump</th></tr>%s</table>"
            % "".join(rows)))

    def missing_pages(self, miss):
        def cand_rows(lst):
            return "".join(
                "<tr><td></td><td>&rarr; %s</td><td>%s</td><td>%s</td></tr>"
                % (" + ".join(self.slink(s, 0) for s in c["species"]),
                   "%.0f" % c["dH"] if c["dH"] is not None else "–",
                   missing.mark(c["dH"])) for c in lst)
        head = ("<tr><th>reactants</th><th>products</th><th>&Delta;H [K]"
                "</th><th></th></tr>")
        intro = ("<p>Candidates have the reactants' elements and charge (one "
                 "or two species of the same phase, no cation-anion pair); "
                 "&Delta;H at 298.15 K from the thermochemical tables, "
                 "negative when exothermic. A lead for the network's author, "
                 "not a verdict: a candidate may be forbidden (spin, "
                 "barriers) or negligible.</p>")
        skipped = "".join("<p class='muted'>%s</p>" % escape(x)
                          for x in miss["skipped"])
        body = ["<h1>Missing reactions</h1>", intro,
                "<p>Reactant sets (one species, or two not of the same "
                "charge sign) that no reaction of the network uses.</p>",
                skipped]
        for ph, d in miss["phases"].items():
            body.append("<h2>%s phase (%d reactant sets)</h2><table>%s"
                        % (ph, len(d["missing"]), head))
            for m in d["missing"]:
                body.append("<tr><td>%s</td><td colspan=3></td></tr>%s" % (
                    " + ".join(self.slink(s, 0) for s in m["reactants"]),
                    cand_rows(m["products"])))
            body.append("</table>")
        self.page("missing.html", "Missing reactions", "\n".join(body))
        body = ["<h1>Missing branches</h1>", intro,
                "<p>Reactant sets of the network with products it does not "
                "have.</p>", skipped]
        for ph, d in miss["phases"].items():
            body.append("<h2>%s phase (%d reactant sets)</h2><table>%s"
                        % (ph, len(d["branches"]), head))
            for m in d["branches"]:
                body.append(
                    "<tr><td>%s</td><td colspan=3 class='small muted'>in the "
                    "network: %s</td></tr>%s" % (
                        " + ".join(self.slink(s, 0) for s in m["reactants"]),
                        "; ".join(" + ".join(p["species"])
                                  for p in m["present"]),
                        cand_rows(m["missing"])))
            body.append("</table>")
        self.page("branches.html", "Missing branches", "\n".join(body))

    def graphs_page(self, graphs):
        done, skipped, rendered = graphs
        body = ["<h1>Graphs</h1><p>One graph per element: the species with "
                "it, an arrow from each reactant to each product with it, "
                "labelled with the other reactants. Graphviz files in "
                "<code>graphs/</code>.</p>"]
        if not rendered:
            body.append("<p class='muted'>Graphviz (<code>dot</code>) was "
                        "not found: the graphs are drawn in the browser by "
                        "viz.js.</p>")
        for el, n in skipped:
            body.append("<p class='muted'>%s: %d species, above %d; .dot "
                        "file not written.</p>" % (el, n, self.max_graph))
        for el, svg, n in done:
            body.append("<h2>%s (%d species)</h2>" % (escape(el), n))
            if svg:
                body.append('<img src="graphs/%s" style="max-width:100%%">'
                            % svg)
            else:
                dot = (self.out / "graphs" / ("%s.dot" % el)).read_text()
                body.append('<div class="viz" data-dot="%s"></div>'
                            % escape(dot, quote=True))
        if not rendered:
            body.append(
                '<script src="https://cdn.jsdelivr.net/npm/@viz-js/viz@3.11.0/'
                'lib/viz-standalone.js"></script><script>'
                'Viz.instance().then(function(viz){document.querySelectorAll('
                '".viz").forEach(function(el){try{el.appendChild('
                'viz.renderSVGElement(el.dataset.dot));}catch(e){el.textContent'
                '="graph: "+e;}});});</script>')
        self.page("graphs.html", "Graphs", "\n".join(body))

    def latex_page(self, files):
        body = ["<h1>LaTeX</h1><p>The reaction table for a paper "
                "(longtable, amsmath): <a href='latex/%s'>%s</a>, to "
                "<code>\\input</code>; <a href='latex/%s'>%s</a>, a document "
                "around it with the definitions.</p>" % (
                    files[0], files[0], files[1], files[1])]
        self.page("latex.html", "LaTeX", "\n".join(body))
