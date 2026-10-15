"""The network as graphs, one per element (as legacy docmake): the species
with that element, an edge from each reactant to each product that both
have it, labelled with the other reactants. Graphviz .dot files, rendered
to SVG when Graphviz's `dot` is installed."""

import shutil
import subprocess

MAX_NODES = 60


def element_graphs(net, max_nodes=MAX_NODES):
    """{element: {"nodes": [...], "edges": {(a, b): {partners}}}}, and the
    elements left out (too many species)"""
    sp = net.species
    elements = sorted({e for s in sp.values() if s["composition"]
                       for e in s["composition"]})
    graphs, skipped = {}, []
    for el in elements:
        def has(x):
            c = sp[x]["composition"]
            return bool(c) and el in c
        edges = {}
        for rx in net.reactions:
            reac = [r["sp"] for r in rx["reactants"]]
            prod = [p["sp"] for p in rx["products"]]
            for a in reac:
                if not has(a):
                    continue
                partners = sorted(set(x for x in reac if x != a))
                if rx["rate"]["kind"] == "photo":
                    partners.append("photon")
                for b in prod:
                    if has(b) and b != a:
                        edges.setdefault((a, b), set()).update(partners)
        nodes = sorted({x for e in edges for x in e})
        if not nodes:
            continue
        if len(nodes) > max_nodes:
            skipped.append((el, len(nodes)))
            continue
        graphs[el] = {"nodes": nodes, "edges": edges}
    return graphs, skipped


def dot_text(el, g):
    ids = {x: "n%d" % i for i, x in enumerate(g["nodes"])}
    L = ["digraph %s {" % ("el_" + el), "  rankdir=LR;",
         '  node [shape=box, style=rounded, fontname="Helvetica", '
         'fontsize=11];', '  edge [fontname="Helvetica", fontsize=8, '
         'color="#555555"];']
    for x in g["nodes"]:
        L.append('  %s [label="%s"];' % (ids[x], x))
    for (a, b), partners in sorted(g["edges"].items()):
        lab = ", ".join(sorted(partners))
        L.append('  %s -> %s [label="%s"];' % (ids[a], ids[b], lab))
    L.append("}")
    return "\n".join(L) + "\n"


def write(net, outdir, max_nodes=MAX_NODES):
    """write graphs/<element>.dot (and .svg with Graphviz); returns
    ([(element, svg or None, nodes)], skipped, rendered)"""
    outdir.mkdir(parents=True, exist_ok=True)
    graphs, skipped = element_graphs(net, max_nodes)
    dot = shutil.which("dot")
    done = []
    for el, g in graphs.items():
        p = outdir / ("%s.dot" % el)
        p.write_text(dot_text(el, g))
        svg = None
        if dot:
            q = outdir / ("%s.svg" % el)
            try:
                subprocess.run([dot, "-Tsvg", str(p), "-o", str(q)],
                               check=True, capture_output=True, timeout=120)
                svg = q.name
            except (subprocess.SubprocessError, OSError):
                svg = None
        done.append((el, svg, len(g["nodes"])))
    return done, skipped, dot is not None
