"""kromex_doc: the documentation of a KROME-X network (the successor of
legacy KROME's docmake), from its .knet file alone (no build). A client of
KROME-X (its .knet parser and data), independent of the code generation:

- an HTML site: the species (formation and destruction channels, mass,
  charge, composition, formation enthalpy), the reactions (the rate as a
  formula, its source line, its T range and reference, k against T and
  the other plotted quantities), the piecewise rates and the jumps at
  their bounds, the missing reactions and branches, one graph per
  element;
- the reaction table for a paper (LaTeX: longtable, amsmath).

    python3 -m kromex_doc NETWORK [-o DIR] [--range Av=0.1,20] [--set crflux=1.3e-17]
    kromex-doc NETWORK ...                     (installed with pip)

In a terminal it then serves the site on http://127.0.0.1:8080/ until
Ctrl-C (--no-serve, --port).

    from kromex_doc import document
    document("primordial", "doc_primordial", ranges={"Tgas": (10, 1e5)})

The rate coefficients are evaluated in numpy from the network's
expressions; photo and surface rates, rates from the database, tables and
extern functions are computed by the generated code only, and are
documented without plots."""

import math
import re
import sys
from pathlib import Path

import kromex  # noqa: F401 (puts the kromex modules on the path)
import balance  # noqa: E402
import codegen  # noqa: E402
import knet  # noqa: E402
import paths  # noqa: E402

from . import chem, graph, latex, missing
from .pages import Site

__all__ = ["document"]

TRANGE = (10.0, 1e5)


def document(network, outdir=None, ranges=None, values=None, tref=100.0,
             plots=True, max_graph=graph.MAX_NODES,
             max_missing=missing.MAX_SPECIES, outside="asym", static=False,
             quiet=False):
    """write the documentation of network (a .knet path, or a shipped
    network by name) into outdir (default doc_<name>).

    ranges  {name: (lo, hi)}: Tgas (default 10-1e5 K) and the parameters
            or state quantities to plot the rates against
    values  {name: value}: the reference values of parameters (default:
            the network's defaults), Tdust, Hnuclei, ntot, J_LW, and the
            densities "n(X)" that rates may use
    tref    the temperature of the k column of the tables [K]
    static  SVG plots only (default: interactive plots by plotly.js, with
            the SVG kept where plotly.js does not load)
    outside the rates outside their T range, as the generation option
            rates.outside: "asym" (zero below, the bound's value above;
            the default), "clamp" or "zero"

    Returns the path of index.html."""
    src = paths.find(str(network), "networks", suffix=".knet")
    net = knet.parse_network(src)
    if net.errors:
        raise ValueError("%s: %s" % (src, "; ".join(net.errors[:5])))
    out = Path(outdir or "doc_%s" % src.stem)
    out.mkdir(parents=True, exist_ok=True)
    ranges = dict(ranges or {})
    ranges.setdefault("Tgas", TRANGE)
    ranges = {"Tgas": ranges.pop("Tgas"), **ranges}
    for k, (lo, hi) in ranges.items():
        if not 0 < lo < hi:
            raise ValueError("range %s: need 0 < lo < hi (log axes)" % k)
    values = dict(values or {})
    if outside not in ("asym", "clamp", "zero"):
        raise ValueError("outside: asym, clamp or zero, not %r" % outside)

    def say(msg):
        if not quiet:
            print(msg, flush=True)

    say("%s: %d species, %d reactions" % (src.name, len(net.species),
                                          len(net.reactions)))
    thermo = chem.Thermo()
    say("NOTE: " + codegen.BURCAT_TERMS)
    say("missing chemistry ...")
    miss = missing.analyse(net, thermo, max_species=max_missing)
    say("detailed balance ...")
    lo, hi = ranges["Tgas"]
    rev = balance.analyse(net, balance.load_tables(paths.ROOT / "data"),
                          tgrid=(lo, hi, balance.TGRID[2]))
    say("graphs ...")
    graphs = graph.write(net, out / "graphs", max_graph)
    say("LaTeX table ...")
    (out / "latex").mkdir(exist_ok=True)
    files = ("%s_table.tex" % src.stem, "%s_reactions.tex" % src.stem)
    latex.table(net, out / "latex" / files[0], out / "latex" / files[1],
                "The reactions of %s" % src.name)
    say("pages%s ..." % (" and plots" if plots else ""))
    site = Site(net, src, out, values, ranges, float(tref), thermo,
                max_graph=max_graph, plots=plots, outside=outside,
                interactive=not static)
    site.build(miss, graphs, files, rev)
    say("wrote %s" % (out / "index.html"))
    return out / "index.html"


def make_server(outdir, port=8080, tries=10):
    """an HTTP server of outdir on 127.0.0.1 (this computer only, unlike
    python -m http.server, which listens on every interface): port, or the
    next free one of the following tries; returns (server, url)"""
    import functools
    import http.server
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=str(outdir))
    handler.log_message = lambda *args: None      # quiet
    for p in range(port, port + max(1, tries)):
        try:
            srv = http.server.ThreadingHTTPServer(("127.0.0.1", p), handler)
        except OSError:
            continue
        return srv, "http://127.0.0.1:%d/" % srv.server_address[1]
    raise OSError("no free port in %d-%d" % (port, port + tries - 1))


def serve(outdir, port=8080):
    """serve outdir (index.html) until Ctrl-C"""
    srv, url = make_server(outdir, port)
    print("serving %s at %s (Ctrl-C to stop)" % (outdir, url), flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print()
    finally:
        srv.server_close()


def parse_assign(item, what):
    """NAME=VALUE (or NAME=LO,HI for a range)"""
    m = re.match(r"^\s*([^=\s]+)\s*=\s*(.+)$", item)
    if not m:
        raise ValueError("--%s %s: expected NAME=VALUE" % (what, item))
    name, raw = m.group(1), m.group(2)
    try:
        v = [float(x) for x in raw.split(",")]
    except ValueError:
        raise ValueError("--%s %s: not a number" % (what, item))
    if what == "range":
        if len(v) != 2:
            raise ValueError("--range %s: expected NAME=LO,HI" % item)
        return name, tuple(v)
    if len(v) != 1 or not math.isfinite(v[0]):
        raise ValueError("--set %s: expected NAME=VALUE" % item)
    return name, v[0]


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        prog="kromex-doc",
        description="The documentation of a KROME-X network (HTML site and "
                    "LaTeX table), from its .knet file.")
    ap.add_argument("network", help="a .knet file, or a shipped network "
                    "by name (primordial)")
    ap.add_argument("-o", "--outdir", help="output folder (default "
                    "doc_<network>)")
    ap.add_argument("--range", action="append", default=[],
                    metavar="NAME=LO,HI", help="plot the rates against "
                    "NAME over LO-HI (Tgas: default 10,1e5; others: e.g. "
                    "Av=0.1,20)")
    ap.add_argument("--set", action="append", default=[],
                    metavar="NAME=VALUE", help="a reference value: a "
                    "parameter (crflux=1.3e-17), Tdust, Hnuclei, or a "
                    "density n(X)=1e-4")
    ap.add_argument("--tref", type=float, default=100.0,
                    help="the temperature of the k column [K] (default 100)")
    ap.add_argument("--outside", default="asym",
                    choices=("asym", "clamp", "zero"),
                    help="the rates outside their T range, as the "
                    "generation option rates.outside (default asym)")
    ap.add_argument("--static", action="store_true",
                    help="SVG plots only (default: interactive plots, by "
                    "plotly.js from a CDN, with the SVG where it does not "
                    "load)")
    ap.add_argument("--no-plots", action="store_true",
                    help="no rate plots (faster for large networks)")
    ap.add_argument("--max-graph", type=int, default=graph.MAX_NODES,
                    help="elements with more species get no graph "
                    "(default %d)" % graph.MAX_NODES)
    ap.add_argument("--max-missing", type=int, default=missing.MAX_SPECIES,
                    help="list the missing reactions up to this many "
                    "species per phase (default %d)" % missing.MAX_SPECIES)
    sv = ap.add_mutually_exclusive_group()
    sv.add_argument("--serve", dest="serve", action="store_true",
                    default=None, help="serve the site afterwards on "
                    "http://127.0.0.1:PORT (the default in an interactive "
                    "terminal)")
    sv.add_argument("--no-serve", dest="serve", action="store_false",
                    help="do not serve it (the default in scripts)")
    ap.add_argument("--port", type=int, default=8080,
                    help="the port to serve on, or the next free one "
                    "(default 8080)")
    a = ap.parse_args(argv)
    try:
        ranges = dict(parse_assign(x, "range") for x in a.range)
        values = dict(parse_assign(x, "set") for x in a.set)
        index = document(a.network, a.outdir, ranges, values, a.tref,
                         plots=not a.no_plots, max_graph=a.max_graph,
                         max_missing=a.max_missing, outside=a.outside,
                         static=a.static)
        # serving blocks until Ctrl-C: by default only in a terminal
        if a.serve if a.serve is not None else \
                (sys.stdin.isatty() and sys.stdout.isatty()):
            serve(index.parent, a.port)
    except (ValueError, OSError) as e:
        print("ERROR: %s" % e)
        return 1
    return 0
