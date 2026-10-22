"""Small log-log line plots as SVG (no matplotlib: KROME-X needs numpy
only, and a site of thousands of reactions stays small)."""

import json
import math
from html import escape

import numpy as np

COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
          "#8c564b", "#e377c2", "#17becf"]
DASH = {"solid": "", "dashed": "6,4", "dotted": "2,3"}
_N = [0]   # unique clip ids for several plots in one page


def _ticks(lo, hi):
    """decade ticks in [lo, hi] (log10 values), with 2 and 5 when the span
    is short"""
    a, b = math.floor(lo), math.ceil(hi)
    step = max(1, int(math.ceil((b - a) / 8)))
    major = [float(e) for e in range(a, b + 1, step) if lo - 1e-9 <= e <= hi + 1e-9]
    minor = []
    if hi - lo < 2.5:
        for e in range(a, b + 1):
            for m in (2, 5):
                v = e + math.log10(m)
                if lo <= v <= hi:
                    minor.append(v)
    return major, minor


def _label(e):
    """10^e as SVG text"""
    if e == 0:
        return "1"
    if e == 1:
        return "10"
    return '10<tspan dy="-6" font-size="9">%d</tspan>' % e


def _prepare(series, xrange, decades):
    """the positive points of every series, their runs (a line breaks where
    the values are not positive), and the axes' log10 ranges; None when
    nothing is positive. A series with "scale": False (context, e.g. a
    formula continued outside its T range) is drawn but does not set the
    y range, unless no other series has positive values"""
    xs, ys, runs, yscale = [], [], [], []
    for s in series:
        x, y = np.asarray(s["x"], float), np.asarray(s["y"], float)
        m = (x > 0) & (y > 0) & np.isfinite(y)
        xs.append(x[m])
        ys.append(y[m])
        if s.get("scale", True):
            yscale.append(y[m])
        idx = np.flatnonzero(m)
        cut = np.flatnonzero(np.diff(idx) > 1) + 1
        runs.append([(x[r], y[r]) for r in np.split(idx, cut) if len(r)])
    if not any(len(y) for y in ys):
        return None
    xall = np.concatenate([x for x in xs if len(x)])
    yall = np.concatenate([y for y in (yscale if any(len(y) for y in yscale)
                                       else ys) if len(y)])
    x0, x1 = (np.log10(xrange[0]), np.log10(xrange[1])) if xrange else \
        (np.log10(xall.min()), np.log10(xall.max()))
    y1 = np.log10(yall.max())
    y0 = max(np.log10(yall.min()), y1 - decades)
    if y1 - y0 < 1.0:
        mid = 0.5 * (y0 + y1)
        y0, y1 = mid - 1.0, mid + 1.0
    else:
        y0, y1 = y0 - 0.05 * (y1 - y0), y1 + 0.05 * (y1 - y0)
    if x1 <= x0:
        x0, x1 = x0 - 1.0, x1 + 1.0
    return xs, ys, runs, (x0, x1, y0, y1)


PLOTLY_DASH = {"solid": "solid", "dashed": "dash", "dotted": "dot"}


def _r(v):
    """a value for the JSON: 5 significant digits, None where the plot
    breaks"""
    return float("%.5g" % v) if v > 0 and math.isfinite(v) else None


def figure(series, xlabel, ylabel="k", title="", xrange=None, note=None,
           decades=20.0, interactive=True):
    """a plot for a page: the SVG, and with interactive the same data for
    plotly.js, which replaces the SVG in the browser when it loads (zoom,
    values on hover); without plotly.js (offline, no JavaScript, print)
    the SVG stays. None when nothing is positive."""
    svg = plot(series, xlabel, ylabel, title, xrange=xrange, note=note,
               decades=decades)
    if svg is None or not interactive:
        return svg
    x0, x1, y0, y1 = _prepare(series, xrange, decades)[3]
    data = []
    for i, s in enumerate(series):
        color = s.get("color", COLORS[i % len(COLORS)])
        x = np.asarray(s["x"], float)
        y = np.broadcast_to(np.asarray(s["y"], float), x.shape)
        data.append({
            "x": [_r(v) for v in x], "y": [_r(v) for v in y],
            "mode": "lines", "name": s.get("label") or "",
            "showlegend": bool(s.get("label")),
            "opacity": s.get("opacity", 1.0),
            "line": {"color": color, "width": s.get("width", 1.6),
                     "dash": PLOTLY_DASH[s.get("style", "solid")]},
            "hovertemplate": "%{x:.4g}, %{y:.4g}<extra></extra>"})
        if s.get("markers"):
            mk = [(a, b) for a, b in s["markers"] if a > 0 and b > 0]
            data.append({"x": [_r(a) for a, _ in mk],
                         "y": [_r(b) for _, b in mk], "mode": "markers",
                         "showlegend": False, "marker": {"color": color},
                         "hovertemplate": "%{x:.4g}, %{y:.4g}<extra></extra>"})
    layout = {
        "height": 360, "margin": {"l": 70, "r": 20, "t": 40 if title else
                                                           20, "b": 50},
        "title": {"text": title} if title else None,
        "xaxis": {"type": "log", "title": {"text": xlabel},
                  "range": [x0, x1], "exponentformat": "power"},
        "yaxis": {"type": "log", "title": {"text": ylabel},
                  "range": [y0, y1], "exponentformat": "power"},
        "legend": {"x": 1, "xanchor": "right", "y": 1},
        "annotations": [{"text": note, "xref": "paper", "yref": "paper",
                         "x": 0.01, "y": 0.02, "showarrow": False,
                         "font": {"size": 10, "color": "#666"}}]
        if note else []}
    fig = json.dumps({"data": data, "layout": layout},
                     separators=(",", ":"))
    return '<div class="kplot" data-fig="%s">%s</div>' % (
        escape(fig, quote=True), svg)


def plot(series, xlabel, ylabel="k", title="", width=560, height=330,
         xrange=None, note=None, decades=20.0):
    """series: [{"x", "y", "style" (solid|dashed|dotted), "color", "width",
    "label", "opacity", "markers"}]; log-log; a line breaks where its
    values are not positive; the y axis spans at most `decades` below the
    largest value. Returns the SVG (a string), or None when nothing is
    positive."""
    L, R, T, B = 64, 14, 26 if title else 12, 42
    pw, ph = width - L - R, height - T - B
    prep = _prepare(series, xrange, decades)
    if prep is None:
        return None
    xs, ys, runs, (x0, x1, y0, y1) = prep

    def px(v):
        return L + (v - x0) / (x1 - x0) * pw

    def py(v):
        return T + (y1 - v) / (y1 - y0) * ph

    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
           'viewBox="0 0 %d %d" font-family="sans-serif" font-size="11">'
           % (width, height, width, height),
           '<rect x="%d" y="%d" width="%d" height="%d" fill="none" '
           'stroke="#888"/>' % (L, T, pw, ph)]
    if title:
        out.append('<text x="%d" y="16" text-anchor="middle" '
                   'font-size="12">%s</text>' % (L + pw / 2, escape(title)))
    for axis, (lo, hi) in (("x", (x0, x1)), ("y", (y0, y1))):
        major, minor = _ticks(lo, hi)
        for v in minor + major:
            big = v in major
            if axis == "x":
                X = px(v)
                out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" '
                           'stroke="#ddd"/>' % (X, T, X, T + ph))
                if big:
                    out.append('<text x="%.1f" y="%d" text-anchor="middle">'
                               '%s</text>' % (X, T + ph + 16,
                                              _label(int(round(v)))))
            else:
                Y = py(v)
                out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" '
                           'stroke="#ddd"/>' % (L, Y, L + pw, Y))
                if big:
                    out.append('<text x="%d" y="%.1f" text-anchor="end">'
                               '%s</text>' % (L - 5, Y + 4,
                                              _label(int(round(v)))))
    out.append('<text x="%d" y="%d" text-anchor="middle">%s</text>'
               % (L + pw / 2, height - 6, escape(xlabel)))
    out.append('<text x="14" y="%d" text-anchor="middle" transform="rotate('
               '-90 14 %d)">%s</text>' % (T + ph / 2, T + ph / 2,
                                          escape(ylabel)))
    _N[0] += 1
    out.append('<clipPath id="c%d"><rect x="%d" y="%d" width="%d" '
               'height="%d"/></clipPath><g clip-path="url(#c%d)">'
               % (_N[0], L, T, pw, ph, _N[0]))
    legend = []
    for i, (s, x, y) in enumerate(zip(series, xs, ys)):
        color = s.get("color", COLORS[i % len(COLORS)])
        for rx, ry in runs[i]:
            # below the axis: drawn at its bottom edge, clipped
            ly = np.maximum(np.log10(ry), y0 - 0.5)
            pts = " ".join("%.1f,%.1f" % (px(a), py(b))
                           for a, b in zip(np.log10(rx), ly))
            out.append('<polyline points="%s" fill="none" stroke="%s" '
                       'stroke-width="%s" stroke-dasharray="%s" '
                       'opacity="%s"/>' % (pts, color, s.get("width", 1.6),
                                           DASH[s.get("style", "solid")],
                                           s.get("opacity", 1.0)))
        if len(x):
            for a, b in s.get("markers", []):
                if a > 0 and b > 0:
                    out.append('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>'
                               % (px(np.log10(a)), py(np.log10(b)), color))
        if s.get("label"):
            legend.append((s["label"], color, s.get("style", "solid"),
                           s.get("opacity", 1.0)))
    out.append("</g>")
    for j, (lab, color, style, op) in enumerate(legend):
        Y = T + 12 + 14 * j
        out.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
                   'stroke-width="2" stroke-dasharray="%s" opacity="%s"/>'
                   '<text x="%d" y="%d" font-size="10">%s</text>'
                   % (L + pw - 150, Y, L + pw - 128, Y, color, DASH[style],
                      op, L + pw - 124, Y + 4, escape(lab)))
    if note:
        out.append('<text x="%d" y="%d" font-size="10" fill="#666">%s</text>'
                   % (L + 6, T + ph - 6, escape(note)))
    out.append("</svg>")
    return "\n".join(out)
