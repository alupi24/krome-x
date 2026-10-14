"""TOML for the KROME-X configuration (docs/config-design.md section 5):
reading with the standard library's tomllib (Python 3.11+) or the bundled
tomli (kromex/_vendor, Python 3.8-3.10); writing with the small writer
below, which needs only strings, numbers, booleans and lists."""

import math

try:
    import tomllib as _toml
except ImportError:  # Python < 3.11
    from _vendor import tomli as _toml

TOMLDecodeError = _toml.TOMLDecodeError


def loads(text):
    return _toml.loads(text)


def load(path):
    with open(path, "rb") as f:
        return _toml.load(f)


def value(v):
    """one value in TOML syntax"""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if not math.isfinite(v):
            raise ValueError("no TOML for %r" % v)
        s = repr(v)
        return s if any(c in s for c in ".en") else s + ".0"
    if isinstance(v, str):
        return '"%s"' % (v.replace("\\", "\\\\").replace('"', '\\"')
                         .replace("\n", "\\n").replace("\t", "\\t"))
    if isinstance(v, (list, tuple)):
        return "[%s]" % ", ".join(value(x) for x in v)
    raise ValueError("no TOML for %r (%s)" % (v, type(v).__name__))


def dumps(sections, comments=None):
    """[section] tables of key = value lines; a None value is written as a
    comment (TOML has no null: the key is left at its default).
    comments: {section: text} written above each table"""
    out = []
    for sec, keys in sections.items():
        if out:
            out.append("")
        if comments and comments.get(sec):
            out += ["# " + ln if ln else "#"
                    for ln in comments[sec].splitlines()]
        out.append("[%s]" % sec)
        for k, v in keys.items():
            out.append("# %s: default" % k if v is None
                       else "%s = %s" % (k, value(v)))
    return "\n".join(out) + "\n"
