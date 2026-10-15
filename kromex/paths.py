"""Where KROME-X finds its files (runtime/, data/, networks/, examples/,
VERSION): the KROME-X folder in a checkout or an editable install, or
kromex/_share in an installed package (setup.py copies them there)."""

from pathlib import Path

HERE = Path(__file__).resolve().parent


def _root():
    for cand in (HERE.parent, HERE / "_share"):
        if (cand / "runtime" / "krome.hpp").is_file() and \
                (cand / "data" / "database").is_dir():
            return cand
    raise RuntimeError("KROME-X: runtime/ and data/ not found next to %s"
                       % HERE)


ROOT = _root()


def find(value, *subdirs, suffix=""):
    """an input file given by the user: as given (absolute, or relative to
    the working directory), else relative to the KROME-X folder or its
    subdirs, with or without suffix (a shipped network by name:
    find("primordial", "networks", suffix=".knet"))"""
    p = Path(value)
    # a file as given; not a folder of the same name (a project generated
    # in the working directory, ./primordial, is not the network)
    if p.is_absolute() or p.is_file():
        return p
    for base in (ROOT,) + tuple(ROOT / s for s in subdirs):
        for q in (base / value,) + ((base / (str(value) + suffix),)
                                     if suffix else ()):
            if q.is_file():
                return q
    return p
