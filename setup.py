"""pip install . (or pip install -e .): the package kromex, and, in a
regular install, the files it generates projects from (runtime/, data/,
networks/, examples/, VERSION, options_template.toml) and the user
documentation (docs/manual.md, docs/network-format.md) copied into
kromex/_share (kromex/paths.py finds them there; an editable install uses
the repository's own). The metadata is in pyproject.toml."""

import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

HERE = Path(__file__).resolve().parent
SHARE = ["runtime", "data", "networks", "examples", "VERSION",
         "options_template.toml", "docs/manual.md", "docs/network-format.md"]
# never shipped: build products, and the ZCIE tables extracted from
# data/coolZ_CIE.tar.xz on first use (77 MB)
SKIP_PARTS = {"__pycache__", ".git"}
SKIP_NAMES = {"coolZ_CIE2005.dat", "coolZ_CIE2012.dat",
              "coolZ_CIE2012NOUV.dat", "KROME_ERROR_REPORT"}


def share_files():
    """the files to ship: the tracked ones in a git checkout, else every
    file of SHARE but the build products"""
    try:
        out = subprocess.run(["git", "ls-files", "-z", "--"] + SHARE,
                             cwd=HERE, capture_output=True, text=True,
                             check=True).stdout
        files = [f for f in out.split("\0") if f]
        if files:
            return files
    except (OSError, subprocess.CalledProcessError):
        pass
    files = []
    for top in SHARE:
        p = HERE / top
        for f in ([p] if p.is_file() else sorted(p.rglob("*"))):
            rel = f.relative_to(HERE)
            if f.is_file() and not (SKIP_PARTS & set(rel.parts)) \
                    and f.name not in SKIP_NAMES and f.suffix not in (
                        ".o", ".a", ".so", ".dylib", ".pyc"):
                files.append(str(rel))
    return files


class BuildWithShare(build_py):
    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return  # the editable install reads the repository itself
        dest = Path(self.build_lib) / "kromex" / "_share"
        for rel in share_files():
            src = HERE / rel
            if src.is_file():
                (dest / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dest / rel)


setup(cmdclass={"build_py": BuildWithShare})
