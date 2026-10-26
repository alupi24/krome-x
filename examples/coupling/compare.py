"""compare the outputs of the coupling examples (abundances.dat of
fortran/, c/, cpp/, gpu/fortran/ and gpu/cuda/, after running each in its
folder): the same cell with the same steps. The C and Fortran examples
call the same library (libkrome) and agree to the printed digits; the C++
example compiles the solver itself (the header-only API), ~1e-6 with
different floating-point contractions; the GPU examples integrate on the
device (its math library, and the energy round trip and consistent
advection of each step): 1e-4."""
import sys
from pathlib import Path

HERE = Path(__file__).parent


def load(path):
    with open(path) as f:
        return [[float(v) for v in ln.split()] for ln in f
                if ln.strip() and not ln.startswith("#")]


runs = {d: load(HERE / d / "abundances.dat")
        for d in ("fortran", "c", "cpp", "gpu/fortran", "gpu/cuda")
        if (HERE / d / "abundances.dat").is_file()}
if len(runs) < 2:
    sys.exit("run at least two of the examples first")
ref = "c" if "c" in runs else next(iter(runs))
ok = True
for d, rows in runs.items():
    if d == ref:
        continue
    dev = max(abs(u - v) / max(abs(v), 1e-300)
              for ra, rb in zip(rows, runs[ref]) for u, v in zip(ra, rb)
              if max(abs(u), abs(v)) > 1e-30)
    same = len(rows) == len(runs[ref]) \
        and dev < (1e-4 if d.startswith("gpu") else 1e-5)
    ok = ok and same
    print("%-11s vs %-8s %d rows, max relative difference %.1e: %s"
          % (d, ref, len(rows), dev, "the same" if same else "DIFFERENT"))
sys.exit(0 if ok else 1)
