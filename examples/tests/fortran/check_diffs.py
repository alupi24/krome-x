"""compare the Fortran-binding run (abundances.dat) with the C++ run
(abundances_cpp.dat). Both call the same solver: built with strict
floating point (make CXXFLAGS="-O2 -std=c++17 -ffp-contract=off" FFLAGS="-O2
-ffp-contract=off") they agree to the printed digits (5e-8); the default
-O3 FMA build compiles the solver twice (libkrome for the binding, the
header-only API for test_cpp) with different contractions, ~1e-6"""
import sys


def load(name):
    with open(name) as f:
        return [[float(v) for v in ln.split()] for ln in f
                if ln.strip() and not ln.startswith("#")]


a, b = load("abundances.dat"), load("abundances_cpp.dat")
dev = max(abs(u - v) / max(abs(v), 1e-300)
          for ra, rb in zip(a, b) for u, v in zip(ra, rb)
          if max(abs(u), abs(v)) > 1e-30)
ok = len(a) == len(b) and dev < 1e-5
print("| abundances.dat, abundances_cpp.dat: %d rows, max relative "
      "difference %.1e: %s"
      % (len(a), dev, "the same" if ok else "DO NOT match!!"))
sys.exit(0 if ok else 1)
