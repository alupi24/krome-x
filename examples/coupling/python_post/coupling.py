"""KROME-X coupling example: Python, for post-processing a simulation
snapshot or prototyping a coupling (numpy arrays of cells).

The "snapshot" is a set of numpy arrays over the cells (gas density,
specific internal energy, species mass fractions X[ncell, nspec]), as read
from an HDF5/npz output. Every call works on all cells at once:
net.solve(n[ncell, nspec], Tgas[ncell], dt, ...) integrates them in one
krome_solve_cells call, and the queries (net.mu, net.gamma, net.hnuclei)
are vectorized the same way.

The same operator-split problem as the Fortran, C and C++ examples:
  1. "hydro": a uniform compression rho -> 1.02 rho, adiabatic heating
     e -> e 1.02^(gamma-1)
  2. chemistry + cooling: Tgas from e with KROME's mu and gamma, net.solve
     over dt, then e from the new Tgas with the new mu and gamma
Failed cells (res.flags != 0, left untouched by KROME) are retried with two,
then four substeps. Then, as a post-processing step, the cooling budget of
the final snapshot (net.thermal: the rate of every configured process).

Output: result.dat, per cell: index, initial nH [cm^-3], final Tgas [K],
n(H2)/nH, n(e-)/nH; and cooling.dat (per cell the total cooling [erg cm^-3
s^-1] and each process's share). Run:
    python3 coupling.py path/to/project/libkrome.so   (or .dylib)
(the project's krome.py is imported from the library's directory)
"""
import os
import sys

import numpy as np

LIB = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else None
sys.path.insert(0, os.path.dirname(LIB) if LIB else ".")
import krome  # noqa: E402  (the project's krome.py)

KB, MP = 1.380649e-16, 1.67262192e-24
NCELL, NSTEP, DT, FCOMP = 16, 20, 3e11, 1.02

net = krome.load(LIB)
opts = net.options(evolve_T=True)  # solver options, once
m = net.masses                     # species masses [g]


def temperature(n, e):
    """Tgas [K] from the specific internal energy e; gamma may depend on
    Tgas (codegen thermal.gamma): two fixed-point passes from 1e4 K"""
    T = np.full(len(e), 1e4)
    for _ in range(2):
        T = e * (net.gamma(n, T) - 1.0) * net.mu(n) * MP / KB
    return T


def energy(n, T):
    return KB * T / ((net.gamma(n, T) - 1.0) * net.mu(n) * MP)


def subcycle(n, T):
    """one failed cell: 2, then 4 substeps of DT"""
    for nsub in (2, 4):
        ns, Ts = n.copy(), T
        for _ in range(nsub):
            r = net.solve(ns, Ts, DT / nsub, opts=opts, check=False)
            if r.ierr:
                break
            ns, Ts = r.n, r.Tgas
        else:
            return ns, Ts
    raise krome.KromeError("cell failed after subcycling")


# the snapshot: nH from 1 to 1e4 cm^-3, 1e4 K, 1% ionized
nH0 = 10.0 ** (4.0 * np.arange(NCELL) / (NCELL - 1))
n = np.outer(nH0, np.full(net.nspec, 1e-20))
for s, f in (("H", 0.99), ("H+", 0.01), ("H2", 1e-6), ("He", 0.0775),
             ("e-", 0.01)):
    n[:, net.idx(s)] = f * nH0
rho = net.rho(n)
X = n * m / rho[:, None]
e = energy(n, np.full(NCELL, 1e4))

for step in range(NSTEP):
    # 1. hydro: uniform compression
    n = rho[:, None] * X / m
    g = net.gamma(n, temperature(n, e))
    rho *= FCOMP
    e *= FCOMP ** (g - 1.0)
    # 2. chemistry and cooling, every cell in one call
    n = rho[:, None] * X / m
    res = net.solve(n, temperature(n, e), DT, opts=opts, check=False)
    n, T = res.n, res.Tgas
    for c in np.flatnonzero(res.flags):
        n[c], T[c] = subcycle(n[c], T[c])
    X = n * m / rho[:, None]
    e = energy(n, T)

n = rho[:, None] * X / m
T = temperature(n, e)
nH = net.hnuclei(n)
with open("result.dat", "w") as f:
    for c in range(NCELL):
        f.write("%4d %19.11E %19.11E %19.11E %19.11E\n"
                % (c, nH0[c], T[c], n[c, net.idx("H2")] / nH[c],
                   n[c, net.idx("e-")] / nH[c]))

# post-processing: the cooling budget of the final snapshot, by process
th = net.thermal(n, T)
names = [k for k in th.cooling]
with open("cooling.dat", "w") as f:
    f.write("# cell Tgas cool " + " ".join(names) + "\n")
    for c in range(NCELL):
        f.write("%4d %12.5E %12.5E" % (c, T[c], th.cool[c])
                + "".join(" %9.2E" % (th.cooling[k][c] / th.cool[c])
                          for k in names) + "\n")
print("wrote result.dat, cooling.dat")
