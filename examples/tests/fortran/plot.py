# plot of the collapse test (abundances.dat: ntot, Tgas, then n_i / ntot)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
fig, (a1, a2) = plt.subplots(2, 1, sharex=True)
for s in ("H", "H2"):
    a1.loglog(d[:, 0], d[:, names.index(s)], label=s)
a1.set_ylabel("fraction")
a1.set_ylim(top=3)
a1.legend()
a2.loglog(d[:, 0], d[:, 1])
a2.set_ylim(bottom=1e2)
a2.set_ylabel("T/K")
a2.set_xlabel(r"density/cm$^{-3}$")
fig.savefig("plot.png")
print("wrote plot.png")
