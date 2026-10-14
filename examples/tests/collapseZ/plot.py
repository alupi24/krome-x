# plot of the collapseZ test (abundances.dat: Z, ntot, Tgas, then n_i / ntot;
# one block per metallicity, blocks separated by blank lines)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
fig, ax = plt.subplots()
for z in np.unique(d[:, 0]):
    b = d[d[:, 0] == z]
    ax.loglog(b[:, 1], b[:, 2], label="Z=-inf" if z < -5 else "Z=%g" % z)
ax.set_xlabel(r"density/cm$^{-3}$")
ax.set_ylabel("Tgas/K")
ax.grid(True)
ax.legend(loc="lower right")
fig.savefig("plot.png")
print("wrote plot.png")
