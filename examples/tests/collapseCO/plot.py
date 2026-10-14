# plot of the collapseCO test (abundances.dat: ntot, Tgas, Tdust, then n_i /
# ntot)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
fig, ax = plt.subplots()
ax.loglog(d[:, 0], d[:, 1], lw=2, label="Tgas/K")
ax.loglog(d[:, 0], d[:, 2], lw=2, label="Tdust/K")
ax.set_xlabel(r"density/cm$^{-3}$")
ax.set_ylabel("T/K")
ax.grid(True)
a2 = ax.twinx()
for s in ("CO", "CO_ice", "H2O", "H2O_ice"):
    a2.loglog(d[:, 0], d[:, names.index(s)], ls="--", label=s)
a2.set_ylabel("n/ntot")
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = a2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize="small")
fig.savefig("plot.png")
print("wrote plot.png")
