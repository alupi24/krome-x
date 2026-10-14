# plot of the earlyUniverse test (abundances.dat: 1+z, density, Tgas, Trad,
# then n_i / ntot)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
z = d[:, 0]
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
groups = [("H", "H+", "H-", "H2", "H2+", "H3+"),
          ("He", "He+", "He++", "HeH+"),
          ("D", "D+", "HD", "HD+", "H2D+")]
for a, g in zip(ax.flat, groups):
    for s in g:
        a.loglog(z, d[:, names.index(s)], label=s)
    a.set_ylim(1e-25, 1e1)
    a.set_ylabel("n/ntot")
    a.legend(fontsize=7)
a = ax[1, 1]
a.loglog(z, d[:, names.index("Tgas")], label="Tgas")
a.loglog(z, d[:, names.index("Trad")], label="Trad")
a.set_ylabel("T/K")
a.legend()
for a in ax[1]:
    a.set_xlabel("1+z")
ax[0, 0].set_xlim(1e4, 1e1)  # decreasing redshift
fig.savefig("plot.png")
print("wrote plot.png")
