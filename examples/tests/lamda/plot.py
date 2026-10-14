# plot of the lamda test: the O2 cooling per ntot [erg s^-1] over the
# (ntot, Tgas) grid (cooling.dat: ntot, Tgas, cooling; one block per ntot)
import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt("cooling.dat")
n = np.unique(d[:, 0])
T = np.unique(d[:, 1])
c = d[:, 2].reshape(len(n), len(T)).T  # rows Tgas, columns ntot
fig, ax = plt.subplots()
m = ax.pcolormesh(np.log10(n), np.log10(T), np.log10(np.maximum(c, 1e-99)),
                  shading="nearest")
fig.colorbar(m, ax=ax, label="log(O2 cooling / ntot / erg s$^{-1}$)")
ax.set_xlabel("log(ntot / cm$^{-3}$)")
ax.set_ylabel("log(T$_{gas}$ / K)")
ax.set_title("KROME TEST: LAMDA")
fig.savefig("plot.png")
print("wrote plot.png")
