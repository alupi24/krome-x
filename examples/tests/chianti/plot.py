# plot of the chianti test: the equilibrium cooling of carbon [erg cm^3
# s^-1, n(C) = n(e) = 1] vs Tgas, one curve per J21 (cooling.dat: J21, Tgas,
# cooling; one block per J21)
import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt("cooling.dat")
fig, ax = plt.subplots()
for j21 in np.unique(d[:, 0]):
    b = d[d[:, 0] == j21]
    ax.loglog(b[:, 1], b[:, 2], label="J21=%g" % j21)
ax.set_ylim(1e-21, 1e-17)
ax.set_xlabel("T$_{gas}$ / K")
ax.set_ylabel("cooling / (erg cm$^{3}$ s$^{-1}$)")
ax.set_title("KROME TEST: CHIANTI")
ax.legend()
fig.savefig("plot.png")
print("wrote plot.png")
