# plot of the collapseUV_Xrays test (abundances.dat: J21 ntot Tgas, one block
# per UV J21, blocks separated by blank lines; J21_X = 0.1)
import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt("abundances.dat")
fig, ax = plt.subplots()
for j in np.unique(d[:, 0]):
    b = d[d[:, 0] == j]
    ax.loglog(b[:, 1], b[:, 2], label="J21x=0.1, J21=%g" % j)
ax.set_xlim(right=1e8)
ax.set_ylim(7e1, 2e4)
ax.set_xlabel(r"density/cm$^{-3}$")
ax.set_ylabel("Tgas/K")
ax.legend(loc="lower right")
fig.savefig("plot.png")
print("wrote plot.png")
