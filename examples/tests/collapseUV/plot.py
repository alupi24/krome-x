# plot of the collapseUV test (abundances.dat: J21 ntot Tgas H2 H, one block
# per J21, blocks separated by blank lines)
import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt("abundances.dat")
fig, (a1, a2) = plt.subplots(2, 1, sharex=True)
for j in np.unique(d[:, 0]):
    b = d[d[:, 0] == j]
    a1.loglog(b[:, 1], b[:, 2], label="J21=%g" % j)
    a2.loglog(b[:, 1], b[:, 3])
a1.set_ylabel("Tgas/K")
a1.legend(loc="lower right")
a2.set_ylabel("xH2")
a2.set_ylim(bottom=1e-8)
a2.set_xlabel(r"density/cm$^{-3}$")
fig.savefig("plot.png")
print("wrote plot.png")
