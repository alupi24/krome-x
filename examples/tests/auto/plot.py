# plot of the auto test (abundances.dat: j21, time/yr, then the species; one
# block per j21 value)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
for j21 in np.unique(d[:, 0]):
    b = d[d[:, 0] == j21]
    for s, ls in (("C", "-"), ("C+", "--")):
        plt.loglog(b[:, 1], b[:, names.index(s)], ls,
                   label="%s, J21=%g" % (s, j21))
plt.ylim(1e-2, 2)
plt.xlabel("time/yr")
plt.ylabel("abundance/cm$^{-3}$")
plt.legend()
plt.savefig("plot.png")
print("wrote plot.png")
