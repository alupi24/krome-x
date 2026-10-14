# plot of the reverse test (abundances.dat: time, then the species in mol)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
for s in ("N2", "O", "NO", "N", "O2"):
    plt.loglog(d[:, 0], d[:, names.index(s)], label=s)
plt.ylim(bottom=1e-6)
plt.xlabel("time/s")
plt.ylabel("concentration/mol")
plt.legend()
plt.savefig("plot.png")
print("wrote plot.png")
