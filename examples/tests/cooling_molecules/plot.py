# plot of the cooling_molecules test (data.out: time [yr], Tgas, then
# n_i / ntot)
import matplotlib.pyplot as plt
import numpy as np

names = open("data.out").readline().lstrip("#").split()
d = np.loadtxt("data.out")
fig, a1 = plt.subplots()
a1.loglog(d[:, 0], d[:, 1], "k", lw=2, label="Tgas")
a1.set_xlabel("time / yr")
a1.set_ylabel("T / K")
a2 = a1.twinx()
for s in ("CO", "OH", "H2O", "HCN"):
    a2.loglog(d[:, 0], d[:, names.index(s)], label=s)
a2.set_ylabel("fraction")
h1, l1 = a1.get_legend_handles_labels()
h2, l2 = a2.get_legend_handles_labels()
a1.legend(h1 + h2, l1 + l2, loc="lower left")
fig.savefig("plot.png")
print("wrote plot.png")
