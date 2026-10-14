# plot of the cloudSurface test (abundances.dat: time [yr], then n_i / 2e4)
# against the Semenov+2010 (ALCHEMIC) curves of semenov.dat
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
s = np.loadtxt("semenov.dat")  # time HCO+ H2O H2S HCN HCO CO
panels = [("CO", 6, (1e-9, 1e-4)), ("HCO+", 1, (1e-13, 1e-8)),
          ("H2O", 2, (1e-11, 5e-6)), ("H2S", 3, (1e-17, 1e-10))]
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
for a, (sp, col, yr) in zip(ax.flat, panels):
    a.loglog(d[:, 0], d[:, names.index(sp)], "o", ms=3, label=sp + " KROME")
    a.loglog(s[:, 0], s[:, col], label=sp + " Semenov2010")
    a.set_xlim(1e2, 1e8)
    a.set_ylim(*yr)
    a.legend(loc="best")
for a in ax[:, 0]:
    a.set_ylabel("fractional abundances")
for a in ax[1]:
    a.set_xlabel("Time [yr]")
fig.savefig("plot.png")
print("wrote plot.png")
