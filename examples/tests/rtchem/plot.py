# plot of the rtchem test: final profiles (profiles_final.dat: time, Av, Tgas,
# then the species) and the attenuated intensity (opacity.dat: Av, E, J)
import matplotlib.pyplot as plt
import numpy as np

names = open("profiles_final.dat").readline().lstrip("#").split()
d = np.loadtxt("profiles_final.dat")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
for s in ("H+", "H", "C", "C+", "CO", "H2", "O"):
    a1.loglog(d[:, 1], d[:, names.index(s)], label=s)
a1.set_xlabel("Av")
a1.set_ylabel("n/cm$^{-3}$")
a1.set_title("chemistry")
a1.legend()
o = np.loadtxt("opacity.dat")
av, e = np.unique(o[:, 0]), np.unique(o[:, 1])
J = o[:, 2].reshape(len(av), len(e))
m = a2.pcolormesh(av, e, np.log10(J.T + 1e-10), shading="auto")
a2.set_xscale("log")
a2.set_xlabel("Av")
a2.set_ylabel("energy/eV")
a2.set_title("log intensity [eV/cm2/s/Hz/sr]")
fig.colorbar(m, ax=a2)
fig.tight_layout()
fig.savefig("plot.png")
print("wrote plot.png")
