# plot of the euler1D test: final Tgas (profiles_final.dat, with a names
# header) against the reference profile of Grassi+2014 (sod_chem.dat)
import matplotlib.pyplot as plt
import numpy as np

names = open("profiles_final.dat").readline().lstrip("#").split()
d = np.loadtxt("profiles_final.dat")
ref = np.loadtxt("sod_chem.dat")
plt.plot(d[:, names.index("x")], d[:, names.index("Tgas")], ".", ms=2,
         label="this test")
plt.plot(ref[:, 1] / 3.08567758128e18, ref[:, 5], label="Grassi+2014")
plt.xlabel("x / pc")
plt.ylabel("Tgas / K")
plt.grid()
plt.legend()
plt.savefig("plot.png")
print("wrote plot.png")
