# plot of the customCooling test (cooling.dat: Tgas, then the cooling terms)
import matplotlib.pyplot as plt
import numpy as np

names = open("cooling.dat").readline().lstrip("#").split()
d = np.loadtxt("cooling.dat", ndmin=2)
plt.loglog(d[:, 0], d[:, names.index("CUSTOM")], label="H2O cooling (CUSTOM)")
plt.xlabel(r"T$_{gas}$/K")
plt.ylabel(r"cooling/(erg cm$^{-3}$ s$^{-1}$)")
plt.title("KROME TEST: custom cooling")
plt.legend()
plt.savefig("plot.png")
print("wrote plot.png")
