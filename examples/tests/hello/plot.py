# plot of the hello test (abundances.dat: time, then the species)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
d = np.loadtxt("abundances.dat")
for s in ("FK1", "FK2", "FK3"):
    plt.plot(d[:, 0], d[:, names.index(s)], label=s)
plt.xlabel("time [arbitrary]")
plt.ylabel("abundance [arbitrary]")
plt.legend()
plt.savefig("plot.png")
print("wrote plot.png")
