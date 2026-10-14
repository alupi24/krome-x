# plot of the map test (map.dat: rho, J21, TgasInit, Tgas, then the
# mass fractions; one block of densities per J21 value): Tgas map and
# x(H2) contours
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

names = open("map.dat").readline().lstrip("#").split()
d = np.loadtxt("map.dat")
rho, j21 = np.unique(d[:, 0]), np.unique(d[:, 1])
shape = (len(j21), len(rho))  # rows: J21, columns: rho (file order)
R, J = d[:, 0].reshape(shape), d[:, 1].reshape(shape)
T = d[:, names.index("Tgas")].reshape(shape)
xh2 = d[:, names.index("H2")].reshape(shape)
fig, ax = plt.subplots()
m = ax.pcolormesh(R, J, T, norm=LogNorm(), shading="auto")
fig.colorbar(m, label="Tgas/K")
c = ax.contour(R, J, xh2, levels=[1e-5, 1e-3, 1e-2, 1e-1, .9], colors="k")
ax.clabel(c, fmt="%g")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"density/(g cm$^{-3}$)")
ax.set_ylabel("J21")
ax.set_title("Tgas, and x(H2) contours")
fig.savefig("plot.png")
print("wrote plot.png")
