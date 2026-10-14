# plot of the shock1Dcool test: final profiles with cooling (hydro_final.dat)
# against the legacy shock1D run without cooling (shock1D.dat); both files
# have a names header
import matplotlib.pyplot as plt
import numpy as np


def load(fname):
    names = open(fname).readline().lstrip("#").split()
    return names, np.loadtxt(fname)


cn, cd = load("hydro_final.dat")
nn, nd = load("shock1D.dat")
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
for i, (col, lab, yr) in enumerate([("rho", r"density/(g cm$^{-3}$)",
                                     (1e-25, 1e-23)),
                                    ("Tgas", "T/K", (5.0, 4e5))]):
    ax[i].semilogy(nd[:, nn.index("radius")], nd[:, nn.index(col)],
                   label="w/out cooling")
    ax[i].semilogy(cd[:, cn.index("radius")], cd[:, cn.index(col)],
                   label="with cooling")
    ax[i].set_ylabel(lab)
    ax[i].set_ylim(*yr)
    ax[i].legend()
ax[1].set_xlabel("radius/cm")
fig.savefig("plot.png")
print("wrote plot.png")
