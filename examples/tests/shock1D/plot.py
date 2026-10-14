# plot of the shock1D test: final profiles (hydro_final.dat: hydro,
# chemistry_final.dat: mass fractions, both with a names header)
import matplotlib.pyplot as plt
import numpy as np


def load(fname):
    names = open(fname).readline().lstrip("#").split()
    return names, np.loadtxt(fname)


hn, hd = load("hydro_final.dat")
cn, cd = load("chemistry_final.dat")
r = hd[:, hn.index("radius")]
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
ax[0].loglog(r, hd[:, hn.index("Tgas")])
ax[0].set_ylabel("Tgas/K")
ax[1].loglog(r, hd[:, hn.index("rho")])
ax[1].set_ylabel(r"density/(g cm$^{-3}$)")
ax[2].loglog(cd[:, cn.index("radius")], cd[:, cn.index("H2")])
ax[2].set_ylabel(r"H$_2$ mass fraction")
ax[2].set_xlabel("radius/cm")
ax[2].set_xlim(7e17, 1e19)
ax[0].set_title("KROME TEST: SHOCK 1D")
fig.savefig("plot.png")
print("wrote plot.png")
