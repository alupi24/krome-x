# plot of the dust test: gas-phase C and Si (gas.dat) and the dust mass
# density of each grain type (dust.dat, summed over its bins) vs Tgas
import matplotlib.pyplot as plt
import numpy as np

gas = np.loadtxt("gas.dat")
dust = np.loadtxt("dust.dat")
T = gas[:, 0]
fig, ax = plt.subplots()
ax.loglog(T, gas[:, 1], label="gas C")
ax.loglog(T, gas[:, 2], label="gas Si")
for typ, lab in ((0, "dust C"), (1, "dust Si")):
    d = dust[dust[:, 0] == typ]
    tot = {}
    for t, m in zip(d[:, 1], d[:, 2]):
        tot[t] = tot.get(t, 0.0) + m
    ts = sorted(tot)
    ax.loglog(ts, [max(tot[t], 1e-40) for t in ts], "--", label=lab)
ax.set_xlabel(r"T$_{gas}$/K")
ax.set_ylabel(r"$\rho$/(g cm$^{-3}$)")
ax.set_ylim(1e-25, None)
ax.set_title("KROME TEST: dust")
ax.legend()
fig.savefig("plot.png")
print("wrote plot.png")
