# plot of the slowmanifold test (abundances.dat: time, then the species; the
# trajectories are separated by blank lines; initial.dat: initial conditions)
import matplotlib.pyplot as plt
import numpy as np

names = open("abundances.dat").readline().lstrip("#").split()
blocks, cur = [], []
for ln in open("abundances.dat"):
    if ln.startswith("#"):
        continue
    if ln.strip():
        cur.append([float(v) for v in ln.split()])
    elif cur:
        blocks.append(np.array(cur))
        cur = []
if cur:
    blocks.append(np.array(cur))
d0 = np.loadtxt("initial.dat", ndmin=2)
iw, i2, io = names.index("H2O"), names.index("H2"), names.index("O")

fig = plt.figure(figsize=(11, 5))
a1 = fig.add_subplot(1, 2, 1)
for b in blocks:
    a1.plot(b[:, iw], b[:, i2], "b-", lw=0.7)
a1.plot(d0[:, iw], d0[:, i2], "kx", ms=4, label="initial conditions")
a1.set_xlabel("H2O")
a1.set_ylabel("H2")
a1.legend()
a1.grid()
a2 = fig.add_subplot(1, 2, 2, projection="3d")
for b in blocks:
    m = b[:, 0] > 3e-4
    a2.plot(b[m, iw], b[m, i2], b[m, io], "b-", lw=0.7)
a2.set_xlabel("H2O")
a2.set_ylabel("H2")
a2.set_zlabel("O")
a2.invert_xaxis()
a2.invert_yaxis()
a2.view_init(20, 70)
fig.savefig("plot.png")
print("wrote plot.png")
