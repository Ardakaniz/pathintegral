import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("out.dat")
T, X, Z = data[:2], data[2:4], data[4:]

T_min, T_max = T[0], T[1]
X_min, X_max = -X[0], X[0]

Z = np.flip(Z.reshape((-1, int(X[1]))), 0)

plt.figure(dpi=180)

plt.xlabel("Position (nm)")
plt.ylabel("Time (fs)")
plt.title("Evolution of a free particle")
plt.imshow(Z, cmap='inferno',aspect='auto', extent=(-X[0], X[0], T[0], T[1]))
plt.colorbar(label='Probability density (nm${}^{-1}$)')

plt.savefig("figs/fig.png")