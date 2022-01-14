import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("out.dat")
T, X, Z = [0.5, 1.5], [2.0, 8.0, 151], data

Z = np.flip(Z.reshape((-1, int(X[2]))), 0)

plt.figure(dpi=180)

plt.xlabel("Position (nm)")
plt.ylabel("Time (fs)")
plt.title("Evolution of a free particle")
plt.imshow(Z, cmap='inferno',aspect='auto', extent=(X[0], X[1], T[0], T[1]))
plt.colorbar(label='Probability density (nm${}^{-1}$)')

plt.savefig("figs/fig.png")