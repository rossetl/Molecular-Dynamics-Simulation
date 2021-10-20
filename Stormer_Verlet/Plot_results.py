import numpy as np
import matplotlib.pyplot as plt

results = np.loadtxt('Result.dat')
fig, ax = plt.subplots(figsize=(12, 6), nrows=1, ncols=2)

ax[0].plot(results[:, 0], results[:, 1], c='coral', lw=2)
ax[0].set_title('Temperature', size=15, fontweight='bold')
ax[0].set_xlabel('time steps', size=15)
ax[0].set_ylabel('temperature', size=15)

ax[1].plot(results[:, 0], results[:, 2], c='dodgerblue', lw=2)
ax[1].set_title('Energy', size=15, fontweight='bold')
ax[1].set_xlabel('time steps', size=15)
ax[1].set_ylabel('energy', size=15)

plt.savefig('energy_temperature.png')
plt.close()

pcf = np.loadtxt('ResultPCF.dat')
plt.figure(figsize=(10,6))
plt.plot(pcf[:,0], pcf[:, 1], c='forestgreen', lw=2)
plt.title('PCA', size=15, fontweight='bold')
plt.xlabel(r'$r/\sigma$', size=15)
plt.ylabel(r'$g(r)$', size=15)
plt.savefig('PCA.png')
