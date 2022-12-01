
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme()


arr_R1 =np.genfromtxt("../R1.csv", delimiter=',')
arr_R2 =np.genfromtxt("../R2.csv", delimiter=',')

#print(arr_R2.shape)
arr_E_L = np.genfromtxt("../E_L.csv", delimiter=',')
arr_xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
arr_theta = np.genfromtxt("../theta.csv", delimiter=',')
print(arr_theta.shape)
print(np.max(arr_xdist),'xd')

arr_r1 = np.square(arr_R1)
arr_r1 = np.sqrt(np.sum(arr_r1, axis=1))
arr_r2 = np.square(arr_R2)
arr_r2 = np.sqrt(np.sum(arr_r2, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho


n_bins = 100
fig_dist, ax_dist = plt.subplots(1,2)
counts_r1, bins_r1 = np.histogram(arr_r1, bins = n_bins, density = True)
counts_r2, bins_r2 = np.histogram(arr_r2, bins = n_bins, density = True)
rvec = np.linspace(0.1, np.max(arr_r2))
#print(rvec)

ax_dist[0].stairs(counts_r1, bins_r1)
ax_dist[0].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--',label= r'$\rho $ optimized')
ax_dist[0].plot(rvec, rho(rvec, 2), color='k', linestyle=':',label=r'$\rho $ unscreened')
ax_dist[0].set_title(r'Distribution for $\mathcal{R}_1 $')

ax_dist[0].set_xlabel(r"Radius [$a_0 $]")
ax_dist[0].set_ylabel("Density [arb. Units]")
ax_dist[1].stairs(counts_r2, bins_r2)
ax_dist[1].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--', label=r'$\rho $ optimized')
ax_dist[0].plot(rvec, rho(rvec, 2), color='k', linestyle=':', label=r'$\rho $ unscreened')
ax_dist[1].set_xlabel(r"Radius [$a_0 $]")
ax_dist[1].set_ylabel("Density [arb. Units]")
ax_dist[1].set_title(r'Distribution for $\mathcal{R}_2 $')
ax_dist[1].legend()
plt.tight_layout()
fig_dist.savefig('plots_python/hist.png')


fig_energy, ax_energy = plt.subplots(1,1)
ax_energy.scatter(arr_r1, arr_E_L)
fig_energy.savefig('plots_python/energy_r1.png')


#for X in x:
#    print(X)

fig_dist, ax_dist = plt.subplots(1,1)
#counts_theta, bins_theta = np.histogram(x, bins = 10, density = True)
#ax_theta.stairs(counts_theta, bins_theta)
counts_xdist, bins_xdist = np.histogram(arr_xdist, bins = n_bins, density = True)
ax_dist.stairs(counts_xdist, bins_xdist)

fig_dist.savefig('plots_python/x_distribution.png')

fig_theta, ax_theta = plt.subplots(1,1)
#counts_theta, bins_theta = np.histogram(x, bins = 10, density = True)
#ax_theta.stairs(counts_theta, bins_theta)
counts_theta, bins_theta = np.histogram(arr_theta[:,1], bins = n_bins, density = True)
ax_theta.stairs(counts_theta, bins_theta)

fig_theta.savefig('plots_python/theta.png')




