
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_theme()

import time
start_time = time.time()

print("Started program")
print("--- %s seconds ---" % (time.time() - start_time))


arr_R1 =np.genfromtxt("../R1.csv", delimiter=',')
arr_R2 =np.genfromtxt("../R2.csv", delimiter=',')
print("Read position arrays")
print("--- %s seconds ---" % (time.time() - start_time))

#print(arr_R2.shape)
arr_E_L = np.genfromtxt("../E_local.csv", delimiter=',')
arr_E_L_Derivative = np.genfromtxt("../E_local_derivative.csv", delimiter=',')
arr_xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
arr_theta_csv = np.genfromtxt("../theta.csv", delimiter=',')
arr_theta = arr_theta_csv[:,1]

params = np.genfromtxt("../alpha_params.csv", delimiter=',')

print("Read arrays")
print("--- %s seconds ---" % (time.time() - start_time))


n_alpha_steps = int(params[0])
n_timesteps = int(params[3])
d_displacement = params[4]

ix = np.linspace(0, n_timesteps, n_timesteps)


arr_steps = arr_theta_csv[:,0]
#print(arr_E_L.shape)
#print(ix.shape)

#print(np.max(arr_xdist),'xd')

arr_r1 = np.square(arr_R1)
arr_r1 = np.sqrt(np.sum(arr_r1, axis=1))
arr_r2 = np.square(arr_R2)
arr_r2 = np.sqrt(np.sum(arr_r2, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho

print("Start plotting")
print("--- %s seconds ---" % (time.time() - start_time))

n_bins = 70
fig_dist, ax_dist = plt.subplots(1,2, figsize=(10,5))
counts_r1, bins_r1 = np.histogram(arr_r1, bins = n_bins, density = True)
counts_r2, bins_r2 = np.histogram(arr_r2, bins = n_bins, density = True)
rvec = np.linspace(0.1, np.max(arr_r2))
#print(rvec)

ax_dist[0].stairs(counts_r1, bins_r1, fill=True, alpha=0.8)
ax_dist[0].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--',label= r'$\rho $ optimized', linewidth=3)
ax_dist[0].plot(rvec, rho(rvec, 2), color='k', linestyle=':',label=r'$\rho $ unscreened', linewidth=3)
ax_dist[0].set_title(r'Distribution for $r_1 $', fontsize=15)
ax_dist[0].set_xlabel(r"Radius [$a_0 $]", fontsize=15)
ax_dist[0].set_ylabel("Probability Density", fontsize=15)

ax_dist[1].stairs(counts_r2, bins_r2, fill=True, alpha=0.8)
ax_dist[1].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--', label=r'$\rho $ optimized', linewidth=3)
ax_dist[1].plot(rvec, rho(rvec, 2), color='k', linestyle=':', label=r'$\rho $ unscreened', linewidth=3)
ax_dist[1].set_xlabel(r"Radius [$a_0 $]", fontsize=15)
ax_dist[1].set_ylabel("Probability Density", fontsize=15)
ax_dist[1].set_title(r'Distribution for $r_2 $', fontsize=15)
ax_dist[1].legend(fontsize =12)
plt.tight_layout()
fig_dist.savefig('plots_python/hist.png')

print("Plotted histogram")
print("--- %s seconds ---" % (time.time() - start_time))
fig_energy, ax_energy = plt.subplots(1,1)
#ax_energy.scatter(ix, arr_E_L_Derivative)
#fig_energy.savefig('plots_python/energy_r1.png')


#for X in x:
#    print(X)

fig_dist, ax_dist = plt.subplots(1,1)
#counts_theta, bins_theta = np.histogram(x, bins = 10, density = True)
#ax_theta.stairs(counts_theta, bins_theta)
counts_xdist, bins_xdist = np.histogram(arr_xdist, bins = n_bins, density = True)
ax_dist.stairs(counts_xdist, bins_xdist)

fig_dist.savefig('plots_python/x_distribution.png')

print("Plotted x-distribution")
print("--- %s seconds ---" % (time.time() - start_time))

fig_theta, ax_theta = plt.subplots(1,2, figsize=(10,5))
#counts_theta, bins_theta = np.histogram(x, bins = 10, density = True)
#ax_theta.stairs(counts_theta, bins_theta)
counts_x, bins_x = np.histogram(arr_theta, bins = n_bins, density = True)

counts_theta, bins_theta = np.histogram(np.arccos(arr_theta), bins = n_bins, density = True)
ax_theta[0].stairs(counts_x, bins_x, fill = True, label='Sampled distribution')
ax_theta[0].axhline(1/2, label='Distribution for uncorrelated electrons', color='r', linewidth=2, linestyle='dashed')
ax_theta[0].set_title('Distribution for x', fontsize=15)
ax_theta[0].set_xlabel('x', fontsize=15)
ax_theta[0].set_ylabel('Probability density', fontsize=15)
ax_theta[0].legend()

Fill = False
ax_theta[1].stairs(counts_theta, bins_theta, fill = True)
ax_theta[1].set_title(r'Distribution for $\theta$', fontsize=15)
ax_theta[1].set_xlabel(r'$\theta$ [rad]', fontsize=15)
ax_theta[1].set_ylabel('Probability Density', fontsize=15)
plt.tight_layout()
fig_theta.savefig('plots_python/theta.png')

print("Plotted theta and finished program")
print("--- %s seconds ---" % (time.time() - start_time))




