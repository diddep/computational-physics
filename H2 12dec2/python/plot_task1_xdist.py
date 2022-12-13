import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_theme()


#arr_xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
arr_theta =np.genfromtxt("../theta.csv", delimiter=',')
arr_xdist = np.cos(arr_theta)
arr_sin = np.sin(arr_theta)
print(arr_xdist.shape)

n_bins =70
fig_xdist, ax_dist = plt.subplots(1,3, figsize =(15,5))

counts_x, bins_x = np.histogram(arr_xdist, bins = n_bins, density = True)

counts_theta, bins_theta = np.histogram(arr_theta, bins = n_bins, density = True)
counts_sin, bins_sin = np.histogram(arr_sin, bins = n_bins, density = True)


ax_dist[0].stairs(counts_x, bins_x, fill = True, label='Sampled distribution')
ax_dist[0].axhline(1/2, label='Distribution for uncorrelated electrons', color='r', linewidth=2, linestyle='dashed')
ax_dist[0].set_title('Distribution for x', fontsize=15)
ax_dist[0].set_xlabel('x', fontsize=15)
ax_dist[0].set_ylabel('Probability density', fontsize=15)
ax_dist[0].legend()



ax_dist[1].stairs(counts_theta, bins_theta, fill = True, label='Sampled distribution')
ax_dist[1].set_title(r'Distribution for $\theta$', fontsize=15)
ax_dist[1].set_xlabel(r'$\theta$ [rad]', fontsize=15)
ax_dist[1].set_ylabel('Probability Density', fontsize=15)


ax_dist[2].stairs(counts_sin, bins_sin, fill = True, label='Sampled distribution')
ax_dist[2].set_title(r'Distribution for $\sin(\theta)$', fontsize=15)
ax_dist[2].set_xlabel(r'$\sin(\theta)$', fontsize=15)
ax_dist[2].set_ylabel('Probability Density', fontsize=15)
plt.tight_layout()

fig_xdist.savefig('xdist_new.png')