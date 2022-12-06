import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

calculation_results = np.genfromtxt("../task1/calculation_results_task1.csv", delimiter=',')
sampled_points = np.genfromtxt("../task2/sampled_points.csv", delimiter=',')
integral_values = calculation_results[:,0]; integral_std = calculation_results[:,1]

number_of_points = len(sampled_points)
point_linspace = np.linspace(0, number_of_points, number_of_points)

n_bins = 100
counts, bins = np.histogram(sampled_points, bins = n_bins, density = True)

fig, ax = plt.subplots(1,1)
ax.stairs(counts, bins, fill=True, label = "sampled points")
#ax.plot(point_linspace, sampled_points, label = sampled_points)
ax.set_title(f"Distribution of {number_of_points:.4g} sinusoidally sampled points")
ax.set_xlabel(f"X-value")
ax.set_ylabel("Number of points")#, fontsize=15)
ax.legend()
plt.tight_layout()
fig.savefig(f'plots/sampled_points.png')