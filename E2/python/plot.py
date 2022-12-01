import numpy as np
import matplotlib.pyplot as plt

energy_array = np.genfromtxt('../verlet_energies.csv', delimiter=',')
position_array = np.genfromtxt('../verlet_positions.csv', delimiter=',')
velocitie_array = np.genfromtxt('../verlet_velocities.csv', delimiter=',')
result_array = np.genfromtxt('../verlet_results.csv', delimiter=',')
param_array = np.genfromtxt('../verlet_params.csv', delimiter=',')
energy_average_array = np.genfromtxt('../verlet_energy_averages.csv', delimiter=',')


t = result_array[:,0]

dt = param_array[0]
alpha = param_array[1]

if(len(t) < 3000):
    str = f"short_alpha_{alpha}"
else:
    str = f"long_alpha_{alpha}"

figp , axp = plt.subplots(1,1)
pos_labels = [f"q({i})" for i in range(5)]

axp.plot(t, position_array[:,0:5], label = pos_labels)

axp.set_xlabel('Time (arb. unit)', fontsize = 15)
axp.set_ylabel('Position (arb. unit)', fontsize = 15)
axp.set_title(f'Position over time with alpha = {alpha}, dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/position_plot_{str}.png')

fige , axe = plt.subplots(1,1)
energy_labels = [f"E({i})" for i in range(5)]

axe.plot(t, energy_array[:,0:5], label = energy_labels)
axe.set_xlabel('Time (arb. unit)', fontsize = 15)
axe.set_ylabel('Energy (arb. unit)', fontsize = 15)
axe.set_title(f'Energy over time with alpha = {alpha}, dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/energy_plot_{str}.png')

figea , axea = plt.subplots(1,1)
#energy_labels = [f"E_average({i})" for i in range()]

axea.plot(t[1:], energy_average_array[1:,:])
axea.set_xlabel('Time (arb. unit)', fontsize = 15)
axea.set_ylabel('Energy (arb. unit)', fontsize = 15)
axea.set_title(f'Time averaged energy with alpha = {alpha}, dt={dt}', fontsize = 15)
plt.yscale('log')
plt.xscale('log')
#plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/energy_average_plot_{str}.png')
