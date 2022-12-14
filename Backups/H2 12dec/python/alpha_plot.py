import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


R1 =np.genfromtxt("../R1.csv", delimiter=',')
R2 =np.genfromtxt("../R2.csv", delimiter=',')

Energy_local = np.genfromtxt("../E_local.csv", delimiter=',')
Energy_local_derivative = np.genfromtxt("../E_local_derivative.csv", delimiter=',')
xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
theta_csv = np.genfromtxt("../theta.csv", delimiter=',')
theta = theta_csv[:,1]
Phi_k =np.genfromtxt("../phi_k.csv", delimiter=',')

alpha_results = np.genfromtxt("../alpha_results.csv", delimiter=',')
alpha_steps = alpha_results[:,0]
average_energy = alpha_results[:,1]
alpha_task4 = alpha_results[:, 2]


params = np.genfromtxt("../alpha_params.csv", delimiter=',')

n_alpha_steps = int(params[0])
n_discarded_steps = int(params[1])
alpha = params[2]
A = params[3]
beta = params[4]
n_steps = int(params[5])
d_displacement = params[6]
is_task1 = params[7]
is_task2 = params[8]
is_task3 = params[9]
is_task4 = params[10]

if is_task1:
    task_str = "task1"
elif is_task2:
    task_str = "task2"
elif is_task3:
    task_str = "task3"
elif is_task4:
    task_str = "task4"

print(alpha_steps.shape)
print(alpha_task4.shape)
print(average_energy.shape)

fig_alpha, ax_alpha = plt.subplots(figsize=(10,5))

color = 'tab:blue'
ax_alpha.set_xlabel("Steps", fontsize=15)
ax_alpha.set_ylabel("alpha", fontsize=15, color=color)
ax_alpha.plot(alpha_steps, alpha_task4, label = "alpha", color = color)
ax_alpha.tick_params(axis='y', labelcolor=color)

ax_energy = ax_alpha.twinx()

color = 'tab:red'
ax_energy.set_ylabel("Energy [a.u.]", fontsize=15)
ax_energy.plot(alpha_steps, average_energy, label = "Average energy", color = color)
ax_energy.tick_params(axis='y', labelcolor=color)

ax_alpha.set_title(f'Evolution of parameter alpha and resulting average energy', fontsize=15)

#ax_alpha.legend()
fig_alpha.tight_layout()
fig_alpha.savefig(f'plots_python/{task_str}/alpha_plot_beta{beta}.png')