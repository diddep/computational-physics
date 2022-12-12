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
average_energy = alpha_results[1]


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

steps = np.linspace(0, n_steps, n_steps)

R1_norm =  np.square(R1)
R1_norm = np.sqrt(np.sum(R1_norm, axis=1))
R2_norm =  np.square(R2)
R2_norm = np.sqrt(np.sum(R2_norm, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho


## Calculate moving average of energy

window_size = 500; i = 0; moving_averages = []
moving_average_length = len(Energy_local) - window_size + 1
moving_average_linspace = np.linspace(0, moving_average_length, moving_average_length)
while i < moving_average_length:
    window = Energy_local[i : i + window_size]
    window_average = round(sum(window) / window_size, 4)
    moving_averages.append(window_average)
    i += 1

## Plot energy

fig_energy, ax_energy = plt.subplots(1,1, figsize=(10,5))
ax_energy.plot(steps, Energy_local, alpha = 0.5, label = "Measured energy")
ax_energy.plot(moving_average_linspace, moving_averages, 'b--', label = "Moving average of energy")
ax_energy.hlines(y=average_energy, xmin=0, xmax=n_steps, linewidth=2, color='r', linestyles='--', label = f"E_average = {average_energy}")
#ax_energy.text(0, average_energy * (1-0.05), f"E_average = {average_energy}")
ax_energy.set_xlabel("Steps [a.u.]", fontsize=15)
ax_energy.set_ylabel("Energy [a.u.]", fontsize=15)
ax_energy.set_title(f'Local energy, averaged with window of {window_size}, alpha = {alpha}', fontsize=15)
ax_energy.legend()
fig_energy.savefig(f'plots_python/{task_str}/energy_alpha{alpha}_nsteps{n_steps}.png')