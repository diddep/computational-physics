import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv

sns.set_theme()
set_plot_style.main()

results, params = unpack_csv.main()
(R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, steps_linspace) = results
(n_alpha_steps, n_discarded_steps, alpha, A, beta, n_steps, d_displacement, is_task1, is_task2, is_task3, is_task4) = params

# R1 =np.genfromtxt("../R1.csv", delimiter=',')
# R2 =np.genfromtxt("../R2.csv", delimiter=',')

# Energy_local = np.genfromtxt("../E_local.csv", delimiter=',')
# Energy_local_derivative = np.genfromtxt("../E_local_derivative.csv", delimiter=',')
# xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
# theta_csv = np.genfromtxt("../theta.csv", delimiter=',')
# theta = theta_csv[:,1]
# Phi_k =np.genfromtxt("../phi_k.csv", delimiter=',')

# alpha_results = np.genfromtxt("../alpha_results.csv", delimiter=',')


# params = np.genfromtxt("../alpha_params.csv", delimiter=',')

# n_alpha_steps = int(params[0])
# n_discarded_steps = int(params[1])
# alpha = params[2]
# A = params[3]
# beta = params[4]
# n_steps = int(params[5])
# d_displacement = params[6]
# is_task1 = params[7]
# is_task2 = params[8]
# is_task3 = params[9]
# is_task4 = params[10]

if is_task1:
    task_str = "task1"
elif is_task2:
    task_str = "task2"
elif is_task3:
    task_str = "task3"
elif is_task4:
    task_str = "task4"

R1_norm =  np.square(R1)
R1_norm = np.sqrt(np.sum(R1_norm, axis=1))
R2_norm =  np.square(R2)
R2_norm = np.sqrt(np.sum(R2_norm, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho

## Calculate moving average of energy
print(E_local_derivative.to_numpy()[:,0])

window_size = 1000; i = 0; moving_averages = []
moving_average_length = len(E_local_derivative) - window_size + 1
moving_average_linspace = np.linspace(0, moving_average_length, moving_average_length)
while i < moving_average_length:
    window = E_local_derivative.to_numpy()[i : i + window_size,0]
    #window = E_local_derivative[i : i + window_size]
    window_average = round(sum(window) / window_size, 4)
    moving_averages.append(window_average)
    i += 1

## Plot energy

fig_energy_deriv, ax_energy_deriv = plt.subplots(1,1)
ax_energy_deriv.plot(steps_linspace, E_local_derivative, alpha = 0.5, label = "Measured energy derivative")
ax_energy_deriv.plot(moving_average_linspace, moving_averages, 'k--', label = "Moving average of energy derivative")
ax_energy_deriv.set_xlabel("Steps [a.u.]")
ax_energy_deriv.set_ylabel("Energy [a.u.]")
ax_energy_deriv.set_title(f'Derivative of local energy, averaged with window of {window_size}, alpha = {alpha}')
ax_energy_deriv.legend()
fig_energy_deriv.savefig(f'plots_python/{task_str}/energy_derivative_alpha{alpha}_nsteps{n_steps}_d{d_displacement}.png')

