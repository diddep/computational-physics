import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str
from scipy.signal import savgol_filter

sns.set_theme()
set_plot_style.main()

results = unpack_csv.main()
(R1, R2, E_local, E_local_derivative, x_distribution, phi_k, steps_linspace, alpha_results, params) = results

task_str = get_task_str.main()

alpha = alpha_results.alpha.values[0]
N_steps = params.N_steps.values[0]
d_displacement = params.d_displacement.values[0]

R1_norm =  np.square(R1.values)
R1_norm = np.sqrt(np.sum(R1_norm, axis=1))
R2_norm =  np.square(R2.values)
R2_norm = np.sqrt(np.sum(R2_norm, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho

# Calculate moving average of energy

window_size = 10000
moving_averages = savgol_filter(E_local_derivative, window_size, 3)

# Plot energy
fig_energy_deriv, ax_energy_deriv = plt.subplots(1,1)
ax_energy_deriv.plot(steps_linspace, E_local_derivative, alpha = 0.5, label = "Measured energy derivative")
ax_energy_deriv.plot(steps_linspace, moving_averages, 'k--', label = f"Savitzky-Golay with w={window_size}")
ax_energy_deriv.set_xlabel("Steps [a.u.]")
ax_energy_deriv.set_ylabel("Energy [a.u.]")
ax_energy_deriv.set_title(f'Derivative of local energy, alpha = {alpha}')
#ax_energy_deriv.set_title(f'Derivative of local energy, alpha = {alpha}')
ax_energy_deriv.legend()
fig_energy_deriv.savefig(f'plots_python/{task_str}/energy_derivative.png')

