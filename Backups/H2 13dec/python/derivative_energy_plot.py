import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
from scipy import signal

def main(results):
    sns.set_theme()
    set_plot_style.main()

    #results = unpack_csv.main()
    (R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params) = results
    
    alpha = alpha_results.alpha.values[0]
    E_average = alpha_results.E_average.values[0]
    N_steps = params.N_steps.values[0]
    d_displacement = params.d_displacement.values[0]

    window_size = 10000; poly_order = 3
    moving_averages = signal.savgol_filter(E_local, window_size, poly_order)
    # Plot energy
    fig_energy_deriv, ax_energy_deriv = plt.subplots(1,1)
    ax_energy_deriv.plot(steps_linspace, E_local_derivative, alpha = 0.5, label = "Measured energy derivative")
    ax_energy_deriv.plot(steps_linspace, moving_averages, 'k--', label = f"Savitzky-Golay with w={window_size}, p={poly_order}")
    ax_energy_deriv.set_xlabel("Steps [a.u.]")
    ax_energy_deriv.set_ylabel("Energy [a.u.]")
    ax_energy_deriv.set_title(f'Derivative of local energy, alpha = {alpha}')
    #ax_energy_deriv.set_title(f'Derivative of local energy, alpha = {alpha}')
    fig_energy_deriv.savefig(f'plots_python/{task_str}/energy_derivative.png')

if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)

