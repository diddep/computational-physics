import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str

def main(results):
    sns.set_theme()
    set_plot_style.main()

    #results = unpack_csv.main()
    (R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params) = results

    task_str = get_task_str.main()

    alpha = alpha_results.alpha.values[0]
    E_average = alpha_results.E_average.values[0]
    N_steps = params.N_steps.values[0]
    d_displacement = params.d_displacement.values[0]

    # Calculate moving average of energy

    window_size = 10000; poly_order = 3
    moving_averages = savgol_filter(E_local, window_size, poly_order)

    # Plot energy
    fig_energy, ax_energy = plt.subplots(1,1)
    ax_energy.plot(steps_linspace, E_local, alpha = 0.5, label = "Measured energy")
    ax_energy.plot(steps_linspace, moving_averages, 'k--', label = f"Savitzky-Golay with w={window_size}, p={poly_order}")
    ax_energy.hlines(y=E_average, xmin=0, xmax=N_steps, linewidth=2, color='r', linestyles='--', label = f"E_average = {E_average}")
    ax_energy.set_xlabel("Steps [a.u.]")
    ax_energy.set_ylabel("Energy [a.u.]")
    ax_energy.set_title(f'Local energy, alpha = {alpha}')
    #ax_energy.set_title(f'Derivative of local energy, alpha = {alpha}')
    ax_energy.legend(loc="lower left")
    fig_energy.savefig(f'plots_python/{task_str}/energy.png')

if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)
