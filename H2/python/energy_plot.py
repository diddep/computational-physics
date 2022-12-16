import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
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
    #print(N_steps)
    #print(E_local)
    #window_size = int(N_steps/10); poly_order = 2
    #moving_averages = savgol_filter(E_local, window_size, poly_order)

    # Plot energy
    fig_energy, ax_energy = plt.subplots(1,1)
    ax_energy.plot(steps_linspace, E_local, alpha = 0.5, label = "Measured energy")
    #ax_energy.plot(steps_linspace, moving_averages, 'k--', label = f"Averaged energy (w={window_size}, p={poly_order})")
    #ax_energy.hlines(y=E_average, xmin=0, xmax=N_steps, linewidth=2, color='r', linestyles='--', label = f"E_average = {E_average}")
    ax_energy.vlines(x=1000, ymin=min(E_local), ymax=max(E_local), linewidth=2, color='r', linestyles='--', label = f"N$_{{eq}}$")
    ax_energy.set_xlabel("Steps [a.u.]")
    ax_energy.set_ylabel("Energy [a.u.]")
    ax_energy.set_title(f'Local energy, alpha = {alpha}')
    ax_energy.legend(loc="lower right")
    fig_energy.savefig(f'plots_python/{task_str}/energy.png')
    print(task_str)


    fig_hist, ax_hist = plt.subplots(1,1)
    n_bins =50
    counts, bins = np.histogram(E_local, bins = n_bins, density = True)
    ax_hist.stairs(counts, bins, fill = True, label='Sampled distribution')
    ax_hist.axvline(E_average, color = 'r', label=rf"$\langle E_L\rangle ={E_average}$")
    ax_hist.legend()
    fig_hist.savefig(f'plots_python/{task_str}/histo_energy.png')
    
    

if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)