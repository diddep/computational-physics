import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str

def main():
    sns.set_theme()
    set_plot_style.main()

    results = unpack_csv.main()
    (R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params) = results

    task_str = get_task_str.main()

    arr_sin = np.sin(theta_distribution)
    arr_x = np.cos(theta_distribution)

    n_bins =70
    fig_xdist, ax_dist = plt.subplots(1,3, figsize =(15,5))

    arr_str = ["x", r"$\theta$ [rad]", r"$\sin(\theta)$"]
    ax_dist[0].axhline(1/2, label='Distribution for uncorrelated electrons', color='r', linewidth=2, linestyle='dashed')

    for idx, arr_dist in enumerate([arr_x, theta_distribution, arr_sin]):
        counts, bins = np.histogram(arr_dist, bins = n_bins, density = True)
        ax_dist[idx].stairs(counts, bins, fill = True, label='Sampled distribution')
        ax_dist[idx].set_title('Distribution for ' + arr_str[idx])
        ax_dist[idx].set_xlabel(f'{arr_str[idx]}')
        ax_dist[idx].set_ylabel('Probability density')
        #ax_dist[idx].legend()

    plt.tight_layout()

    fig_xdist.savefig(f'plots_python/{task_str}/xdist_new.png')

main()