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

    def rho(rvec, z):
        rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

        return rho

    ## Plot histogram

    n_bins = 70
    fig_dist, ax_dist = plt.subplots(1,2, figsize=(10,5))

    for idx, R_chain in enumerate([R1, R2]):
        R_norm =  np.square(R_chain)
        R_norm = np.sqrt(np.sum(R_norm, axis=1))
        
        rvec = np.linspace(0.1, np.max(R_norm))
        counts_r, bins_r = np.histogram(R_norm, bins = n_bins, density = True)
        
        ax_dist[idx].stairs(counts_r, bins_r, fill=True)
        ax_dist[idx].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--',label= r'$\rho $ optimized', linewidth=3)
        ax_dist[idx].plot(rvec, rho(rvec, 2), color='k', linestyle=':',label= r'$\rho $ unscreened', linewidth=3)
        
        ax_dist[idx].set_title(f"Distribution for R$_{idx}$, alpha = {alpha}")
        ax_dist[idx].set_xlabel(f"Radius [$a_0$]")
        ax_dist[idx].set_ylabel("Probability density")
        ax_dist[idx].legend()
    plt.tight_layout()
    fig_dist.savefig(f'plots_python/{task_str}/histogram_alpha.png')
    
if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)
