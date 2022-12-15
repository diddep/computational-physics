import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str

def main(results):
    sns.set_theme()
    set_plot_style.main()

    # results = unpack_csv.main()
    (R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params) = results

    task_str = get_task_str.main()
    alpha_steps = alpha_results.ix.values[:]
    average_energy = alpha_results.E_average.values[:]
    alpha_task4 = alpha_results.alpha.values[:]
    beta_task4 = alpha_results.beta.values[:]
    N_alpha_steps = int(params.N_alpha_steps.values[0])

    values_unique_beta, counts_per_unique_beta = np.unique(beta_task4, return_counts=True)
    
    n = N_alpha_steps
    sliced_alphas = [alpha_task4[i * n:(i + 1) * n] for i in range((len(alpha_task4) + n - 1) // n )]
    sliced_energy = [average_energy[i * n:(i + 1) * n] for i in range((len(average_energy) + n - 1) // n )]
    
    slice_linspace = np.linspace(0, n, n, endpoint=False)

    change = np.empty(shape=(len(values_unique_beta), N_alpha_steps-1))
    for idx in range(N_alpha_steps-1):
        for jdx in range(len(values_unique_beta)):
            change[jdx][idx] = sliced_alphas[jdx][idx+1] - sliced_alphas[jdx][idx]
    

    fig_beta, ax_beta = plt.subplots(1,1)
    fig_change, ax_change = plt.subplots(1,1)
    for idx in range(len(values_unique_beta)):
        #print(sliced_alphas[idx][-1])
        #print("alpha after 200 iterations", sliced_alphas[idx][-1])
        ax_beta.plot(slice_linspace, sliced_alphas[idx][:], label = rf"$\beta$={values_unique_beta[idx]}")
        #print("Average change for last 100 iterations", np.mean(change[idx][-100:]))
        #print(f"{np.mean(change[idx][-100:]):.1e}")
        print(sliced_energy[idx][-1])
        ax_change.plot(slice_linspace[:-1], change[idx][:])

    ax_beta.set_xlabel("Steps",)
    ax_beta.set_ylabel(r"$\alpha$")
    ax_beta.set_title(r'Evolution of $\alpha$ for different values of $\beta$')
    ax_beta.set_xlim(0,200)
    #ax_beta.set_ylim(0.12,0.14)
    ax_beta.legend()
    #ax_beta.set_yscale("log")
    fig_beta.tight_layout()
    fig_beta.savefig(f'plots_python/{task_str}/beta_plot.png')

    ax_change.set_xlabel("Steps",)
    ax_change.set_ylabel("Change in alpha")
    ax_change.set_title(f'Change of  alpha for different values of beta')
    ax_change.set_xlim(150,200)
    ax_change.set_ylim(-0.0005,0.0005)
    #ax_change.set_yscale("log")
    fig_change.tight_layout()
    fig_change.savefig(f'plots_python/{task_str}/change_plot.png')

if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)