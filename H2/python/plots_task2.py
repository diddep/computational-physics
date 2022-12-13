import numpy as np
import pandas as pd
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

    lag_vec = np.arange(0,len(phi_k))

    block_average = pd.read_csv("../csv/block_avg_vec.csv", engine="pyarrow", names= ["block_average"])
    block_average = block_average.values[:,0].astype(float)

    t_relaxation = np.argmin(np.abs(phi_k-np.exp(-2)))

    #stat_ineff_cor = 2*np.sum(phi_k[:t_relaxation]) #factor 2 from fact that it is symmetric, -1 because we count k=0 twice
    stat_ineff_cor = 2*np.sum(phi_k)-1 #factor 2 from fact that it is symmetric, -1 because we count k=0 twice
    
    stat_ineff_block_av = np.average(block_average[150:])
    print(f"------\nstatistical inefficiency calculated from correlation funcion = {stat_ineff_cor}\n----")
    print(f"------\nstatistical inefficiency calculated from block averaging = {stat_ineff_block_av}\n----")
    print("relaxation time =",t_relaxation)


    fig_phi_k, ax_phi = plt.subplots(1,1)

    ax_phi.plot(lag_vec, phi_k, label = "phi_k")
    
    # dom här är enhetslösa
    ax_phi.set_xlabel(r"$k$")
    ax_phi.set_ylabel(r"$\Phi_k$")
    ax_phi.set_title(r'correlation function $\Phi_k$')
    ax_phi.legend()

    fig_phi_k.savefig(f"plots_python/{task_str}/phi_k.png")

    fig_block_avg, ax_blav = plt.subplots(1,1)
    

    block_size_vec = np.arange(1, len(block_average)+1)

    #ax_blav.plot(block_size_vec, block_average)
    ax_blav.set_xlabel("Block size")
    ax_blav.set_ylabel(r"$n_s$")
    ax_blav.set_title(r"Statistical inefficiency $n_s$ calculated from block averaging")
    ax_blav.axhline(stat_ineff_cor, label =r"$n_s$ calculated from correlation function", color="r", linewidth=4)
    ax_blav.legend()
    ax_blav.scatter(block_size_vec, block_average, facecolor ="none", edgecolor="k", alpha=0.8)

    fig_block_avg.savefig(f"plots_python/task2/block_avg.png") 
    

    


if(__name__ == "__main__"):
    results = unpack_csv.main()
    main(results)