import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str


sns.set_theme()
set_plot_style.main()

results = unpack_csv.main()
(R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params) = results

task_str = get_task_str.main()

lag_vec = np.arange(0,len(phi_k))

fig_phi_k, ax_phi = plt.subplots(1,1)

ax_phi.plot(lag_vec, phi_k, label = "phi_k")
ax_phi.set_xlabel("lag_vec [a.u.]")
ax_phi.set_ylabel("Phi [a.u.]")
ax_phi.set_title(f'Phi_k')
ax_phi.legend()

fig_phi_k.savefig(f"plots_python/{task_str}/phi_k.png")