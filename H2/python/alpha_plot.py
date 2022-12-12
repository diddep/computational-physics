import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import set_plot_style
import unpack_csv
import get_task_str

sns.set_theme()
set_plot_style.main()

results = unpack_csv.main()
(R1, R2, E_local, E_local_derivative, x_distribution, phi_k, steps_linspace, alpha_results, params) = results


task_str = get_task_str.main()
#TODO: implement more datapoints to plot sequence
alpha_steps = alpha_results.ix.values[0]
average_energy = alpha_results.E_average.values[0]
alpha_task4 = alpha_results.alpha.values[0]


fig_alpha, ax_alpha = plt.subplots(figsize=(10,5))

color = 'tab:blue'
ax_alpha.set_xlabel("Steps",)
ax_alpha.set_ylabel("alpha", color=color)
ax_alpha.plot(alpha_steps, alpha_task4, label = "alpha", color = color)
ax_alpha.tick_params(axis='y', labelcolor=color)

ax_energy = ax_alpha.twinx()

color = 'tab:red'
ax_energy.set_ylabel("Energy [a.u.]")
ax_energy.plot(alpha_steps, average_energy, label = "Average energy", color = color)
ax_energy.tick_params(axis='y', labelcolor=color)

ax_alpha.set_title(f'Evolution of parameter alpha and resulting average energy')

fig_alpha.tight_layout()
fig_alpha.savefig(f'plots_python/{task_str}/alpha_plot_beta.png')