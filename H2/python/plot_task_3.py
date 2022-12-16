# import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# sns.set_theme()


# E_of_alpha = pd.read_csv("../csv/E_of_alpha.csv", engine="pyarrow", names= ["E_of_alpha"])


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import set_plot_style
sns.set_theme()
set_plot_style.main()
from scipy.optimize import curve_fit

E_of_alpha = pd.read_csv("../csv/E_of_alpha.csv", engine="pyarrow", names= ["E_of_alpha"])

Variance_vec = pd.read_csv("../csv/variance_alpha.csv", engine="pyarrow", names= ["variance"])

alpha_max = 0.25; alpha_min=0.05

def polynomial(x, a, b, c):
    return a * x**2 + b * x + c


E_of_alpha= E_of_alpha.values.astype(float)
Variance_vec = Variance_vec.values.astype(float)

alpha_step = (alpha_max-alpha_min)/len(E_of_alpha[0])

alpha_vec = np.arange(1,len(E_of_alpha[0])+1)* alpha_step + alpha_min

fig_task3, ax_task3 = plt.subplots(1,1)

upper_bound = E_of_alpha[0]+np.sqrt(Variance_vec[0]/(1e7))
lower_bound = E_of_alpha[0]-np.sqrt(Variance_vec[0]/(1e7))
fs = 15

ax_task3.plot(alpha_vec, E_of_alpha[0], label = "Average energy", color='r', linewidth=2)
#ax_task3.errorbar(alpha_vec, E_of_alpha[0], yerr = np.sqrt(Variance_vec[0]/1e6), color = 'r')
ax_task3.fill_between(alpha_vec,lower_bound,upper_bound, label=r"$\pm\sigma$", alpha=0.5)

params, covariance = curve_fit(polynomial, alpha_vec, E_of_alpha[0])

alpha_vals = np.linspace(alpha_min, alpha_max, 50)


ax_task3.legend(fontsize = fs)
ax_task3.set_title(r"Averag local energy as function of parameter $\alpha$", fontsize=fs)
ax_task3.set_xlabel(r"$\alpha$", fontsize = fs+4)
ax_task3.set_ylabel(r"Energy [$E_h$]", fontsize = fs)
opt_e_ind = np.argmin(E_of_alpha[0])
opt_alpha = alpha_vec[opt_e_ind]
min_E =np.min(E_of_alpha[0])

#ax_task3.plot(alpha_vals, polynomial(alpha_vals, *params), linestyle=':',linewidth=4, label=f'Curve fit, \n y = {params[0]:.2} x^2 + {params[1]:.2} x + {params[2]:.2}')


st =r"Optimal $\alpha =$" +f"{opt_alpha:.4f}" + r" with $\langle E_L\rangle$= "+ f"{min_E:.4f}"

ax_task3.scatter(opt_alpha, min_E, color='k',label = st, marker ="+", s =300)

plt.tight_layout()


ax_task3.legend(fontsize = fs-2)
fig_task3.savefig("plots_python/task3/E_of_alpha.png")

print((Variance_vec[0]))
print(alpha_vec.shape)

opt_e_ind = np.argmin(E_of_alpha[0])

opt_alpha = alpha_vec[opt_e_ind]

print("optimal alpha= ", opt_alpha)
print("with energy", np.min(E_of_alpha[0]))


fig_hist, ax_hist = plt.subplots(1,1)

