import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme()


E_of_alpha = pd.read_csv("../csv/E_of_alpha.csv", engine="pyarrow", names= ["E_of_alpha"])

Variance_vec = pd.read_csv("../csv/variance_alpha.csv", engine="pyarrow", names= ["variance"])

alpha_max = 0.25; alpha_min=0.1



E_of_alpha= E_of_alpha.values.astype(float)
Variance_vec = Variance_vec.values.astype(float)

alpha_step = (alpha_max-alpha_min)/len(E_of_alpha[0])

alpha_vec = np.arange(1,len(E_of_alpha[0])+1)* alpha_step + alpha_min

fig_task3, ax_task3 = plt.subplots(1,1)

ax_task3.plot(alpha_vec, E_of_alpha[0])
ax_task3.errorbar(alpha_vec, E_of_alpha[0], yerr = np.sqrt(11*11*Variance_vec[0]/1e6), color = 'r')

fig_task3.savefig("plots_python/task3/E_of_alpha.png")

print((Variance_vec[0]))
print(alpha_vec.shape)
