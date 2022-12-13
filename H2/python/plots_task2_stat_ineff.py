import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


arr_phi_k = np.genfromtxt("../csv/phi_k.csv", delimiter=',')

print(arr_phi_k[:,0])
print(len(arr_phi_k[:,0]))


lag_vec = np.arange(0,len(arr_phi_k[:,0]))

fig_phi_k, ax_phi = plt.subplots(1,1)

ax_phi.plot(lag_vec,arr_phi_k[:,0])

fig_phi_k.savefig("plots_python/phi_k.png")