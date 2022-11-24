import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('../eq_of_motion.csv', delimiter=',')#.reshape(-1,7)

t = array[:,0]
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

dt = t[0]
print(dt)

fig2 , axT = plt.subplots(1,1)

axT.plot(t, temp, label='Temperature')


plt.legend(fontsize=10)
axT.set_xlabel('Time (ps)', fontsize = 15)
axT.set_ylabel('Temperature [K]', fontsize = 15)
axT.set_title(f'Temperature over time with timestep dt={dt}, tau=50', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Temp_plots/plot_T_dt{dt}tau50.png')