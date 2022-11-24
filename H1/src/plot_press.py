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

fig3 , axP = plt.subplots(1,1)
axP.plot(t, press, label='Pressure')
axP.set_xlabel('Time (ps)', fontsize = 15)
axP.set_ylabel('Pressure [?]', fontsize = 15)
axP.set_title(f'Pressure over time with timestep dt={dt}, tau50', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Press_plots/plot_P_dt{dt}_tau50.png')
#plt.show()

