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
str = ["TRUE"]
fig , ax = plt.subplots(1,1)

ax.plot(t, e_pot, label='potential')
ax.plot(t, e_kin, label='kinetic')
ax.plot(t, e_tot, label='total')

ax.set_xlabel('Time (ps)', fontsize = 15)
ax.set_ylabel('Energy (eV/unit cell)', fontsize = 15)
ax.set_title(f'Energy over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_energy_dt{dt}_{str}.png')

fig2 , axT = plt.subplots(1,1)

axT.plot(t, temp, label='Temperature')


plt.legend(fontsize=10)
axT.set_xlabel('Time (ps)', fontsize = 15)
axT.set_ylabel('Temperature [K]', fontsize = 15)
axT.set_title(f'Temperature over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_T_dt{dt}_{str}.png')

fig3 , axP = plt.subplots(1,1)
axP.plot(t, press, label='Pressure')
axP.set_xlabel('Time (ps)', fontsize = 15)
axP.set_ylabel('Pressure [GPa]', fontsize = 15)
axP.set_title(f'Pressure over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_P_dt{dt}_{str}.png')
#plt.show()

