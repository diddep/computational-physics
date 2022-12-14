import numpy as np
import matplotlib.pyplot as plt

str = "eq"
#str = "prod"

array = np.genfromtxt(f'../vel_verlet_{str}.csv', delimiter=',')
parameters = np.genfromtxt(f'../parameters_{str}.csv', delimiter=',')
    
end_time = parameters[0]
dt = parameters[1]
lattice_param = parameters[2]
temp_scaling = parameters[3]
press_scaling = parameters[4]
temp_eq = parameters[5]
press_eq = parameters[6]
tau_T = parameters[7]
tau_P = parameters[8]

t = array[:,0]
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

fig , ax = plt.subplots(1,1)

ax.plot(t, e_pot, label='potential')
ax.plot(t, e_kin, label='kinetic')
ax.plot(t, e_tot, label='total')

ax.set_xlabel('Time (ps)', fontsize = 15)
ax.set_ylabel('Energy (eV/unit cell)', fontsize = 15)
if(str == "eq"):
    ax.set_title(f'Energy over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
else:
    ax.set_title(f'Energy over time with: dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_energy_dt{dt}_{str}.png')

fig2 , axT = plt.subplots(1,1)

axT.plot(t, temp, label='Temperature')


plt.legend(fontsize=10)
axT.set_xlabel('Time (ps)', fontsize = 15)
axT.set_ylabel('Temperature [K]', fontsize = 15)
if(str == "eq"):
    axT.set_title(f'Temperature over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
else:
    axT.set_title(f'Temperature over time with: dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_T_dt{dt}_{str}.png')

fig3 , axP = plt.subplots(1,1)
axP.plot(t, press, label='Pressure')
axP.set_xlabel('Time (ps)', fontsize = 15)
axP.set_ylabel('Pressure [Bar]', fontsize = 15)
if(str == "eq"):
    axP.set_title(f'Pressure over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
else:
    axP.set_title(f'Pressure over time with: dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Energy_plots/plot_P_dt{dt}_{str}.png')
#plt.show()

