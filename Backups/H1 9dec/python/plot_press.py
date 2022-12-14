import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

str = "eq"
#str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')
array_prod = np.genfromtxt(f'../csv/vel_verlet_prod.csv', delimiter=',', skip_header=1)
parameters_prod = np.genfromtxt(f'../csv/parameters_prod.csv', delimiter=',')

array = array_prod

#t = array[:,0]
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

end_time = parameters[-1,0]
dt = parameters[-1,1]
lattice_param = parameters[-1,2]
temp_scaling = parameters[-1,3]
press_scaling = parameters[-1,4]
temp_eq = parameters[-1,5]
press_eq = parameters[-1,6]
tau_T = parameters[-1,7]
tau_P = parameters[-1,8]

t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))

fig3 , axP = plt.subplots(1,1)
axP.plot(t, press, label='Pressure')

average_pressure = np.mean(press)

# add dashed line for the average energy
axP.axhline(
    average_pressure,
    linestyle='--',
    color='k',
    linewidth=2
)
axP.annotate(f'P$_{{average}}$ = {average_pressure:.2f}', xy=(t[-1], average_pressure), xytext=(-150, 150),
            textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

axP.set_xlabel('Time (ps)', fontsize = 15)
axP.set_ylabel('Pressure [Bar]', fontsize = 15)
axP.set_title(f'Pressure, dt={dt}', fontsize = 15)
# if(press_scaling):
#     axP.set_title(f'Pressure, dt={dt:.2f}, tau_P={tau_P:.2f}', fontsize = 15)
# else:
#     axP.set_title(f'Temperature, dt={dt:.2f}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/pressure.png')
#plt.show()

