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

#array = array_prod

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


average_temperature = np.mean(array[:,5])
print(average_temperature)

t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))

fig2 , axT = plt.subplots(1,1)

axT.plot(t, temp, label='Temperature')
print(temp_scaling)

# add dashed line for the average energy
axT.axhline(
    average_temperature,
    linestyle='--',
    color='k',
    linewidth=2
)
axT.annotate(f'T$_{{average}}$ = {average_temperature:.2f}', xy=(t[-1], average_temperature), xytext=(-110, 200),
            textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

plt.legend(fontsize=10)
axT.set_xlabel('Time (ps)', fontsize = 15)
axT.set_ylabel('Temperature [K]', fontsize = 15)
axT.set_title(f'Temperature, dt={dt}', fontsize = 15)
# if(temp_scaling):
#     axT.set_title(f'Temperature, dt={dt}, tau_T={tau_T:.2f}', fontsize = 15)
# else:
#     axT.set_title(f'Temperature, dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/temperature.png')