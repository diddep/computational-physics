import numpy as np
import matplotlib.pyplot as plt

# set default figure size
plt.rcParams["figure.figsize"] = [8, 6]

#str = "eq"
str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

end_time = parameters[-1,0]
dt = parameters[-1,1]
lattice_param = parameters[-1,2]
temp_scaling = parameters[-1,3]
press_scaling = parameters[-1,4]
temp_eq = parameters[-1,5]
press_eq = parameters[-1,6]
tau_T = parameters[-1,7]
tau_P = parameters[-1,8]

# extract columns from array
dt = array[0,0]
t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

# create figure and axes for energy plot
fig , ax = plt.subplots(1,1)

# plot energy data
ax.plot(t, e_pot, label='potential')
ax.plot(t, e_kin, label='kinetic')
ax.plot(t, e_tot, label='total')

# set labels and title
ax.set_xlabel('Time (ps)', fontsize = 15)
ax.set_ylabel('Energy (eV/unit cell)', fontsize = 15)
if(str == "eq"):
    ax.set_title(f'Energy over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
else:
    ax.set_title(f'Energy over time with: dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/energy.png')