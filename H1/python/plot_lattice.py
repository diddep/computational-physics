import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# set default figure size
plt.rcParams["figure.figsize"] = [8, 6]

str = "eq"
#str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')
array_prod = np.genfromtxt(f'../csv/vel_verlet_prod.csv', delimiter=',', skip_header=1)
parameters_prod = np.genfromtxt(f'../csv/parameters_prod.csv', delimiter=',')

array= array_prod

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
#dt = array[0,0]
t = dt * np.linspace(0,len(array[:,0]), len(array[:,0]))

cell_length = array[:,1]
lattice_length = cell_length/4


# create figure and axes for energy plot
figlattice , ax_lattice = plt.subplots(1,1)

# plot energy data
ax_lattice.plot(t, lattice_length, label='Lattice parameter')

average_lattice = np.mean(lattice_length)

# add dashed line for the average energy
ax_lattice.annotate(f'a$_{{0}}$ = {lattice_length[-1]:.2f}', xy=(t[-1], lattice_length[-1]), xytext=(-50, -50),
            textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

# set labels and title
ax_lattice.set_xlabel('Time (ps)', fontsize = 15)
ax_lattice.set_ylabel('Lattice parameter (Å)', fontsize = 15)
ax_lattice.set_title(f'Lattice parameter, dt={dt}', fontsize = 15)
# if(str == "eq"):
#     ax_lattice.set_title(f'Lattice parameter, dt={dt}, $tau_T$={tau_T}, $tau_P$={tau_P}', fontsize = 15)
# else:
#     ax_lattice.set_title(f'Lattice parameter, dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/lattice.png')