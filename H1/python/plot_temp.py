import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# str = "eq"
str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

#t = array[:,0]
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

dt = parameters[-1,1]

#dt = t[0]
#dt = array[0,0]
t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))

fig2 , axT = plt.subplots(1,1)

axT.plot(t, temp, label='Temperature')


plt.legend(fontsize=10)
axT.set_xlabel('Time (ps)', fontsize = 15)
axT.set_ylabel('Temperature [K]', fontsize = 15)
axT.set_title(f'Temperature over time with timestep dt={dt}, tau=50', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/temperature.png')