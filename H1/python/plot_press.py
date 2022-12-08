import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

#str = "eq"
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
t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))

fig3 , axP = plt.subplots(1,1)
axP.plot(t, press, label='Pressure')
axP.set_xlabel('Time (ps)', fontsize = 15)
axP.set_ylabel('Pressure [Bar]', fontsize = 15)
axP.set_title(f'Pressure over time with timestep dt={dt}, tau50', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/pressure.png')
#plt.show()

