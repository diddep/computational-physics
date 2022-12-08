import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

str = "eq"
#str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')
array_prod = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters_prod = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

q1x = array[:,1]
q1y = array[:,2]
q1z = array[:,3]
q2x = array[:,4]
q2y = array[:,5]
q2z = array[:,6]
q3x = array[:,7]
q3y = array[:,8]
q3z = array[:,9]

end_time = parameters[-1,0]
dt = parameters[-1,1]

#dt = array[0,0]
t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))
str = "calibration_melt"

fig , xax = plt.subplots(1,1)

xax.plot(t, q1x, label='q1')
xax.plot(t, q2x, label='q2')
xax.plot(t, q3x, label='q3')

xax.set_xlabel('Time (ps)', fontsize = 15)
xax.set_ylabel('X-coordinate (Å)', fontsize = 15)
xax.set_title(f'X-coordinate over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Position_plots/plot_X{dt}_{str}.png')

fig , yax = plt.subplots(1,1)

yax.plot(t, q1y, label='q1')
yax.plot(t, q2y, label='q2')
yax.plot(t, q3y, label='q3')

yax.set_xlabel('Time (ps)', fontsize = 15)
yax.set_ylabel('Y-coordinate (Å)', fontsize = 15)
yax.set_title(f'Y-coordinate over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Position_plots/plot_Y{dt}_{str}.png')

fig , zax = plt.subplots(1,1)

zax.plot(t, q1z, label='q1')
zax.plot(t, q2z, label='q2')
zax.plot(t, q3z, label='q3')

zax.set_xlabel('Time (ps)', fontsize = 15)
zax.set_ylabel('Z-coordinate (Å)', fontsize = 15)
zax.set_title(f'Z-coordinate over time with timestep dt={dt}', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'../Position_plots/plot_Z{dt}_{str}.png')
