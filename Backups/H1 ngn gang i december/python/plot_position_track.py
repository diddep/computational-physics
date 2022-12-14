import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('../position_track_eq.csv', delimiter=',')#.reshape(-1,7)

t = array[:,0]
q1x = array[:,1]
q1y = array[:,2]
q1z = array[:,3]
q2x = array[:,4]
q2y = array[:,5]
q2z = array[:,6]
q3x = array[:,7]
q3y = array[:,8]
q3z = array[:,9]

dt = t[0]
print(dt)
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
