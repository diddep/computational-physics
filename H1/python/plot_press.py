import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('../csv/eq_of_motion.csv', delimiter=',', skip_header=1)#.reshape(-1,7)

#t = array[:,0]
e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

dt = array[0,0]
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

