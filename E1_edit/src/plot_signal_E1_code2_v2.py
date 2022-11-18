#!/usr/bin/env python
###############################################################################
# E1code4
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt=1e-5; kappa = 1000 /16.0218; m=12.01/9649

pos_array = np.genfromtxt('../Saved_trajectory2.csv', delimiter=',', skip_header=0)
vel_array = np.genfromtxt('../Saved_velocity2.csv', delimiter=',', skip_header=0)

kinetic_energy = m*np.square(vel_array)/2
print(np.average(kinetic_energy, axis=0))

pos1 = pos_array[:,0]; pos2 = pos_array[:,1]; pos3 = pos_array[:,2];
pot_e1 = kappa * np.square(pos1)/2
pot_e2 = kappa*(np.square(pos1-pos2)+np.square(pos3-pos2))/2
pot_e3 = kappa * np.square(pos3)/2

potential_energy = np.array([pot_e1, pot_e2, pot_e3]).reshape(-1,3) #kappa * np.square(np.sum(pos_array))/2
print(potential_energy,'pot')
#print(kinetic_energy[:-1,:] + potential_energy, 'sum')

energy_cart = kinetic_energy[:-1,:] + potential_energy

total_energy = np.sum( energy_cart, axis=1)#np.sum( np.abs(kinetic_energy[:-1,:] + potential_energy), axis=1)
print(total_energy)


fig, ax = plt.subplots(figsize=[12,6])

#t = np.linspace(0,len(array)*dt,len(array))
t = np.arange(0,len(pos_array)*dt, dt)
#print(t.shape)
ax.plot(t, pos_array)
#ax.plot(t, total_energy)
ax.set_title(f"Trajectory with dt={dt}")
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
ax.set_xlim(0,np.max(t))
ax.legend(["Particle 1","Particle 2","Particle 3"])
plt.show()
fig.savefig('Saved_trajectory.pdf')

figenergy, axenergy = plt.subplots(figsize=[12,6])
axenergy.plot(t, total_energy)
#axenergy.plot(t, kinetic_energy[:-1,0]+ kinetic_energy[:-1,1]+ kinetic_energy[:-1,2], color='r')
#axenergy.plot(t, pot_e1, color='k' )
#axenergy.plot(t, pot_e3, color='r' )
axenergy.plot(t,  pot_e2,color='g' )
axenergy.plot(t,  kinetic_energy[:-1,:] )

plt.show()