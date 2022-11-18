#!/usr/bin/env python
###############################################################################
# E1code4
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt=1e-5; kappa = 1600 /16.0218; m=12.01/9649; mo=15.994/9649

pos_array = np.genfromtxt('../Task1_6/Saved_trajectory_co2_small_dt.csv', delimiter=',', skip_header=0)
vel_array = np.genfromtxt('../Task1_6/Saved_velocity_co2_small_dt.csv', delimiter=',', skip_header=0)

kinetic_energy = np.square(vel_array)/2
kinetic_energy[:,0] *= mo
kinetic_energy[:,1] *= m
kinetic_energy[:,2] *= mo


pos1 = pos_array[:,0]; pos2 = pos_array[:,1]; pos3 = pos_array[:,2]
pot_e1 = 0 #kappa * np.square(pos1)/2
pot_e2 = kappa*(np.square(pos1-pos2)+np.square(pos3-pos2))/2
pot_e3 = 0 # kappa * np.square(pos3)/2
#potential_energy = np.array([pot_e1, pot_e2, pot_e3]).reshape(-1,3)
potential_energy = np.transpose(np.array([pot_e2])).reshape(-1,1)

#energy_cart = kinetic_energy[:-1,0] + kinetic_energy[:-1,1] + kinetic_energy[:-1,2] + potential_energy
#total_energy = np.sum(energy_cart, axis=1)
total_energy = np.sum(kinetic_energy, axis=1).reshape(-1,1)
print(total_energy[:-1].shape)
print(potential_energy.shape)
total_energy = total_energy[:-1] + potential_energy


fig, ax = plt.subplots(figsize=[12,6])
t = np.arange(0,len(pos_array)*dt, dt)
X_Lim = np.max(t)

ax.plot(t, pos_array)

ax.set_title(f"Trajectory with dt={dt}, carbon dioxide")
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
ax.set_xlim(0,X_Lim)
ax.legend(["Particle 1","Particle 2","Particle 3"])
plt.show()
fig.savefig('Saved_trajectory_co2.pdf')

figenergy, axenergy = plt.subplots(figsize=[12,6])

axenergy.plot(t, np.sum(kinetic_energy[:-1,:], axis=1), color='r', label='total kinetic energy')
axenergy.plot(t, potential_energy, color='b', label='total potential energy')
axenergy.plot(t, total_energy, color='k', label=' total energy' )
axenergy.legend()
axenergy.set_title(f" Time dependence of energy with dt={dt}, carbon dioxide")
axenergy.set_xlim(0,X_Lim)

plt.show()
figenergy.savefig('energy_cons_co2.pdf')