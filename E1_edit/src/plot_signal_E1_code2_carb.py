#!/usr/bin/env python
###############################################################################
# E1code4
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt=1e-5; kappa = 1000 /16.0218; m=12.01/9649

pos_array = np.genfromtxt('../task1_5/Saved_trajectory_carb_small_dt.csv', delimiter=',', skip_header=0)
vel_array = np.genfromtxt('../task1_5/Saved_velocity_carb_small_dt.csv', delimiter=',', skip_header=0)

kinetic_energy = m*np.square(vel_array)/2

pos1 = pos_array[:,0]; pos2 = pos_array[:,1]; pos3 = pos_array[:,2];
pot_e1 = kappa * np.square(pos1)/2
pot_e2 = kappa*(np.square(pos1-pos2)+np.square(pos3-pos2))/2
pot_e3 = kappa * np.square(pos3)/2
potential_energy = np.array([pot_e1, pot_e2, pot_e3]).reshape(-1,3)

energy_cart = kinetic_energy[:-1,:] + potential_energy
total_energy = np.sum( energy_cart, axis=1)


fig, ax = plt.subplots(figsize=[12,6])
t = np.arange(0,len(pos_array)*dt, dt)
X_Lim = np.max(t)

ax.plot(t, pos_array)

ax.set_title(f"Trajectory with dt={dt}, carbon")
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
ax.set_xlim(0,X_Lim)
ax.legend(["Particle 1","Particle 2","Particle 3"])
plt.show()
fig.savefig('Saved_trajectory_carb.pdf')

figenergy, axenergy = plt.subplots(figsize=[12,6])

axenergy.plot(t, kinetic_energy[:-1,0]+ kinetic_energy[:-1,1]+ kinetic_energy[:-1,2], color='r', label='total kinetic energy')
axenergy.plot(t, pot_e1+ pot_e2+ pot_e3, color='b', label='total potential energy')
axenergy.plot(t, pot_e1+ pot_e2+ pot_e3+kinetic_energy[:-1,0]+ kinetic_energy[:-1,1]+ kinetic_energy[:-1,2], color='k', label=' total energy' )
axenergy.legend()
axenergy.set_title(f" Time dependence of energy with dt={dt}, carbon")
axenergy.set_xlim(0,X_Lim)

plt.show()
figenergy.savefig('energy_cons_carb.pdf')