#!/usr/bin/env python
###############################################################################
# E1code4
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
dt=1e-3/2
array = np.genfromtxt('../Saved_trajectory.csv', delimiter=',', skip_header=0)
fig, ax = plt.subplots(figsize=[12,6])

#t = np.linspace(0,len(array)*dt,len(array))
t = np.arange(0,len(array)*dt, dt)
#print(t.shape)
ax.plot(t, array)
ax.set_title(f"Trajectory with dt={dt}")
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
ax.set_xlim(0,0.25)
ax.legend(["Particle 1","Particle 2","Particle 3"])

fig.savefig('Saved_trajectory.pdf')