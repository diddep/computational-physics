#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
array = np.genfromtxt('../signal.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()
ax.plot(array[:, 0], array[:, 1])
l_sig = len(array[:, 0])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.set_title('dt=0.1/6')
ax.grid()

fig.savefig('signalsuperposedt16.pdf')


array_pow = np.genfromtxt('../powerspectrum_shift.csv', delimiter=',', skip_header=1)
print(array_pow)

fig_pow, ax_pow = plt.subplots()


ax_pow.plot(array_pow[:,1], array_pow[:,0])
ax_pow.set_xlabel('frequency (arb.unit)')
ax_pow.set_ylabel('intensity (arb.unit)')
ax_pow.set_title('dt=0.1/6')
fig_pow.savefig('pow_spectrum_shift_superposdt16.pdf')


plt.show()

