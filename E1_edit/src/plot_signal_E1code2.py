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
ax.set_title('N=')
ax.grid()

fig.savefig('signalsuperpose1.pdf')


array_pow = np.genfromtxt('../powerspectrum_shift.csv', delimiter=',', skip_header=1)
print(array_pow)

fig_pow, ax_pow = plt.subplots()


ax_pow.plot(array_pow[:,1], array_pow[:,0])
ax_pow.set_xlabel('frequency (arb.unit)')
ax_pow.set_ylabel('intensity (arb.unit)')
ax_pow.set_title('N=750')
fig_pow.savefig('pow_spectrum_shift_superpos1.pdf')


plt.show()

