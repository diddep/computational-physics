#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


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

#fig.savefig('signalsuperposedt16.pdf')


array_pow = np.genfromtxt('../powerspectrum_shift_p1.csv', delimiter=',', skip_header=1)
array_pow2 = np.genfromtxt('../powerspectrum_shift_p2.csv', delimiter=',', skip_header=1)
array_pow3 = np.genfromtxt('../powerspectrum_shift_p3.csv', delimiter=',', skip_header=1)



fig_pow, ax_pow = plt.subplots()


#ax_pow.plot(array_pow[:,1], array_pow[:,0])
#ax_pow.plot(array_pow2[:,1], array_pow2[:,0])
ax_pow.plot(array_pow3[:,1], array_pow3[:,0]+ array_pow[:,0]+array_pow2[:,0])

print(array_pow3[:,0])

ax_pow.set_xlabel('frequency (arb.unit)')
ax_pow.set_ylabel('intensity (arb.unit)')
ax_pow.set_title('Power spectrum for co2')
ax_pow.set_xticks([-75,-39,39,75])
ax_pow.set_xlim(-150,150)

fig_pow.savefig('pow_spectrum_co2.pdf')


plt.show()

