#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
array = np.genfromtxt('../Task 1_1/signal.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()
ax.plot(array[:, 0], array[:, 1])
l_sig = len(array[:, 0])

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.set_title('dt=0.1/6')
ax.grid()

fig.savefig('signalsuperposedt16.pdf')