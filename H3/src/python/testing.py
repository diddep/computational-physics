import numpy as np
import matplotlib.pyplot as plt

xcoord = np.linspace(-5,5, 200)
deltatau=0.02
ET=0.5


def weight(xcoord, delta_tau, ET):

    V = 0.5*(1- np.exp(-xcoord))**2

    W = np.exp(-(V-ET)*delta_tau)

    return W


Wvec = weight(xcoord, deltatau, ET)
Wvec2 = weight(xcoord, deltatau, ET/10)
Wvec3 = weight(xcoord, deltatau, ET=15)


fig, ax = plt.subplots(1,1)
ax.plot(xcoord, Wvec)
ax.plot(xcoord, Wvec2)
ax.plot(xcoord, Wvec3)

fig.savefig("weights.png")
