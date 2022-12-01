
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns

#sns.set_theme()


arr_R1 =np.genfromtxt("../R1.csv", delimiter=',')
arr_R2 =np.genfromtxt("../R2.csv", delimiter=',')

#print(arr_R2.shape)
arr_E_L = np.genfromtxt("../E_L.csv", delimiter=',')
arr_xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
print(np.max(arr_xdist),'xd')

arr_r1 = np.square(arr_R1)
arr_r1 = np.sqrt(np.sum(arr_r1, axis=1))
arr_r2 = np.square(arr_R2)
arr_r2 = np.sqrt(np.sum(arr_r2, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho

def theta(R1,R2):
    n,d = np.shape(R1)
    arg = np.zeros(n)

    for st in range(0,n):

        r1 = R1[st,:]; r2 = R1[st,:]
        #print(np.dot(r1,r2)/ (np.linalg.norm(r1)* np.linalg.norm(r2)))
        arg[st] = np.dot(r1,r2) / (np.linalg.norm(r1)* np.linalg.norm(r2))
    return arg

n_bins = 1000
fig_dist, ax_dist = plt.subplots(1,2)
counts_r1, bins_r1 = np.histogram(arr_r1, bins = n_bins, density = True)
counts_r2, bins_r2 = np.histogram(arr_r2, bins = n_bins, density = True)
rvec = np.linspace(0.1, np.max(arr_r2))
#print(rvec)

ax_dist[0].stairs(counts_r1, bins_r1)
ax_dist[0].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--')
ax_dist[0].set_xlabel("Radius")
ax_dist[0].set_ylabel("Counts")
ax_dist[1].stairs(counts_r2, bins_r2)
ax_dist[1].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--')
ax_dist[1].set_xlabel("Radius")
ax_dist[1].set_ylabel("Counts")

fig_dist.savefig('hist.png')


fig_energy, ax_energy = plt.subplots(1,1)
ax_energy.scatter(arr_r1, arr_E_L)
fig_energy.savefig('energy_r1.png')


#for X in x:
#    print(X)

fig_dist, ax_dist = plt.subplots(1,1)
#counts_theta, bins_theta = np.histogram(x, bins = 10, density = True)
#ax_theta.stairs(counts_theta, bins_theta)
counts_xdist, bins_xdist = np.histogram(arr_xdist, bins = n_bins, density = True)
ax_dist.stairs(counts_xdist, bins_xdist)


fig_dist.savefig('x_distribution.png')




