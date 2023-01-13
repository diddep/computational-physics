import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme()


plt.rcParams["figure.figsize"] = [8, 6]

SMALL_SIZE = 15
MEDIUM_SIZE = 18
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 


#coordinate_array= pd.read_csv("../../csv/task1_distribution.csv", engine="pyarrow", names= ["coordinates"])
#ET_array= pd.read_csv("../../csv/task1_ET_vec.csv", engine="pyarrow", names= ["energy"])
#walker_array= pd.read_csv("../../csv/evolution_walkers.csv", engine="pyarrow", names= ["energy"])


dtau = 0.01

#final ET= -2.920467

ET_array= np.genfromtxt("../../csv/task2_ET_vec.csv", delimiter=',')
walker_array =  np.genfromtxt("../../csv/task2_evolution_walkers.csv", delimiter=',')


#coordinate_array =coordinate_array.values[0,:].astype(float)

#ET_array = ET_array.values[0,:].astype(float)
#walker_array = walker_array.values[0,:].astype(int)

iteration_array = np.arange(0, len(ET_array))*dtau


def wavefunction(xarray):

    psi0 = np.sqrt(2)*np.exp(- np.exp(-xarray) -xarray/2)

    return np.square(psi0)

cutoff = 10000

averageE = np.average(ET_array[-cutoff:])
avvec = np.zeros(len(ET_array)-1000)
std_E = np.std(ET_array)
print("std E", std_E)

#for i in range(1000, len(ET_array)):

 #   avvec[i] = np.average(ET_array[-i:])

print("average ET=", averageE)
fig_distribution, ax_distribution = plt.subplots(1,1)


fig_energy, ax_energy = plt.subplots(1,1)

ax_energy.plot(iteration_array, ET_array, label = r"Energy $E_T$")
ax_energy.set_title(r"Energy $E_T$ as a function of $\tau$")
ax_energy.set_ylabel(r"Energy [$E_h$]")
ax_energy.set_xlabel(r"$\tau$")
ax_energy.axhline(averageE, color="r", label ="Average of $E_T =-2.901$")
ax_energy.legend()

fig_energy.savefig("plots_task2/energy_evolution.png")

fig_walkers, ax_walkers = plt.subplots(1,1)

ax_walkers.plot(iteration_array, walker_array, label = r"Number of walkers")
ax_walkers.set_title(r"Number of walkers as function of $\tau$")
ax_walkers.set_ylabel(r"Number of walkers")
ax_walkers.set_xlabel(r"$\tau$")
ax_walkers.legend()

fig_convergence , ax_convergence = plt.subplots(1,2)
dtau=1

ax_convergence[0].plot(iteration_array, ET_array, label = r"Energy $E_T$")
ax_convergence[0].set_title(r"Energy $E_T$ as a function of number of iterations")
ax_convergence[0].set_ylabel(r"Energy [$E_h$]")
ax_convergence[0].set_xlabel(r"Iteration")
ax_convergence[0].axhline(averageE, color="r", label ="Average of $E_T =-2.901$")
ax_convergence[0].legend()

ax_convergence[1].plot(iteration_array, walker_array, label = r"Number of walkers")
ax_convergence[1].set_title(r"Number of walkers as function of iteraion")
ax_convergence[1].set_ylabel(r"Number of walkers")
ax_convergence[1].set_xlabel(r"$Iteration")
ax_convergence[1].legend()


fig_convergence.savefig("plots_task2/convergence.png")
        


        