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


dtau = 0.02
#dtau=1

#final ET= 0.380387
#ET= 0.371350
#final ET= 0.377123
#final ET= 0.373763
#final ET= 0.379698
#final ET= 0.376354 200k walkers
#final ET= 0.357426

coordinate_array = np.genfromtxt("../../csv/task1_distribution.csv", delimiter=',')
ET_array= np.genfromtxt("../../csv/task1_ET_vec.csv", delimiter=',')
walker_array =  np.genfromtxt("../../csv/evolution_walkers.csv", delimiter=',')

print(np.max(coordinate_array), "max")

#coordinate_array =coordinate_array.values[0,:].astype(float)

#ET_array = ET_array.values[0,:].astype(float)
#walker_array = walker_array.values[0,:].astype(int)

iteration_array = np.arange(0, len(ET_array))*dtau


def wavefunction(xarray):

    psi0 = np.sqrt(2)*np.exp(- np.exp(-xarray) -xarray/2)

    return np.square(psi0)


cutoff =50000
xcoord = np.linspace(-3,18, 1000)
wavefunc = wavefunction(xcoord)

averageE = np.average(ET_array[-cutoff:])
std_E = np.std(ET_array[-cutoff:])
print(std_E, "std")

averageWalker = np.average(walker_array[-cutoff:])


print("average ET=", averageE)
fig_distribution, ax_distribution = plt.subplots(1,1)

n_bins =118

counts, bins = np.histogram(coordinate_array, bins = n_bins, density = False)

counts = np.square(counts)
counts = counts/np.linalg.norm(counts)
ax_distribution.stairs(counts, bins, fill = True, label='Sampled distribution')
ax_distribution.plot(xcoord, wavefunc,linewidth=3,color="r", label= r"Analytical distribution $|\Phi_0|^2$")
ax_distribution.set_title(R"Sampled distribution compared "+"\n"+ " to analytical distribution $|\Phi_0|^2$")
ax_distribution.set_xlabel(r"$x$ coordinate [$a_0$]")
ax_distribution.set_ylabel("Probability density")
ax_distribution.set_xlim(-5,18)

ax_distribution.legend()

fig_distribution.savefig("plots_task1/distribution.png")


fig_energy, ax_energy = plt.subplots(1,1)

ax_energy.plot(iteration_array, ET_array, label = r"Energy $E_T$")
ax_energy.set_title(r"Energy $E_T$ as a function of $\tau$")
ax_energy.set_ylabel(r"Energy [$E_h$]")
ax_energy.set_xlabel(r"$\tau$ [atomic units]")
ax_energy.axhline(averageE, color="r",linewidth =2, label =rf"Average $E_T ={averageE:.3f}$")
ax_energy.legend()

fig_energy.savefig("plots_task1/energy_evolution.png")

fig_walkers, ax_walkers = plt.subplots(1,1)

ax_walkers.plot(iteration_array, walker_array, label = r"Energy $E_T$")
ax_walkers.set_title(r"Number of walkers as function of $\tau$")
ax_walkers.set_ylabel(r"Number of walkers")
ax_walkers.set_xlabel(r"$\tau $ [atomic units]")
ax_walkers.axhline(averageWalker, color="r",linewidth =2, label =rf"Average number of walkers $ ={averageWalker:.0f}$")
ax_walkers.legend()

fig_walkers.savefig("plots_task1/walker_evolution.png")
        
fig_convergence , ax_convergence = plt.subplots(2,1, sharex=True)
dtau=1

ax_convergence[0].plot(iteration_array, ET_array,linewidth =2, label = r"Energy $E_T$")
ax_convergence[0].set_title(r"Energy $E_T$ as a function of number of iterations")
ax_convergence[0].set_ylabel(r"Energy [$E_h$]")
ax_convergence[0].set_xlabel(r"Iteration")
ax_convergence[0].axhline(averageE, color="r",linewidth =2, label =rf"Average $E_T ={averageE:.3f}$")
ax_convergence[0].legend()

ax_convergence[1].plot(iteration_array, walker_array,linewidth =2, label = r"Number of walkers")
ax_convergence[1].axhline(averageWalker,linewidth =2, color="r", label =rf"Average number of walkers ={averageWalker:.0f}$")

ax_convergence[1].set_title(r"Number of walkers as function of iteraion")
ax_convergence[1].set_ylabel(r"Number of walkers")
ax_convergence[1].set_xlabel(r"Iteration")
ax_convergence[1].set_yticks((150, 200,300, 400))
ax_convergence[1].legend()
plt.tight_layout()

fig_convergence.savefig("plots_task1/convergence.png")


        