import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


R1 =np.genfromtxt("../R1.csv", delimiter=',')
R2 =np.genfromtxt("../R2.csv", delimiter=',')

Energy_local = np.genfromtxt("../E_local.csv", delimiter=',')
Energy_local_derivative = np.genfromtxt("../E_local_derivative.csv", delimiter=',')
xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
theta_csv = np.genfromtxt("../theta.csv", delimiter=',')
theta = theta_csv[:,1]
Phi_k =np.genfromtxt("../phi_k.csv", delimiter=',')

alpha_results = np.genfromtxt("../alpha_results.csv", delimiter=',')


params = np.genfromtxt("../alpha_params.csv", delimiter=',')

n_alpha_steps = int(params[0])
n_discarded_steps = int(params[1])
alpha = params[2]
A = params[3]
beta = params[4]
n_steps = int(params[5])
d_displacement = params[6]
is_task1 = params[7]
is_task2 = params[8]
is_task3 = params[9]
is_task4 = params[10]

if is_task1:
    task_str = "task1"
elif is_task2:
    task_str = "task2"
elif is_task3:
    task_str = "task3"
elif is_task4:
    task_str = "task4"

steps = np.linspace(0, n_steps, n_steps)

R1_norm =  np.square(R1)
R1_norm = np.sqrt(np.sum(R1_norm, axis=1))
R2_norm =  np.square(R2)
R2_norm = np.sqrt(np.sum(R2_norm, axis=1))

def rho(rvec, z):
    rho = z**3 *4 * rvec**2 * np.exp(-2*z*rvec)

    return rho



## Plot histogram
n_bins = 70
fig_dist, ax_dist = plt.subplots(1,2, figsize=(10,5))

for idx, R_chain in enumerate([R1, R2]):
    R_norm =  np.square(R_chain)
    R_norm = np.sqrt(np.sum(R_norm, axis=1))
    
    rvec = np.linspace(0.1, np.max(R_norm))
    counts_r, bins_r = np.histogram(R_norm, bins = n_bins, density = True)
    
    ax_dist[idx].stairs(counts_r, bins_r, fill=True)
    ax_dist[idx].plot(rvec, rho(rvec, 27/16), color='r', linestyle='--',label= r'$\rho $ optimized', linewidth=3)
    ax_dist[idx].plot(rvec, rho(rvec, 2), color='k', linestyle=':',label= r'$\rho $ unscreened', linewidth=3)
    
    ax_dist[idx].set_title(f"Distribution for R$_{idx}$, alpha = {alpha}")#, fontsize=15)
    ax_dist[idx].set_xlabel(f"Radius [$a_0$]")#, fontsize=15)
    ax_dist[idx].set_ylabel("Probability density")#, fontsize=15)
    ax_dist[idx].legend()
plt.tight_layout()
fig_dist.savefig(f'plots_python/{task_str}/histogram_alpha{alpha}_nsteps{n_steps}.png')