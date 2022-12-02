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


## Plot xdistribution
# fig_xdist, ax_xdist = plt.subplots(1,1)
# counts_xdist, bins_xdist = np.histogram(xdist, bins = n_bins, density = True)
# ax_xdist.stairs(counts_xdist, bins_xdist)
# fig_xdist.savefig(f'plots_python/{task_str}/x_distribution_alpha{alpha}_nsteps{n_steps}.png')


n_bins = 70

## Plot theta
fig_theta, ax_theta = plt.subplots(1,2, figsize=(10,5))
counts_x, bins_x = np.histogram(theta, bins = n_bins, density = True)
counts_theta, bins_theta = np.histogram(np.arccos(theta), bins = n_bins, density = True)

ax_theta[0].stairs(counts_x, bins_x, fill = True, label='Sampled distribution')
ax_theta[0].axhline(1/2, label='Distribution for uncorrelated electrons', color='r', linewidth=2, linestyle='dashed')
ax_theta[0].set_title('Distribution for x', fontsize=15)
ax_theta[0].set_xlabel('x', fontsize=15)
ax_theta[0].set_ylabel('Probability density', fontsize=15)
ax_theta[0].legend()

Fill = False
ax_theta[1].stairs(counts_theta, bins_theta, fill = True)
ax_theta[1].set_title(r'Distribution for $\theta$', fontsize=15)
ax_theta[1].set_xlabel(r'$\theta$ [rad]', fontsize=15)
ax_theta[1].set_ylabel('Probability density', fontsize=15)
plt.tight_layout()
fig_theta.savefig(f'plots_python/{task_str}/Distributions_alpha{alpha}_nsteps{n_steps}.png')


## Plot phi_k
#plt.plot(np.arange(0,len(Phi_k)), Phi_k)
#plt.savefig(f"plots_python/{task_str}/phi_k_alpha{alpha}_nsteps{n_steps}.png")