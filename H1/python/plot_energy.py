import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# set default figure size
plt.rcParams["figure.figsize"] = [8, 6]

#str = "eq"
str = "prod"

# load data from file
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

end_time = parameters[-1,0]
dt = parameters[-1,1]
lattice_param = parameters[-1,2]
temp_scaling = parameters[-1,3]
press_scaling = parameters[-1,4]
temp_eq = parameters[-1,5]
press_eq = parameters[-1,6]
tau_T = parameters[-1,7]
tau_P = parameters[-1,8]

# extract columns from array
#dt = array[0,0]
t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))

e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]
temp = array[:,5]
press = array[:,6]

# create figure and axes for energy plot
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

# plot potential energy on first subplot
axes[0].plot(t, e_pot, label='Potential energy', color="b")
axes[0].set_title(f'Potential Energy, dt={dt} (ps)', fontsize = 15)

# plot kinetic energy on second subplot
axes[1].plot(t, e_kin, label='Kinetic energy', color="r")
axes[1].set_title(f'Kinetic Energy, dt={dt} (ps)', fontsize = 15)

# plot total energy on third subplot
axes[2].plot(t, e_tot, label='Total energy', color="g")
axes[2].set_title(f'Total Energy, dt={dt} (ps)', fontsize = 15)

for idx in range(3):
    axes[idx].set_xlabel('Time (ps)', fontsize = 15)
    axes[idx].set_ylabel('Potential Energy (eV/unit cell)', fontsize = 15)
    axes[idx].legend(fontsize=10)

# # create figure and axes for energy plot
# fig , ax = plt.subplots(1,1)

# # plot energy data
# ax.plot(t, e_pot, label='potential')
# ax.plot(t, e_kin, label='kinetic')
# ax.plot(t, e_tot, label='total')

# # set labels and title
# ax.set_xlabel('Time (ps)', fontsize = 15)
# ax.set_ylabel('Energy (eV/unit cell)', fontsize = 15)
# if(str == "eq"):
#     ax.set_title(f'Energy over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
# else:
#     ax.set_title(f'Energy over time with: dt={dt}', fontsize = 15)
plt.tight_layout()
plt.savefig(f'plots/energy.png')

# # calculate x and y coordinates for vertical arrow
# x1a = t[len(t)//5]
# y1a = e_tot[len(t)//5]
# y2a = (e_tot[0] - e_kin[0])/1.2

# x1b = t[len(t)-1]
# x2b = t[len(t)-200]
# y1b = e_tot[len(t)-1]
# y2b = (e_tot[0] - e_kin[0])/1.2

# # add vertical arrow to the plot
# ax.annotate(
#     f'E$_{{tot}}$ = {y1a:.6}', 
#     xy=(x1a, y1a),
#     xytext=(x1a, y2a), 
#     arrowprops=dict(arrowstyle='-|>', linewidth=2, color='k')
# )
# ax.annotate(
#     f'E$_{{tot}}$ = {y1b:.6}', 
#     xy=(x1b, y1b),
#     xytext=(x2b, y2b), 
#     arrowprops=dict(arrowstyle='-|>', linewidth=2, color='k')
# )



# # create figure and axes for energy plot
# fig , ax = plt.subplots(1,1)

# # plot energy data
# ax.plot(t, e_pot, label='potential')
# ax.plot(t, e_kin, label='kinetic')
# ax.plot(t, e_tot, label='total')

# # set labels and title
# ax.set_xlabel('Time (ps)', fontsize = 15)
# ax.set_ylabel('Energy (eV/unit cell)', fontsize = 15)
# if(str == "eq"):
#     ax.set_title(f'Energy over time: dt={dt}, $\tau_T$={tau_T}, $\tau_P$={tau_P}', fontsize = 15)
# else:
#     ax.set_title(f'Energy over time with: dt={dt}', fontsize = 15)
# plt.legend(fontsize=10)
# plt.tight_layout()
# plt.savefig(f'plots/energy.png')