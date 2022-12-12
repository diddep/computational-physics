import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# set default figure size
plt.rcParams["figure.figsize"] = [8, 6]

SMALL_SIZE = 15
MEDIUM_SIZE = 18
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 

for idx, val in enumerate(["eq", "prod"]):
    print(idx)
    print(val)
    str = val

    # load data from file
    array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',', skip_header=1)
    parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')
    pos_array = np.genfromtxt(f'../csv/position_track_{str}.csv', delimiter=',', skip_header=1)

    end_time = parameters[-1,0]
    dt = parameters[-1,1]
    lattice_param = parameters[-1,2]
    temp_scaling = parameters[-1,3]
    press_scaling = parameters[-1,4]
    temp_eq = parameters[-1,5]
    press_eq = parameters[-1,6]
    tau_T = parameters[-1,7]
    tau_P = parameters[-1,8]

    t = dt * np.linspace(0,len(array[:,1]), len(array[:,1]))
    cell_length = array[:,1]
    lattice_length = cell_length/4
    e_pot = array[:,2]
    e_kin = array[:,3]
    e_tot = array[:,4]
    temp = array[:,5]
    press = array[:,6]

    q1x = pos_array[:,1]
    q1y = pos_array[:,2]
    q1z = pos_array[:,3]
    q2x = pos_array[:,4]
    q2y = pos_array[:,5]
    q2z = pos_array[:,6]
    q3x = pos_array[:,7]
    q3y = pos_array[:,8]
    q3z = pos_array[:,9]

    # create figure and axes for energy plot
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(14, 5))

    # plot potential energy on first subplot
    axes[0].plot(t, e_pot, label='Potential energy', color="b")
    axes[0].set_title(f'Potential Energy, dt={dt} [ps]')

    # plot kinetic energy on second subplot
    axes[1].plot(t, e_kin, label='Kinetic energy', color="r")
    axes[1].set_title(f'Kinetic Energy, dt={dt} [ps]')

    # plot total energy on third subplot
    axes[2].plot(t, e_tot, label='Total energy', color="g")
    axes[2].set_title(f'Total Energy, dt={dt} [ps]')
    axes[2].ticklabel_format(useOffset=False)

    for idx in range(3):
        axes[idx].set_xlabel('Time [ps]')
        axes[idx].set_ylabel('Energy [eV/unit cell]')

    plt.tight_layout()
    plt.savefig(f'plots/energy_{str}.png')

    figlattice , ax_lattice = plt.subplots(1,1)

    ############# LATTICE ################
    # plot energy data
    ax_lattice.plot(t, lattice_length, label='Lattice parameter')

    average_lattice = np.mean(lattice_length)

    # add dashed line for the average energy
    ax_lattice.annotate(f'a$_{{0}}$ = {lattice_length[-1]:.3f}', xy=(t[-1], lattice_length[-1]), xytext=(-50, -50),
                textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

    # set labels and title
    ax_lattice.set_xlabel('Time [ps]')
    ax_lattice.set_ylabel('Lattice parameter [Å]')
    ax_lattice.set_title(f'Lattice parameter, dt={dt}')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/lattice_{str}.png')



    ####### POSITIONS ###################
    fig_pos, ax_pos = plt.subplots(1, 3, figsize=(14, 5))

    # Iterate over dimensions
    for i, dim in enumerate(["X", "Y", "Z"]):
        # Plot each dimension in a subplot
        ax_pos[i].plot(t, eval(f"q1{dim.lower()}"), label="q1")
        ax_pos[i].plot(t, eval(f"q2{dim.lower()}"), label="q2")
        ax_pos[i].plot(t, eval(f"q3{dim.lower()}"), label="q3")
        
        # Set x- and y-labels and title
        ax_pos[i].set_xlabel("Time [ps]")
        ax_pos[i].set_ylabel(f"{dim}-coordinate [Å]")
        ax_pos[i].set_title(f"{dim}-coordinate, dt={dt}")
        
        # Add legend
        ax_pos[i].legend(loc='lower right')

    # Adjust layout and save figure
    plt.tight_layout()
    plt.savefig(f"plots/position_track_{str}.png")


    ########### PRESSURE #############
    figP , axP = plt.subplots(1,1)
    axP.plot(t, press, label='Pressure')

    average_pressure = np.mean(press)

    # add dashed line for the average energy
    if(str == "prod"):
        print("Hej")
        axP.axhline(
            average_pressure,
            linestyle='--',
            color='r',
            linewidth=4,
            label=f'P$_{{average}}$ = {average_pressure:.2f} [Bar]'
        )
    
        #axP.annotate(f'P$_{{average}}$ = {average_pressure:.2f}', xy=(t[-1], average_pressure), xytext=(-150, 150),
        #            textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

    axP.set_xlabel('Time [ps]')
    axP.set_ylabel('Pressure [Bar]')
    axP.set_title(f'Pressure, dt={dt}')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/pressure_{str}.png')
    #plt.show()



    ########### TEMPTERATURE #############

    average_temperature = np.mean(array[:,5])

    fig2 , axT = plt.subplots(1,1)
    axT.plot(t, temp, label='Temperature')

    # add dashed line for the average energy
    if(str == "prod"):
        axT.axhline(
            average_temperature,
            linestyle='--',
            color='r',
            linewidth=4,
            label=f'T$_{{average}}$ = {average_temperature:.2f} [K]'
        )
    
        #axT.annotate(f'T$_{{average}}$ = {average_temperature:.2f}', xy=(t[-1], average_temperature), xytext=(-110, 200),
        #            textcoords='offset pixels', arrowprops=dict(arrowstyle='->', color='k'))

    plt.legend()
    axT.set_xlabel('Time [ps]')
    axT.set_ylabel('Temperature [K]')
    axT.set_title(f'Temperature, dt={dt}')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/temperature_{str}.png')





