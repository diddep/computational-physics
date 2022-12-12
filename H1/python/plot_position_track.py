import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set the seaborn plotting style
sns.set()

str = "eq"
#str = "prod"

# load data from file
pos_array = np.genfromtxt(f'../csv/position_track_{str}.csv', delimiter=',', skip_header=1)
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')
pos_array_prod = np.genfromtxt(f'../csv/position_track_prod.csv', delimiter=',', skip_header=1)
parameters_prod = np.genfromtxt(f'../csv/parameters_prod.csv', delimiter=',')

pos_array = pos_array_prod

#t = pos_array[1:,0]
q1x = pos_array[:,1]
q1y = pos_array[:,2]
q1z = pos_array[:,3]
q2x = pos_array[:,4]
q2y = pos_array[:,5]
q2z = pos_array[:,6]
q3x = pos_array[:,7]
q3y = pos_array[:,8]
q3z = pos_array[:,9]

parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

end_time = parameters[-1,0]
dt = parameters[-1,1]

#dt = pos_array[0,0]
t = dt * np.linspace(0,len(pos_array[:,1]), len(pos_array[:,1]))

# Create a figure with three subplots
# Set up figure and axes
fig_pos, ax_pos = plt.subplots(1, 3, figsize=(12, 4))

# Iterate over dimensions
for i, dim in enumerate(["X", "Y", "Z"]):
    # Plot each dimension in a subplot
    ax_pos[i].plot(t, eval(f"q1{dim.lower()}"), label="q1")
    ax_pos[i].plot(t, eval(f"q2{dim.lower()}"), label="q2")
    ax_pos[i].plot(t, eval(f"q3{dim.lower()}"), label="q3")
    
    # Set x- and y-labels and title
    ax_pos[i].set_xlabel("Time (ps)", fontsize=15)
    ax_pos[i].set_ylabel(f"{dim}-coordinate (Å)", fontsize=15)
    ax_pos[i].set_title(f"{dim}-coordinate, dt={dt}", fontsize=15)
    
    # Add legend
    ax_pos[i].legend(fontsize=10)

# Adjust layout and save figure
plt.tight_layout()

str = ""
#str = "_calibration_melt"
plt.savefig(f"plots/position_track{str}.png")


