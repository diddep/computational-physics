import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set the seaborn plotting style
sns.set()

array = np.genfromtxt('../position_track_eq.csv', delimiter=',')

t = array[:,0]
q1x = array[:,1]
q1y = array[:,2]
q1z = array[:,3]
q2x = array[:,4]
q2y = array[:,5]
q2z = array[:,6]
q3x = array[:,7]
q3y = array[:,8]
q3z = array[:,9]

dt = t[0]


# Create a figure with three subplots
# Set up figure and axes
fig, ax = plt.subplots(1, 3, figsize=(12, 4))

# Iterate over dimensions
for i, dim in enumerate(["X", "Y", "Z"]):
    # Plot each dimension in a subplot
    ax[i].plot(t, eval(f"q1{dim.lower()}"), label="q1")
    ax[i].plot(t, eval(f"q2{dim.lower()}"), label="q2")
    ax[i].plot(t, eval(f"q3{dim.lower()}"), label="q3")
    
    # Set x- and y-labels and title
    ax[i].set_xlabel("Time (ps)", fontsize=15)
    ax[i].set_ylabel(f"{dim}-coordinate (Å)", fontsize=15)
    ax[i].set_title(f"{dim}-coordinate, dt={dt}", fontsize=15)
    
    # Add legend
    ax[i].legend(fontsize=10)

# Adjust layout and save figure
plt.tight_layout()

str = ""
#str = "_calibration_melt"
plt.savefig(f"plots/position_track{str}.png")


