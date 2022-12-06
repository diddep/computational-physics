import numpy as np
import matplotlib.pyplot as plt

# Read in data from file
array = np.genfromtxt('../csv/try_lattice_constants.csv', delimiter=',')

# Fit a polynomial model to the data
model = np.poly1d(np.polyfit(array[:, 2], array[:, 4], 2))

# Create x values for the fitted curve
x_vals = np.linspace(np.min(array[:, 2]), np.max(array[:, 2]), 50)

# Initialize figure and axes
fig, ax = plt.subplots(1, 1)

# Plot data points and fitted curve
ax.scatter(array[:, 2], array[:, 4], label='Data points')
ax.plot(x_vals, model(x_vals), linestyle=':', label=f'Curve fit, \n y = {model[2]:.2} x^2 + {model[1]:.2} x + {model[0]:.2}')

# Set axes labels and font sizes
ax.set_ylabel('Volume Å^3', fontsize=15)
ax.set_xlabel('Energy (eV/unit cell)', fontsize=15)

# Add legend
plt.legend(fontsize=10)

# Adjust figure layout and save
plt.tight_layout()
plt.savefig('plots/task1_plot.png')
