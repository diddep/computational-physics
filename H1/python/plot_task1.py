import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Read in data from file
array = np.genfromtxt('../csv/try_lattice_constants.csv', delimiter=',')
print(array[:, 2].shape)

# Define a function that represents the polynomial model
def polynomial(x, a, b, c):
    return a * x**2 + b * x + c
    #return d * x**3 + a * x**2 + b * x + c

# Fit the polynomial model to the data
params, covariance = curve_fit(polynomial, array[:, 2], array[:, 4])

# Create x values for the fitted curve
x_vals = np.linspace(np.min(array[:, 2]), np.max(array[:, 2]), 50)

# Get the y-coordinate of the smallest y-value for the model
y_min = np.min(polynomial(x_vals, *params))

# Find the index of the x-coordinate that corresponds to y_min
x_min_idx = np.argmin(polynomial(x_vals, *params))

# Get the x-coordinate that corresponds to y_min
x_min = x_vals[x_min_idx]


# Initialize figure and axes
fig, ax = plt.subplots(1, 1)

# Plot data points and fitted curve
ax.scatter(array[:, 2], array[:, 4], label='Data points')
ax.plot(x_vals, polynomial(x_vals, *params), linestyle=':', label=f'Curve fit, \n y = {params[0]:.2} x^2 + {params[1]:.2} x + {params[2]:.2}')

# Add an arrow pointing to the smallest y-value for the model
ax.annotate(f'a$_{{min}}$ = {np.cbrt(x_min):.2f}', xy=(x_min, y_min), xytext=(50, -20),
            textcoords='offset pixels', arrowprops=dict(arrowstyle='->'))
            
# Set axes labels and font sizes
ax.set_ylabel(f'Volume Å$^3$', fontsize=15)
ax.set_xlabel(f'Energy (eV/unit cell)', fontsize=15)

# Add legend
plt.legend(fontsize=10)

# Adjust figure layout and save
plt.tight_layout()
plt.savefig('plots/task1_plot.png')



# import numpy as np
# import matplotlib.pyplot as plt

# # Read in data from file
# array = np.genfromtxt('../csv/try_lattice_constants.csv', delimiter=',')

# # Fit a polynomial model to the data
# model = np.poly1d(np.polyfit(array[:, 2], array[:, 4], 2))

# # Create x values for the fitted curve
# x_vals = np.linspace(np.min(array[:, 2]), np.max(array[:, 2]), 50)

# # Initialize figure and axes
# fig, ax = plt.subplots(1, 1)

# # Plot data points and fitted curve
# ax.scatter(array[:, 2], array[:, 4], label='Data points')
# ax.plot(x_vals, model(x_vals), linestyle=':', label=f'Curve fit, \n y = {model[2]:.2} x^2 + {model[1]:.2} x + {model[0]:.2}')

# # Set axes labels and font sizes
# ax.set_ylabel('Volume Å^3', fontsize=15)
# ax.set_xlabel('Energy (eV/unit cell)', fontsize=15)

# # Add legend
# plt.legend(fontsize=10)

# # Adjust figure layout and save
# plt.tight_layout()
# plt.savefig('plots/task1_plot.png')
