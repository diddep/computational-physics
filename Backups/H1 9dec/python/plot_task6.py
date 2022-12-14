import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

Cell_length = 4*4.03; number_particles = 256

array = np.genfromtxt(f'../csv/radial_distribution.csv', delimiter=',', skip_header=0)

number_of_bins = len(array)
print(number_of_bins)

def N_ideal(bin_arr, Cell_length, Number_of_particles):
    Number_of_bins = len(bin_arr)

    delta_r = Cell_length/Number_of_bins
    V= Cell_length**3

    N_ideal = (Number_of_particles-1)*4*np.pi/(3*Cell_length**3) * (np.square(bin_arr) -3*bin_arr+1)*delta_r**3
    
    return N_ideal


bin_array = np.arange(number_of_bins) 

ideal = N_ideal(bin_array, Cell_length, number_particles)

delta_r = Cell_length/number_particles
fig, ax_radial = plt.subplots(1,1)
array = (array/ideal)

max =np.max(array)
conv = 1#3/max

ax_radial.plot(bin_array*delta_r, array*conv)
ax_radial.set_xlim(0,4)

minima = array[array > 0].min()

# Annotate the plot with an arrow pointing to the minimum
ax_radial.annotate(
    'First non-zero minimum',  # text to display
    xy=(minima, 0),           # coordinates of the arrow's tail
    xytext=(minima, 1),       # coordinates of the arrow's head
    arrowprops=dict(           # arrow style
        facecolor='black',
        shrink=0.05,
        width=2,
        headwidth=8
    )
)

fig.savefig(f'plots/radial.png')