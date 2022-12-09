import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
sns.set()

Cell_length = 4*4.03; number_particles = 256

array = np.genfromtxt(f'../csv/radial_distribution.csv', delimiter=',', skip_header=0)

number_of_bins = len(array)

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

max =np.argmax(array)
conv = 1/1.8#3/max

ax_radial.plot(2*bin_array*delta_r, array*conv)
ax_radial.set_xlim(0,8)



max = np.argmax(array)
print(max)
minima = np.argmin(array[max:])
print(minima*delta_r)
#print(array[15:80])

ax_radial.axvline(minima*delta_r)

print(array[:minima])

integrand =[]
#34

for i,r in enumerate(2*bin_array*delta_r):
    print(r,i, array[i])
    integrand.append(r**2* 4*np.pi* array[i]*conv)

#print(integrand)

#coordination_number = sp.integrate.trapezoid(2*bin_array[max]*delta_r, integrand)

#print(coordination_number)

# Annotate the plot with an arrow pointing to the minimum


fig.savefig(f'plots/radial.png')