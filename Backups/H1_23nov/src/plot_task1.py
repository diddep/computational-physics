import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('../try_lattice_constants.csv', delimiter=',').reshape(-1,5)


fig_vol_en , ax = plt.subplots(1,1)

model = np.poly1d(np.polyfit(array[:,2],array[:,4], 2))
print(model)

lin = np.linspace(np.min(array[:,2]), np.max(array[:,2]), 50)
ax.scatter(array[:,2],array[:,4], label='datapoints')
ax.plot(lin, model(lin), linestyle=':', label=f'Curvefit, \n y ={model[2]:.2} x^2+ {model[1]:.2} x +{model[0]:.2}')
ax.set_ylabel('Volume Å^3', fontsize = 15)
ax.set_xlabel('Energy (eV/unit cell)', fontsize = 15)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('task1_plot.pdf')


plt.show()

