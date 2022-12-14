import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

arr_Phi_k =np.genfromtxt("../phi_k.csv", delimiter=',')

print(arr_Phi_k)

plt.plot(np.arange(0,len(arr_Phi_k)), arr_Phi_k)

plt.savefig("plots_python/phi_k.png")