import numpy as np
import pandas as pd
import time
start_time = time.time()

print("Starting program")
print("--- %s seconds ---" % (time.time() - start_time))

#array_names = ["R1", "R2", "E_local", "E_local_derivative", "x_distribution", "theta", "phi_k", "alpha_results", "alpha_params"]
array_names = ["R1"]
for array_name in array_names:
    print(array_name)
    array_name =np.genfromtxt(f"../{array_name}.csv", delimiter=',')



#TODO: Change array name in plot task1, multiple names of arr_r1/R1
# R1 =np.genfromtxt("../R1.csv", delimiter=',')
# R2 =np.genfromtxt("../R2.csv", delimiter=',')
# arr_E_L = np.genfromtxt("../E_local.csv", delimiter=',')
# arr_E_L_Derivative = np.genfromtxt("../E_local_derivative.csv", delimiter=',')
# arr_xdist = np.genfromtxt("../x_distribution.csv", delimiter=',')
# arr_theta_csv = np.genfromtxt("../theta.csv", delimiter=',')
# Phi_k =np.genfromtxt("../phi_k.csv", delimiter=',')
# alpha_results = np.genfromtxt("../alpha_results.csv", delimiter=',')
# params = np.genfromtxt("../alpha_params.csv", delimiter=',')
# arr_theta = arr_theta_csv[:,1]
#TODO: arr_steps from alpha_params

#alpha_param_vector[] = {n_alpha_steps, N_discarded_steps, alpha, A, beta, N_steps, d_displacement, is_task1, is_task2, is_task3, is_task4};

print("Read with numpy")
print("--- %s seconds ---" % (time.time() - start_time))

#import pandas as pd
#R1 = pd.read_csv("../R1.csv", engine="pyarrow")

print("Read with pandas pyArrow")
print("--- %s seconds ---" % (time.time() - start_time))