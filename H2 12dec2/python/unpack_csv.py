import numpy as np
import pandas as pd
import time

def main():
    start_time = time.time()

    R1 = pd.read_csv(f"../csv/R1.csv", engine="pyarrow", names = ["x", "y", "z"])
    R2 = pd.read_csv("../csv/R2.csv", engine="pyarrow", names = ["x", "y", "z"])
    E_local = pd.read_csv("../csv/E_local.csv", engine="pyarrow", names= ["E_local"])
    E_local_derivative = pd.read_csv("../csv/E_local_derivative.csv", engine="pyarrow", names= ["E_local_derivative"])
    x_distribution = pd.read_csv("../csv/x_distribution.csv", engine="pyarrow", names= ["x_distribution"])
    #arr_theta_csv = pd.read_csv("../csv/theta.csv", engine="pyarrow", names= ["E_local"])
    
    phi_k =pd.read_csv("../csv/phi_k.csv", engine="pyarrow", names= ["phi_k"])
    alpha_results = pd.read_csv("../csv/alpha_results.csv", engine="pyarrow", names= ["ix", "E_average", "alpha", "gamma", "E_PD_average"])
    params = pd.read_csv("../csv/params.csv", names = ["n_alpha_steps", "N_discarded_steps", "alpha", "A", "beta", "N_steps", "d_displacement", "is_task1", "is_task2", "is_task3", "is_task4"])
    
    
    steps_linspace = np.linspace(0,int(params.N_steps), int(params.N_steps), endpoint=False)
    array_tuple = (R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, steps_linspace, params)

    # n_alpha_steps = int(params[0])
    # n_discarded_steps = int(params[1])
    # alpha = params[2]
    # A = params[3]
    # beta = params[4]
    # n_steps = int(params[5])
    # d_displacement = params[6]
    # is_task1 = params[7]
    # is_task2 = params[8]
    # is_task3 = params[9]
    # is_task4 = params[10]

    # param_tuple = (n_alpha_steps, n_discarded_steps, alpha, A, beta, n_steps, d_displacement, is_task1, is_task2, is_task3, is_task4)
    
    
    print("Read CSV:s with pyarrow in")
    print("--- %s seconds ---" % (time.time() - start_time))
    # TODO: Change array names
    # TODO: Change array name in plot task1, multiple names of arr_r1/R1
    # TODO: alpha-params -> params filename change


    return array_tuple

# results = main()

# (R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, params, steps_linspace) = results

