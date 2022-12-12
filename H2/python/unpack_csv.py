import numpy as np
import pandas as pd
import time

def main():
    start_time = time.time()

    R1 = pd.read_csv(f"../csv/R1.csv", engine="pyarrow", header=0)
    R2 = pd.read_csv("../csv/R2.csv", engine="pyarrow", header=0)
    E_local = pd.read_csv("../csv/E_local.csv", engine="pyarrow", header=0)
    E_local_derivative = pd.read_csv("../csv/E_local_derivative.csv", engine="pyarrow", header=0)
    x_distribution = pd.read_csv("../csv/x_distribution.csv", engine="pyarrow", header=0)
    #arr_theta_csv = pd.read_csv("../csv/theta.csv", engine="pyarrow", header=0)
    phi_k =pd.read_csv("../csv/phi_k.csv", engine="pyarrow", header=0)
    alpha_results = pd.read_csv("../csv/alpha_results.csv", engine="pyarrow", header=0)
    params = np.genfromtxt("../csv/params.csv", delimiter=',')
    steps_linspace = np.linspace(0,int(params[5])+1, int(params[5])+1, endpoint=False)
  

    array_tuple = (R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, params, steps_linspace)

    print("Read CSV:s with pyarrow in")
    print("--- %s seconds ---" % (time.time() - start_time))
    # TODO: Change array names
    # TODO: Change array name in plot task1, multiple names of arr_r1/R1
    # TODO: alpha-params -> params filename change



    return array_tuple

# results = main()

# (R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, params, steps_linspace) = results

