import numpy as np
import pandas as pd
import time

def main():
    start_time = time.time()

    R1 = pd.read_csv(f"../csv/R1.csv", engine="pyarrow", names = ["R1x", "R1y", "R1z"])
    R2 = pd.read_csv("../csv/R2.csv", engine="pyarrow", names = ["R1x", "R1y", "R1z"])
    E_local = pd.read_csv("../csv/E_local.csv", engine="pyarrow", names= ["E_local"])
    E_local_derivative = pd.read_csv("../csv/E_local_derivative.csv", engine="pyarrow", names= ["E_local_derivative"])
    x_distribution = pd.read_csv("../csv/x_distribution.csv", engine="pyarrow", names= ["x_distribution"])
    theta_distribution = pd.read_csv("../csv/theta.csv", engine="pyarrow", names= ["theta"])
    phi_k =pd.read_csv("../csv/phi_k.csv", engine="pyarrow", names= ["phi_k"])
    alpha_results = pd.read_csv("../csv/alpha_results.csv", engine="pyarrow", names= ["ix", "E_average", "alpha", "gamma", "E_PD_average", "beta"])
    params = pd.read_csv("../csv/params.csv", names = ["N_alpha_steps", "N_discarded_steps", "alpha", "A", "beta", "N_steps", "d_displacement", "is_task1", "is_task2", "is_task3", "is_task4"])
    
    steps_linspace = np.linspace(0,int(params.N_steps), int(params.N_steps), endpoint=False)

    E_local=E_local.values[:,0].astype(float)
    E_local_derivative=E_local_derivative.values[:,0].astype(float)
    x_distribution=x_distribution.values[:,0].astype(float)
    theta_distribution = theta_distribution.values[:,0].astype(float)
    phi_k=phi_k.values[:,0].astype(float)

    array_tuple = (R1, R2, E_local, E_local_derivative, x_distribution, theta_distribution, phi_k, steps_linspace, alpha_results, params)


    return array_tuple