import unpack_csv

def main():

    results = unpack_csv.main()
    (R1, R2, E_local, E_local_derivative, x_distribution, phi_k, alpha_results, steps_linspace, params) = results

    if params.is_task1.values[0]:
        task_str = "task1"
    elif params.is_task2.values[0]:
        task_str = "task2"
    elif params.is_task3.values[0]:
        task_str = "task3"
    elif params.is_task4.values[0]:
        task_str = "task4"

    return task_str