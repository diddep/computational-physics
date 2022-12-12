
import unpack_csv

results = unpack_csv.main()

(R1, R2, E_local, E_local_derivative, x_distribution, phi_k, steps_linspace, alpha_results, params) = results

print(R1.R1x)
print(E_local.shape)