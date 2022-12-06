
import numpy as np
import matplotlib.pyplot as plt



str = "eq"
array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',')
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

t = array[-1:,0]
e_pot = array[-1:,2]
e_kin = array[-1:,3]
e_tot = array[-1:,4]
temp_eq = parameters[-1,5]


k_B= 8.61733 * 1e-5; N = 256



def C_v(epsilon_vec, T):

    variance = np.var(e_kin)
    term = 1-2*variance/(3*N*(k_B*T)**2)

    C_v = 3*N*k_B/term

    return

C_v_from_kinetic = C_v(e_kin,temp_eq)
C_v_from_potential = C_v(e_pot,temp_eq)

print("C_v from kinetic fluctuations=",C_v_from_kinetic,"\n-------------")
print("C_v from potential fluctuations=", C_v_from_potential,"\n--------")

