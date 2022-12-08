
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()



# str = "eq"
 str = "prod"

array = np.genfromtxt(f'../csv/vel_verlet_{str}.csv', delimiter=',')
parameters = np.genfromtxt(f'../csv/parameters_{str}.csv', delimiter=',')

print(array.shape)
t = array[-1:,0]
#e_pot = array[-1:,2]
#e_kin = array[-1:,3]
#e_tot = array[-1:,4]

e_pot = array[:,2]
e_kin = array[:,3]
e_tot = array[:,4]


temp_eq = parameters[-1,5]

print(e_kin.shape, np.var(e_kin))
k_B= 8.61733 * 1e-5; N = 256



def C_v(epsilon_vec, T):

    variance = np.var(epsilon_vec)
    print(variance,'var')
    term = 1-2*variance/(3*N*(k_B*T)**2)

    C_v = 3*N*k_B/term

    return C_v

C_v_from_kinetic = C_v(e_kin,temp_eq)
C_v_from_potential = C_v(e_pot,temp_eq)

print("C_v from kinetic fluctuations=",C_v_from_kinetic,"\n-------------")
print("C_v from potential fluctuations=", C_v_from_potential,"\n--------")

