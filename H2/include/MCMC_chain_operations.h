//
// Created by didri on 2022-11-28.
//

void
Energy(
        double *E_local,
        double alpha,
        int N,
        double **R1,
        double **R2
        );
        
void
partialEnergyDerivative(
        double alpha,
        int N,
        double **R1,
        double **R2
        );
void
x_distribution(
        double *x_chain, int N, double **R1_chain, double **R2_chain
);

void theta_fun(double *theta_chain, int N_steps, double **R1_chain, double **R2_chain);

double theta_fun_vec(double *R1, double *R2);