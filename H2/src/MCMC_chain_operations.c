#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


#include "tools.h"
//number of spacial dimensions
#define NDIM 3



/*Function that calculates energy at each step of MCMC-chain from two position vectors at each step
    args:
        double *E_local = N_steps long array to save values of energy in
        double alpha = parameter value 
        int N_steps = number of steps in mcmc-chain
        double **R1 = 2D-array [3][N_steps] that saves x,y,z for electron 1 in each step of mcmc-chain
        double **R2 = 2D-array [3][N_steps] that saves x,y,z for electron 2 in each step of mcmc-chain
*/

void Energy(double *E_local, double alpha, int N_steps, double **R1, double **R2){

    double r12;
    double *r1_nrm = malloc(sizeof(double) * NDIM), *r2_nrm = malloc(sizeof(double) * NDIM);
    double *diff_vec = malloc(sizeof(double) * NDIM), *diff_nrm = malloc(sizeof(double) * NDIM);
    double prod = 0., div = 0.;

    for (int ix = 0; ix < N_steps; ++ix){
        for (int dim=0; dim<NDIM; ++dim){
            r1_nrm[dim] = R1[ix][dim];
            r2_nrm[dim] = R2[ix][dim];
        }

        normalize_vector(r1_nrm, NDIM);
        normalize_vector(r2_nrm, NDIM);
        r12 = distance_between_vectors(R1[ix], R2[ix], NDIM);

        elementwise_subtraction(diff_vec, R1[ix], R2[ix], NDIM);
        elementwise_subtraction(diff_nrm, r1_nrm, r2_nrm, NDIM);
        prod = dot_product(diff_vec, diff_nrm,NDIM);
        div = (1. + alpha * r12);

        E_local[ix] = - 4.0 \
                     + prod/(r12 * pow(div, 2.0)) \
                     - 1.0/(r12* pow(div, 3.0)) \
                     - 1./(4.0* pow(div,4.0)) \
                     + 1. / r12;

    }
    free(r1_nrm), free(r2_nrm), free(diff_nrm), free(diff_vec);
}

/*Function that calculates gradient with respect to parameter alpha
for performing damped gradient descent.
    args:
        double *E_local_derivative = N_steps long array to save values of derivative
        double alpha = parameter value 
        int N_steps = number of steps in mcmc-chain
        double **R1 = 2D-array [3][N_steps] that saves x,y,z for electron 1 in each step of mcmc-chain
        double **R2 = 2D-array [3][N_steps] that saves x,y,z for electron 2 in each step of mcmc-chain
    returns:
       double gradient = Gradient for performing damped gradient descent
*/


double partialEnergyDerivative(double *E_local_derivative, double alpha, int N_steps, double **R1, double **R2)
{
    double *E_local_chain = malloc(sizeof(double)*N_steps);
    double r12=0, ln_d_psi=0;
    double average_E=0, average_derivative =0, average_mix=0;
    double gradient=0, gradient_step = 0;
    
    Energy(E_local_chain, alpha, N_steps, R1, R2);

    for(int step=0; step<N_steps; ++step)
    {
        r12 = distance_between_vectors(R1[step], R2[step], NDIM);
        ln_d_psi = - r12*r12/(2*pow((r12*alpha +1),2));
        average_derivative += ln_d_psi;
        average_mix += ln_d_psi*E_local_chain[step];
        average_E += E_local_chain[step];
        gradient_step = E_local_chain[step] * ln_d_psi;
        E_local_derivative[step] = gradient_step;
    }
    average_derivative /= N_steps;
    average_mix /= N_steps;
    average_E /= N_steps; 

    gradient = 2*(average_mix - average_E*average_derivative);
    
    free(E_local_chain);

    return gradient;
}


/*function calculating distribution of x from two position mcmc chains
 * args:
 *      x_chain = array of length N_steps to save distribution in
 *      N_steps = number of steps in mcmc chain
 *      R1_chain = mcmc chain for particle 1, NX3 matrix
 *      R2_chain = mcmc chain for particle 1, NX3 matrix
 */
void x_distribution(double *x_chain, int N_steps, double **R1_chain, double **R2_chain){

    // initializing values used to calculate instance of x
    double length_r1=0, length_r2=0, dot_prod=0;

    // stepping through mcmc chain calculating x at every step and saving in x_chain
    for(int step=0; step<N_steps; ++step)
    {
        length_r1 = vector_norm(R1_chain[step], NDIM);
        length_r2 = vector_norm(R2_chain[step], NDIM);
        dot_prod = dot_product(R1_chain[step], R1_chain[step], NDIM);
        x_chain[step] = dot_prod/(length_r2*length_r1);
    }
}

/*
function that calculates the angle between arrays along entire mcmc-chain
 args:
        double *theta_chain = N_steps long array to save values of angle in
        double alpha = parameter value 
        int N_steps = number of steps in mcmc-chain
        double **R1 = 2D-array [3][N_steps] that saves x,y,z for electron 1 in each step of mcmc-chain
        double **R2 = 2D-array [3][N_steps] that saves x,y,z for electron 2 in each step of mcmc-chain
*/
void theta_fun(double *theta_chain, int N_steps, double **R1_chain, double **R2_chain)
{
    // initializing values used to calculate instance of x
    double length_R1=0, length_R2=0, dot_prod=0;
    for(int step = 0; step < N_steps; ++step)
    {  
        length_R1 = vector_norm(R1_chain[step], NDIM);
        length_R2 = vector_norm(R2_chain[step], NDIM);
        dot_prod = dot_product(R1_chain[step], R1_chain[step], NDIM);
        theta_chain[step] = acos( dot_prod / (length_R1 * length_R2));
    }
}

/*
function that calculates angle between position vectors for electron
    args:
        double *R1 = position array for electron 1 [x,y,z]
        double *R2 = position array for electron 2 [x,y,z]
    returns:
        double theta = angle between the position vectors of the two electrons
*/
double theta_fun_vec(double *R1, double *R2)
{
    double R1_abs = vector_norm(R1,NDIM);
    double R2_abs = vector_norm(R2,NDIM);
    double R1_R2_dot = dot_product(R1, R2, NDIM);
    double theta = acos( R1_R2_dot / (R1_abs * R2_abs));

    return theta;
}
