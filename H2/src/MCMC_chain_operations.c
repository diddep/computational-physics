#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


#include "tools.h"

#define NDIM 3


//easier to use R1, R2 as matrices or input x1,x2,y1,y2,z1,z2?

//R1, R2 on form [[x1,y1,z1],...,[xN,yN,zN]]
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

// double partialEnergyDerivative(double *E_local_derivative, double alpha, int N_steps, double **R1, double **R2){

//     double r12;
//     double *r1_nrm = malloc(sizeof(double) * NDIM), *r2_nrm = malloc(sizeof(double) * NDIM);
//     double *diff_vec = malloc(sizeof(double) * NDIM), *diff_nrm = malloc(sizeof(double) * NDIM);
//     //double *E_local_derivative;// = malloc(sizeof(double) * N_steps);
//     double prod = 0., div = 0., E_derivative_sum = 0;

//     for (int ix = 0; ix < N_steps; ++ix){
//         for (int dim=0; dim<NDIM; ++dim){
//             r1_nrm[dim] = R1[ix][dim];
//             r2_nrm[dim] = R2[ix][dim];
//         }

//         normalize_vector(r1_nrm, NDIM);
//         normalize_vector(r2_nrm, NDIM);
//         r12 = distance_between_vectors(R1[ix], R2[ix], NDIM);

//         elementwise_subtraction(diff_vec, R1[ix], R2[ix], NDIM);
//         elementwise_subtraction(diff_nrm, r1_nrm, r2_nrm, NDIM);
//         prod = dot_product(diff_vec, diff_nrm,NDIM);
//         div = (1. + alpha * r12);

//         E_local_derivative[ix] = - 0.0 \
//                                 - 2 * r12 * prod/(r12 * pow(div, 3.0)) \
//                                 + 3 * r12 * 1.0/(r12* pow(div, 4.0)) \
//                                 + 4 * r12 * 1./(4.0* pow(div, 5.0)) 
//                                 + 0.0;

//         E_derivative_sum += E_local_derivative[ix];
//     }

//     free(r1_nrm), free(r2_nrm), free(diff_nrm), free(diff_vec);

//     return E_derivative_sum/N_steps;
// }

//TODO kan göras snabbare om vi skickar in E_local chain istället
double partialEnergyDerivative(double *E_local_derivative, double alpha, int N_steps, double **R1, double **R2)
{
    double *E_local_chain = malloc(sizeof(double)*N_steps);
    double r12=0, ln_d_psi=0;
    double average_E=0, average_derivative =0, average_mix=0;
    double gradient=0;
    
    Energy(E_local_chain, alpha, N_steps, R1, R2);

    for(int step=0; step<N_steps; ++step)
    {
        r12 = distance_between_vectors(R1[step], R2[step], NDIM);
        ln_d_psi = - r12*r12/(2*pow((r12*alpha +1),2));
        average_derivative += ln_d_psi;
        average_mix += ln_d_psi*E_local_chain[step];
        average_E += E_local_chain[step];
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

/*theta(*R1, *R2) =  arccos( dot(R1,R2) / (|R1|*|R2|)*/
double theta_fun_vec(double *R1, double *R2)
{
    double R1_abs = vector_norm(R1,NDIM);
    double R2_abs = vector_norm(R2,NDIM);
    double R1_R2_dot = dot_product(R1, R2, NDIM);
    double theta = acos( R1_R2_dot / (R1_abs * R2_abs));

    return theta;
}
