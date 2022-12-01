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
void Energy(double *E_L, double alpha, int N_steps, double **R1, double **R2){

    double r12;
    double *r1_nrm = malloc(sizeof(double) * NDIM), *r2_nrm = malloc(sizeof(double) * NDIM);
    double *diff_vec = malloc(sizeof(double) * NDIM), *diff_nrm = malloc(sizeof(double) * NDIM);
    double prod = 0., div = 0.;

    for (int i = 0; i < N_steps; ++i){
        for (int dim=0; dim<NDIM; ++dim){
            r1_nrm[dim] = R1[i][dim];
            r2_nrm[dim] = R2[i][dim];
        }

        normalize_vector(r1_nrm, NDIM);
        normalize_vector(r2_nrm, NDIM);
        r12 = distance_between_vectors(R1[i], R2[i], NDIM);

        elementwise_subtraction(diff_vec, R1[i], R2[i], NDIM);
        elementwise_subtraction(diff_nrm, r1_nrm, r2_nrm, NDIM);
        prod = dot_product(diff_vec, diff_nrm,NDIM);
        div = (1. + alpha * r12);

        E_L[i] = - 4.0 + prod/(r12 * pow(div, 2.0)) -1.0/(r12* pow(div, NDIM))\
                 - 1./(4.0* pow(div,4.0))+ 1. / r12;

    }
    //r12 = NULL;
    //prod = NULL;
    //div = NULL;
    free(r1_nrm);
    free(r2_nrm);
    free(diff_nrm);
    free(diff_vec);
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

/*theta(*R1, *R2) =  arccos( dot(R1,R2) / (|R1|*|R2|)*/
double theta_fun(double *R1, double *R2)
{
    double R1_abs = vector_norm(R1,NDIM);
    double R2_abs = vector_norm(R2,NDIM);
    double R1_R2_dot = dot_product(R1, R2, NDIM);
    double theta = acos( R1_R2_dot / (R1_abs * R2_abs));

    return theta;
}