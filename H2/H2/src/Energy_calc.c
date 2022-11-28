#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.c"


//easier to use R1, R2 as matrices or input x1,x2,y1,y2,z1,z2?

int Energy(double *E_L, double alpha, int N, double **R1, double **R2){

    double r12;
    double *r1_nrm = malloc(sizeof(double)*3);
    double *r2_nrm = malloc(sizeof(double)*3);
    double *diff_vec = malloc(sizeof(double)*3);
    double *diff_nrm = malloc(sizeof(double)*3);
    double prod=0.;
    double div = 0.;

    for (int i=0;i<N; ++i){

        for (int dim=0; dim<3; ++dim){
            r1_nrm[dim] = R1[i][dim];
            r2_nrm[dim] = R2[i][dim];
        }

        normalize_vector(r1_nrm);
        normalize_vector(r2_nrm);
        r12 = distance_between_vectors(R1, R2, 3);
        elementwise_subtraction( diff_vec, R1[i], R2[i], 3);
        elementwise_subtraction( diff_nrm, r1_nrm, r2_nrm, 3);
        prod = dot_product(diff_vec, diff_nrm);
        div = (1. +alpha*r12)

        E_L[i]=-4.0 + prod/(r12 * pow(div,2.0)) -1.0/(r12* pow(div,3.0)) - 1/(4.0* pow(div,4.0))+ 1/r12;

    }
    r12 = NULL;
    prod = NULL;
    div = NULL;
    free(r1_nrm);
    free(r2_nrm);
    free(diff_nrm)
    free(diff_vec);
}