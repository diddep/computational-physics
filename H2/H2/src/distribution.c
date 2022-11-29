//
// Created by didri on 2022-11-28.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


#include "tools.h"




double distribution(double *R1, double *R2, double alpha){

    double r1 = vector_norm(R1,3);
    double r2 = vector_norm(R2, 3);
    double r12 = distance_between_vectors(R1, R2, 3);
    double psi = exp(-2.0*r1) * exp(-2.0*r2) *exp(r12/(2.0*(1.0+alpha)));
    double w = psi*psi;

    return w;
}
