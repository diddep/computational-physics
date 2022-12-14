/*
H1lattice.c
Program that arranges atoms on a fcc lattice. 
Created by Anders Lindman on 2013-03-15.
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "tools.h"

/* Function takes a matrix of size [4*N*N*N][3] as input and stores a fcc lattice in it. N is the number of unit cells in each dimension and lattice_param is the lattice parameter. */
void init_fcc(double positions[][3], int N, double lattice_param)
{
    int i, j, k;
    int xor_value;
    
    for (i = 0; i < 2 * N; i++){
        for (j = 0; j < 2 * N; j++){
            for (k = 0; k < N; k++){
                if (j % 2 == i % 2 ){
                    xor_value = 0;
                }
                else {
                    xor_value = 1;
                }
                positions[i * N * 2 * N + j * N + k][0] = lattice_param * (0.5 * xor_value + k);
                positions[i * N * 2 * N + j * N + k][1] = lattice_param * (j * 0.5);
                positions[i * N * 2 * N + j * N + k][2] = lattice_param * (i * 0.5);
            }
        }
    }
}

void
displace_fcc(double positions[][3], int N, double lattice_param)
{
    int xor_value;
    double abs_displacement = 0.065, random_number; 

    //Initializing random number generator
    gsl_rng * r;
    r = init_random_num_generator();
    
    for (int ix = 0; ix < 2 * N; ix++){
        for (int jx = 0; jx < 2 * N; jx++){
            for (int kx = 0; kx < N; kx++){
                for(int lx = 0; lx < 3; lx++){
                    random_number = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                    positions[ix * N * 2 * N + jx * N + kx][lx] += random_number*lattice_param;
                }
            }
        }
    };
    gsl_rng_free(r);
}