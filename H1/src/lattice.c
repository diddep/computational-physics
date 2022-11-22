/*
H1lattice.c
Program that arranges atoms on a fcc lattice. 
Created by Anders Lindman on 2013-03-15.
*/

#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

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
    int i, j, k;
    int xor_value;
    
    double abs_displacement = 0.065;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed); 
    
    /*int n = 10;
    for (int ix = 0; ix < n; ix++){
	    //double u = gsl_rng_uniform(r);
        double displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
	    printf ("%.5f\n", displacement_factor);
    }*/
    
    double displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);

    double disp; srand(42);
    for(int ix = 0; ix <10; ix++){
        //disp = (2 * (rand() / RAND_MAX) - 1);// * 0.065;
        disp = (double)0.065*(rand() % 2000- 1000)/1000.0; //(2*rand()/RAND_MAX- 1) * 0.065;
        printf("%f\n", disp);
    }
    

    for (i = 0; i < 2 * N; i++){
        for (j = 0; j < 2 * N; j++){
            for (k = 0; k < N; k++){
                if (j % 2 == i % 2 ){
                    xor_value = 0;
                }
                else {
                    xor_value = 1;
                }
                displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                positions[i * N * 2 * N + j * N + k][0] = lattice_param * (0.5 * xor_value + k) + 0.;// displacement_factor * lattice_param;
                //positions[i * N * 2 * N + j * N + k][0] += displacement_factor * lattice_param;

                displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                positions[i * N * 2 * N + j * N + k][1] = lattice_param * (j * 0.5) + 0.; //displacement_factor * lattice_param;
                //positions[i * N * 2 * N + j * N + k][1] += displacement_factor * lattice_param;

                displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                positions[i * N * 2 * N + j * N + k][2] = lattice_param * (i * 0.5) + 0.; // displacement_factor * lattice_param;
                //positions[i * N * 2 * N + j * N + k][2] += displacement_factor * lattice_param;
                //printf ("%.5f\n", displacement_factor);
                
                /*for(int kx = 0; kx <3; kx++){
                    displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                    positions[i * N * 2 * N + j * N + k][kx] += displacement_factor * lattice_param;
                    printf ("%.5f\n", displacement_factor);
                }*/
                
                //positions[i * N * 2 * N + j * N + k][0] += displacement_factor * positions[i * N * 2 * N + j * N + k][0];
                //positions[i * N * 2 * N + j * N + k][1] *= 1; //displacement_factor;
                //positions[i * N * 2 * N + j * N + k][2] *= 1; //displacement_factor;

                
                //double displacement_factor = gsl_ran_flat(r, -abs_displacement, abs_displacement);
                //positions[i * N * 2 * N + j * N + k][1] += displacement_factor * lattice_param;
                //positions[i * N * 2 * N + j * N + k][2] += displacement_factor * lattice_param;
            
            }
        }
    }
    gsl_rng_free (r);
}