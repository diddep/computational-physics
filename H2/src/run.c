#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include <stdbool.h>

#include "tools.h"
#include "distribution.h"
#include "Energy_calc.h"

#define NDIM 3

// Function running the markov-chain
int
run(
    int argc,
    char *argv[]
   )
{
    
    // MCMC Parameters
    // N_steps is number of steps in loop, alpha --- , d_displacement offsets particle positions
    int N_steps = 1e5; double alpha = 0.1, d_displacement = 0.1; 

    // Filenames for saving in csv
    char filename_R1[] = {"R1.csv"}, filename_R2[] = {"R2.csv"};
    char filename_energy[] = {"E_L.csv"}, filename_xdist[] = {"x_distribution.csv"};
    char filename_results[] = {"MCMC_results.csv"}, filename_params[] = {"MCMC_params.csv"};

     // Initializing arrays
    double **R1 = create_2D_array(N_steps,NDIM), **R2 = create_2D_array(N_steps,NDIM);
    double *E_L = malloc(sizeof(double) * N_steps);
    double *x_chain = malloc(sizeof(double) * N_steps);

    double R1_test[NDIM], R2_test[NDIM];
    double random_number1 = 0, random_number2 = 0;
    gsl_rng * r;
    r = init_random_num_generator();

    // Initializing random starting values for all particles
    printf("Initializing random positions:\n");
    for (int kx = 0; kx < NDIM; kx++)
    {
        // -0.5 to 0.5 because it's length 1 and symmetric around zero
        random_number1 = gsl_ran_flat(r, -0.5, 0.5);
        random_number2 = gsl_ran_flat(r, -0.5, 0.5);

        R1[0][kx] = d_displacement * random_number1;
        R2[0][kx] = d_displacement * random_number2;
        printf("R1[0][%d]: %f, R2[0][%d]: %f\n", kx, R1[0][kx], kx, R2[0][kx]);
    }

    int accept_count = 0;
    for(int ix = 0; ix < N_steps - 1; ++ix)
    {
        // Get proposal position, d_displacement step size
        for (int kx = 0; kx < NDIM; ++kx)
        {
            random_number1 = gsl_ran_flat(r,-0.5,0.5);
            random_number2 = gsl_ran_flat(r,-0.5,0.5);

            R1_test[kx] = R1[ix][kx] + d_displacement * random_number1;
            R2_test[kx] = R2[ix][kx] + d_displacement * random_number2;
        }

        // Probability for particle occupying new and old positions
        double prob_test = distribution(R1_test, R2_test, alpha);
        double prob_old = distribution(R1[ix], R2[ix], alpha);

        // If new prob > old, make step OR take exploration step
        if(prob_test > prob_old || prob_test / prob_old > gsl_ran_flat(r,0.0,1.0))
        {
            // Counter for how many accepts
            accept_count = accept_count + 1;

            // If accepted, save new position in next row.
            for (int kx=0; kx<NDIM; ++kx)
            {
                R1[ix+1][kx] = R1_test[kx];
                R2[ix+1][kx] = R2_test[kx];
            }
        }
        // If not accepted save old position in next row
        else{
            for (int kx=0; kx<NDIM; ++kx)
            {
                R1[ix+1][kx] = R1[ix][kx];
                R2[ix+1][kx] = R2[ix][kx];
            }
        }
    }

    // Calculate energies of all positions in chain
    Energy(E_L, alpha, N_steps, R1, R2);
    x_distribution(x_chain, N_steps, R1,R2);
    printf("Accept_count = %d \n", accept_count);
    
    // Save in csv:s
    save_matrix_to_csv( R1, N_steps, NDIM, filename_R1);
    save_matrix_to_csv( R2, N_steps, NDIM, filename_R2);

    bool open_with_write = true;
    double param_vector[] = {N_steps, alpha, d_displacement};
    save_vector_to_csv(param_vector, 3, filename_params, true);
    save_vector_to_csv(E_L, N_steps, filename_energy, open_with_write);
    save_vector_to_csv(x_chain, N_steps, filename_xdist, open_with_write);

    // Destroy and free arrays
    destroy_2D_array(R1,N_steps); destroy_2D_array(R2,N_steps);
    gsl_rng_free(r); free(E_L), free(x_chain);

    return 0;
}
