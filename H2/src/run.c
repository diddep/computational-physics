#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include <stdbool.h>

#include "tools.h"
#include "distribution.h"
#include "MCMC_chain_operations.h"
#include "statistical_ineff.h"

#define NDIM 3

// Function running the markov-chain
void initialize_positions(double **R1, double **R2, double d_displacement)
{
    double random_number = 0;
    gsl_rng * r;
    r = init_random_num_generator();

    printf("Initializing random positions:\n");
    for (int kx = 0; kx < NDIM; kx++)
    {
        // -0.5 to 0.5 because it's length 1 and symmetric around zero
        random_number = gsl_ran_flat(r, -0.5, 0.5);
        R1[0][kx] = d_displacement * random_number;

        random_number = gsl_ran_flat(r, -0.5, 0.5);
        R2[0][kx] = d_displacement * random_number;
        printf("R1[0][%d]: %f, R2[0][%d]: %f\n", kx, R1[0][kx], kx, R2[0][kx]);
    }
    gsl_rng_free(r);
}

void MCMC(int N_steps, double alpha, double d_displacement, double **R1, double **R2)
{
    // Filenames for saving in csv
    char filename_R1[] = {"R1.csv"}, filename_R2[] = {"R2.csv"};
    char filename_energy[] = {"E_local.csv"}, filename_xdist[] = {"x_distribution.csv"}, filename_theta[] = {"theta.csv"};
    char filename_results[] = {"MCMC_results.csv"};
    bool open_with_write;
    
     // Initializing arrays
    double *E_local = malloc(sizeof(double) * N_steps);
    double *theta_chain = malloc(sizeof(double) * N_steps);
    double *x_chain = malloc(sizeof(double) * N_steps);

    double R1_test[NDIM], R2_test[NDIM];
    
    gsl_rng * r;
    r = init_random_num_generator();
    double random_number = 0;

    int accept_count = 0;
    for(int ix = 0; ix < N_steps - 1; ++ix)
    {
        // get proposal positons
        for (int kx = 0; kx < NDIM; ++kx)
        {
            random_number = gsl_ran_flat(r,-0.5,0.5);
            R1_test[kx] = R1[ix][kx] + d_displacement * random_number;
            random_number = gsl_ran_flat(r,-0.5,0.5);
            R2_test[kx] = R2[ix][kx] + d_displacement * random_number;
        }
        
        // Probability for particle occupying new and old positions
        double prob_test = distribution(R1_test, R2_test, alpha);
        double prob_old = distribution(R1[ix], R2[ix], alpha);

        // If new prob > old, make step OR take exploration step
        if(prob_test > prob_old || prob_test / prob_old > gsl_ran_flat(r, 0.0, 1.0))
        {
            // If accepted, save new position in next row.
            for (int kx=0; kx < NDIM; ++kx)
            {
                R1[ix+1][kx] = R1_test[kx];
                R2[ix+1][kx] = R2_test[kx];
            }
            
            accept_count = accept_count + 1;
        } else {
            // If not accepted save old position in next row
            for (int kx=0; kx < NDIM; ++kx)
            {
                R1[ix+1][kx] = R1[ix][kx];
                R2[ix+1][kx] = R2[ix][kx];
            }
        }

        double theta_ix = theta_fun_vec(R1[ix], R2[ix]);
        double x_cos = cos(theta_ix);

        double result_vec[] = {ix, x_cos};
        if(ix == 0){ open_with_write = true; } else { open_with_write = false; }
        save_vector_to_csv(result_vec, 2, filename_theta, open_with_write);
    }

    // Calculate energies of all positions in chain
    Energy(E_local, alpha, N_steps, R1, R2);
    theta_fun(theta_chain, N_steps, R1, R2);
    x_distribution(x_chain, N_steps, R1,R2);
    printf("Accept_count = %d \n", accept_count);
    
    // Save in csv:s
    save_matrix_to_csv(R1, N_steps, NDIM, filename_R1);
    save_matrix_to_csv(R2, N_steps, NDIM, filename_R2);

    open_with_write = true;
    save_vector_to_csv(E_local, N_steps, filename_energy, open_with_write);
    save_vector_to_csv(x_chain, N_steps, filename_xdist, open_with_write);

    // Destroy and free arrays
    free(E_local), free(x_chain), free(theta_chain);
    gsl_rng_free(r);
}

int
run(
    int argc,
    char *argv[]
   )
{
    // MCMC Parameters
    // N_steps is number of steps in loop, alpha --- , d_displacement offsets particle positions
    int N_steps = 1e5; int N_discarded_steps = 1e4; double alpha = 0.1, d_displacement = 0.1; 

    double **R1 = create_2D_array(N_steps,NDIM), **R2 = create_2D_array(N_steps,NDIM);
    double E_PD_average;

    initialize_positions((double **) R1, (double **) R2, (double) d_displacement);

    MCMC(N_discarded_steps, alpha, d_displacement, R1, R2);

    E_PD_average = partialEnergyDerivative(alpha, N_discarded_steps, R1, R2);
    printf("E_PD discarded: %f\n", E_PD_average);

    MCMC(N_steps, alpha, d_displacement, R1, R2);

    E_PD_average = partialEnergyDerivative(alpha, N_steps, R1, R2);
    printf("E_PD rest: %f\n", E_PD_average);

    char filename_params[] = {"MCMC_params.csv"};
    double param_vector[] = {N_steps, alpha, d_displacement};
    save_vector_to_csv(param_vector, 3, filename_params, true);

    destroy_2D_array(R1, N_steps); destroy_2D_array(R2, N_steps);

    return 0;
}
