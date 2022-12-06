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
#define M_C 1000

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

void MCMC_burn_in(int N_steps, double alpha, double d_displacement, double **R1, double **R2)
{
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
            random_number = gsl_ran_flat(r, -0.5,0.5);
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
    }
    gsl_rng_free(r);
}

double MCMC(int N_steps, double alpha, double d_displacement, double **R1, double **R2)
{
    // Filenames for saving in csv
    char filename_R1[] = {"R1.csv"}, filename_R2[] = {"R2.csv"};
    char filename_energy[] = {"E_local.csv"}, filename_xdist[] = {"x_distribution.csv"}, filename_theta[] = {"theta.csv"};
    char filename_energy_derivative[] = {"E_local_derivative.csv"};
    char filename_results[] = {"MCMC_results.csv"}, filename_phi_k[] ={"phi_k.csv"}, filename_block_avg[] ={"block_avg_vec.csv"};
    bool open_with_write;
    
     // Initializing arrays
    int n_phi_rows = 2*M_C+10;
    int number_of_blocks = 3;
    double *E_local = malloc(sizeof(double) * N_steps);
    double *E_local_derivative = malloc(sizeof(double) * N_steps);
    double *Phi_k_vec = malloc(sizeof(double) *n_phi_rows);
    double *block_average_vec = malloc(sizeof(double) *number_of_blocks);
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
            random_number = gsl_ran_flat(r, -0.5, 0.5);
            R1_test[kx] = R1[ix][kx] + d_displacement * random_number;
            random_number = gsl_ran_flat(r, -0.5, 0.5);
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

        //double result_vec[] = {ix, E_PD_average};
        //if(ix == 0){ open_with_write = true; } else { open_with_write = false; }
        //save_vector_to_csv(result_vec, 2, filename_results, open_with_write);
        //printf("MCMC step: %d\n", ix);
    }

    // Calculate energies of all positions in chain
    Energy(E_local, alpha, N_steps, R1, R2); 
    
    double E_PD_average = partialEnergyDerivative(E_local_derivative, alpha, N_steps, R1, R2);

    double average_E_local = 0;
    for(int ix = 0; ix < N_steps - 1; ++ix)
    {
        average_E_local += E_local[ix]/N_steps;
    }

    theta_fun(theta_chain, N_steps, R1, R2);
    x_distribution(x_chain, N_steps, R1,R2);
    double statistical_inefficiency = correlation_function(Phi_k_vec, E_local, N_steps, M_C);
    printf("statistical inefficiency from correlation function= %f\n", statistical_inefficiency);
    statistical_inefficiency=0;
    statistical_inefficiency = block_average(block_average_vec,E_local, N_steps, number_of_blocks);

    //printf("Accept_count = %d \n", accept_count);
    printf("statistical inefficiency from block averaging= %f\n", statistical_inefficiency);
    
    // Save in csv:s
    save_matrix_to_csv(R1, N_steps, NDIM, filename_R1);
    save_matrix_to_csv(R2, N_steps, NDIM, filename_R2);

    open_with_write = true;
    save_vector_to_csv(E_local_derivative, N_steps, filename_energy_derivative, open_with_write);
    save_vector_to_csv(E_local, N_steps, filename_energy, open_with_write);
    save_vector_to_csv(x_chain, N_steps, filename_xdist, open_with_write);
    save_vector_to_csv(Phi_k_vec, n_phi_rows, filename_phi_k, open_with_write);

    // Destroy and free arrays
    free(E_local), free(E_local_derivative), free(x_chain), free(theta_chain), free(Phi_k_vec), free(block_average_vec);
    gsl_rng_free(r);

    return average_E_local;
}

int
run(
    int argc,
    char *argv[]
   )
{
    // MCMC Parameters
    int N_steps; int N_discarded_steps; double alpha, d_displacement; 
    // alpha Parameters
    int n_alpha_steps; double A, beta, E_average;

    bool is_task1 = false, is_task2 = false, is_task3 = true, is_task4 = false;

    if(is_task1)
    {
        N_steps = 1e5; N_discarded_steps = 0; alpha = 0.1, d_displacement = 0.1; 
        n_alpha_steps = 1; A = 0.; beta = 0.; 
    }
    if(is_task2)
    {
        N_steps = 1e5; N_discarded_steps = 1e4; alpha = 0.1, d_displacement = 0.1; 
        n_alpha_steps = 1; A = 0.; beta = 0.; 
    }
    if(is_task3)
    {
        N_steps = 1e5; N_discarded_steps = 1e4; alpha = 0.05, d_displacement = 0.1; 
        n_alpha_steps = 1; A = 0.; beta = 0.;
    }
    if(is_task4)
    {
        N_steps = 1e5; N_discarded_steps = 1e4; alpha = 0.1, d_displacement = 0.1; 
        n_alpha_steps = 50; A = 1.; beta = 1.; // beta from 0.5 to 1
    }

    double **R1 = create_2D_array(N_steps, NDIM), **R2 = create_2D_array(N_steps,NDIM), E_PD_average;
    double *E_local_derivative = malloc(sizeof(double) * N_steps);
    char filename_alpha_results[] = {"alpha_results.csv"}, filename_alpha_params[] = {"alpha_params.csv"};
    bool open_with_write;

    initialize_positions((double **) R1, (double **) R2, (double) d_displacement);
    
    if(is_task3 || is_task4)
    {
        MCMC_burn_in(N_discarded_steps, alpha, d_displacement, R1, R2);
    }

    E_average = 0;
    for(int ix = 1; ix < n_alpha_steps + 1; ix++)
    {
        double gamma = A*pow(ix, (double) - beta);

        E_average = MCMC(N_steps, alpha, d_displacement, R1, R2);

        E_PD_average = partialEnergyDerivative(E_local_derivative, alpha, N_steps, R1, R2);

        alpha -= gamma * E_PD_average;

        double alpha_result_vector[] = {ix, E_average, alpha, gamma, E_PD_average};
        if(ix == 1){ open_with_write = true; } else { open_with_write = false; }
        save_vector_to_csv(alpha_result_vector, 5, filename_alpha_results, open_with_write);
        printf("Iteration: %d\n", ix);
    }

    double alpha_param_vector[] = {n_alpha_steps, N_discarded_steps, alpha, A, beta, N_steps, d_displacement, is_task1, is_task2, is_task3, is_task4};
    save_vector_to_csv(alpha_param_vector, 11, filename_alpha_params, true);

    destroy_2D_array(R1, N_steps); destroy_2D_array(R2, N_steps);
    free(E_local_derivative);

    return 0;
}
