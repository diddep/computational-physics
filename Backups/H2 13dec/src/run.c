#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>

#include <stdbool.h>

#include "tools.h"
#include "distribution.h"
#include "MCMC_chain_operations.h"
#include "statistical_ineff.h"
#include "restructured_stat_ineff.h"

#define NDIM 3
#define M_C 1000

void initialize_positions(double **R1, double **R2, double d_displacement);
void MCMC_burn_in(int N_steps, double alpha, double d_displacement, double **R1, double **R2);
double MCMC(int N_steps, double alpha, double d_displacement, double **R1, double **R2, bool is_save);


int
run(
    int argc,
    char *argv[]
   )
{
    bool is_task1 = false, is_task2 = false, is_task3 = false, is_task4 =  false;
    int task_num = 1;
    if (argc < 2) {
        task_num = 1;
        //printf("Usage: %s TASK\n", argv[0]);
        //exit(1);
    } else {
        task_num = atol(argv[1]);

        if (task_num > 5) {
            printf("Task should be between 1 and 5\n");
            exit(1);
        }
    }


    // MCMC Parameters
    int N_steps; int N_discarded_steps; double alpha, d_displacement; 
    // alpha Parameters
    int N_alpha_steps; double A, beta, E_average;
    bool is_save = true;

    if(task_num == 1)
    {
        N_steps = 1e6; N_discarded_steps = 0; alpha = 0.1, d_displacement = 1.24; 
        N_alpha_steps = 1; A = 0.; beta = 0.; 
        is_task1 = true;
}
    if(task_num == 2)
    {
        N_steps = 1e6; N_discarded_steps = 1e4; alpha = 0.1, d_displacement = 1.24; 
        N_alpha_steps = 1; A = 0.; beta = 0.; 
        is_task2 = true;
    }
    if(task_num == 3)
    {
        N_steps = 1e6; N_discarded_steps = 1e4; alpha = 0.05, d_displacement = 1.24; 
        N_alpha_steps = 1; A = 0.; beta = 0.; 
        is_task3 = true;
    }
    if(task_num == 4)
    {
        N_steps = 1e6; N_discarded_steps = 1e4; alpha = 0.1, d_displacement = 1.24; 
        N_alpha_steps = 50; A = 1.; beta = 1.; is_save = false; // beta from 0.5 to 1
        is_task4 = true;
    }

    double **R1 = create_2D_array(N_steps, NDIM), **R2 = create_2D_array(N_steps, NDIM);
    double E_PD_average;
    double *E_local_derivative = malloc(sizeof(double) * N_steps);

    char filename_alpha_results[200], filename_params[200];
    char cwd[200], buf[200];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        printf("Current working directory %s\n", cwd);
    } else {
        perror("getcwd() error");
        return 1;
    }
    
    int written = snprintf(buf, 200, "%s", cwd);
    snprintf(buf + written, 200 - written, "/csv/alpha_results.cvs");
    strcpy(filename_alpha_results, buf);
    snprintf(buf + written, 200 - written, "/csv/params.csv");
    strcpy(filename_params, buf);
    


    bool open_with_write;

    initialize_positions((double **) R1, (double **) R2, (double) d_displacement);
    
    if(is_task3 || is_task4)
    {
        MCMC_burn_in(N_discarded_steps, alpha, d_displacement, R1, R2);
    }

    E_average = 0;
    for(int ix = 1; ix < N_alpha_steps + 1; ix++)
    {
        if (ix == N_alpha_steps) {is_save = true; }
        E_average = MCMC(N_steps, alpha, d_displacement, R1, R2, is_save);

        E_PD_average = partialEnergyDerivative(E_local_derivative, alpha, N_steps, R1, R2);

        double gamma = 0;
        if(is_task3 || is_task4)
        {
            gamma = A*pow(ix, (double) - beta);
        }

        alpha -= gamma * E_PD_average;

        // TODO: Eventuell speed up att spara dessa i en matris en gång istället för en vector n_alpha_steps gånger
        double alpha_result_vector[] = {ix, E_average, alpha, gamma, E_PD_average};
        if(ix == 1){ open_with_write = true; } else { open_with_write = false; }
        save_vector_to_csv(alpha_result_vector, 5, filename_alpha_results, open_with_write);
        printf("Iteration: %d\n", ix);
    }

    double param_vector[] = {N_alpha_steps, N_discarded_steps, alpha, A, beta, N_steps, d_displacement, is_task1, is_task2, is_task3, is_task4};
    save_vector_to_csv(param_vector, 11, filename_params, true);
    destroy_2D_array(R1, N_steps); destroy_2D_array(R2, N_steps);
    free(E_local_derivative);

    return 0;
}


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
        if(prob_test / prob_old > gsl_ran_flat(r, 0.0, 1.0))
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

double MCMC(int N_steps, double alpha, double d_displacement, double **R1, double **R2, bool is_save)
{
    char filename_R1[200], filename_R2[200], filename_energy[200], filename_xdist[200], filename_theta[200], \
     filename_energy_derivative[200], filename_results[200], filename_phi_k[200], filename_block_avg[200];
    
    char cwd[200], buf[200];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        perror("getcwd() error");
        return 1;
    }
    
    int written = snprintf(buf, 200, "%s", cwd);
    snprintf(buf + written, 200 - written, "/csv/R1.cvs");
    strcpy(filename_R1, buf);
    snprintf(buf + written, 200 - written, "/csv/R2.csv");
    strcpy(filename_R2, buf);
    snprintf(buf + written, 200 - written, "/csv/E_local.cvs");
    strcpy(filename_energy, buf);
    snprintf(buf + written, 200 - written, "/csv/x_distribution.csv");
    strcpy(filename_xdist, buf);
    snprintf(buf + written, 200 - written, "/csv/theta.cvs");
    strcpy(filename_theta, buf);
    snprintf(buf + written, 200 - written, "/csv/E_local_derivative.cvs");
    strcpy(filename_energy_derivative, buf);
    snprintf(buf + written, 200 - written, "/csv/filename_results.csv");
    strcpy(filename_results, buf);
    snprintf(buf + written, 200 - written, "/csv/phi_k.cvs");
    strcpy(filename_phi_k, buf);
    snprintf(buf + written, 200 - written, "/csv/block_avg_vec.csv");
    strcpy(filename_block_avg, buf);

    bool open_with_write;
    
     // Initializing arrays
    int n_phi_rows = 2*M_C+10;
    int number_of_blocks = 100;
    double *E_local = malloc(sizeof(double) * N_steps);
    double *E_local_derivative = malloc(sizeof(double) * N_steps);
    double *Phi_k_vec = malloc(sizeof(double) *n_phi_rows);
    double *block_average_vec = malloc(sizeof(double) *number_of_blocks);
    double *theta_chain = malloc(sizeof(double) * N_steps);
    double *x_chain = malloc(sizeof(double) * N_steps);

    // testing corr func

    int max_lag = 2*1e2;
    double phi_k=0;
    double *phi_k_vec = malloc(sizeof(double)*max_lag);
    

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
        if(prob_test / prob_old > gsl_ran_flat(r, 0.0, 1.0))
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
        theta_chain[ix] = theta_ix;

    }

    // Calculate energies of all positions in chain
    Energy(E_local, alpha, N_steps, R1, R2); 
    
    for (int lag=0; lag< max_lag; ++lag)
    {
        double phi_inst =phi_lag(E_local, N_steps, lag); 
        phi_k_vec[lag] = phi_inst;
        //printf("phi_k=%f\n", phi_inst);
    }


    double E_PD_average = partialEnergyDerivative(E_local_derivative, alpha, N_steps, R1, R2);

    double average_E_local = 0;
    for(int ix = 0; ix < N_steps - 1; ++ix)
    {
        average_E_local += E_local[ix]/N_steps;
    }

    if(is_save)
    {
        printf("Accept ratio = %f\n", (double) accept_count/N_steps);

        //theta_fun(theta_chain, N_steps, R1, R2);
        x_distribution(x_chain, N_steps, R1,R2);
        double statistical_inefficiency = correlation_function(Phi_k_vec, E_local, N_steps, M_C);
        printf("Statistical inefficiency from correlation function= %f\n", statistical_inefficiency);
        statistical_inefficiency=0;
        statistical_inefficiency = block_average(block_average_vec,E_local, N_steps, number_of_blocks);

        //printf("Accept_count = %d \n", accept_count);
        printf("Statistical inefficiency from block averaging= %f\n", statistical_inefficiency);

        // Save in csv:s
        save_matrix_to_csv(R1, N_steps, NDIM, filename_R1);
        save_matrix_to_csv(R2, N_steps, NDIM, filename_R2);

        double result_vec[] = {N_steps, accept_count};
        open_with_write = true;
        save_vector_to_csv(result_vec, 2, filename_results, open_with_write);
        save_transposedvector_to_csv(E_local_derivative, N_steps, filename_energy_derivative, open_with_write);
        save_transposedvector_to_csv(E_local, N_steps, filename_energy, open_with_write);
        save_transposedvector_to_csv(x_chain, N_steps, filename_xdist, open_with_write);
        save_transposedvector_to_csv(theta_chain, N_steps, filename_theta, open_with_write);
        save_transposedvector_to_csv(phi_k_vec, max_lag, filename_phi_k, open_with_write);
    }
    // Destroy and free arrays
    free(E_local), free(E_local_derivative), free(x_chain), free(theta_chain), free(Phi_k_vec), free(block_average_vec);
    gsl_rng_free(r);

    return average_E_local;
}
