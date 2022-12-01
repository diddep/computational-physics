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

// Function running the markov-chain
int
run(
    int argc,
    char *argv[]
   )
{
    // N number of steps, acceptan
    int N = 1e5, accept_count = 0;
    double alpha = 0.1, d = 0.1; 
    double **R1 = create_2D_array(N,3), **R2 = create_2D_array(N,3);
    double *E_L = malloc(sizeof(double) * N);
    double *x_distribution = malloc(sizeof(double) * N);

    double R1_test[3], R2_test[3];
    double random_number1 = 0, random_number2 = 0;
    char filename_R1[] = {"R1.csv"}, filename_R2[] = {"R2.csv"}, filename_energy[] = {"E_L.csv"}, filename_xdist[] = {"x_distribution.csv"};

    gsl_rng * r;
    r = init_random_num_generator();

    // Initializing random starting values for all particles
    printf("Initializing random positions:\n");
    for (int kx = 0; kx < 3; kx++)
    {
        // -0.5 to 0.5 because it's length 1 and symmetric around zero
        random_number1 = gsl_ran_flat(r, -0.5, 0.5);
        random_number2 = gsl_ran_flat(r, -0.5, 0.5);

        R1[0][kx] = d * random_number1;
        R2[0][kx] = d * random_number2;
        printf("R1[0][%d]: %f, R2[0][%d]: %f\n", kx, R1[0][kx], kx, R2[0][kx]);
    }

    for(int i = 0; i < N - 1; ++i)
    {
        // Get proposal position, d step size
        for (int kx = 0; kx < 3; ++kx)
        {
            random_number1 = gsl_ran_flat(r,-0.5,0.5);
            random_number2 = gsl_ran_flat(r,-0.5,0.5);

            R1_test[kx] = R1[i][kx] + d * random_number1;
            R2_test[kx] = R2[i][kx] + d * random_number2;
        }

        // Probability for particle occupying new and old positions
        double prob_test = distribution(R1_test, R2_test, alpha);
        double prob_old = distribution(R1[i], R2[i], alpha);
        //printf("prob old= %f\n", prob_old);
        //printf("prob test= %f\n", prob_test);
        //printf("prob kvot= %f\n", prob_test/prob_old);


        // If new prob > old, make step OR take exploration step
        if(prob_test > prob_old || prob_test / prob_old > gsl_ran_flat(r,0.0,1.0))
        {
            // Counter for how many accepts
            accept_count = accept_count + 1;

            // If accepted, save new position in next row.
            for (int kx=0; kx<3; ++kx)
            {
                R1[i+1][kx] = R1_test[kx];
                R2[i+1][kx] = R2_test[kx];
                //printf("new value: %f,   %f\n", R1[i+1][kx], R2[i+1][kx]);
            }
        }
        // If not accepted save old position in next row
        else{
            for (int kx=0; kx<3; ++kx)
            {
                R1[i+1][kx] = R1[i][kx];
                R2[i+1][kx] = R2[i][kx];
            }
        }
    }

    // Calculate energies of all positions in chain
    Energy(E_L, alpha, N, R1, R2);
    x_distribution(x_distribution, N, R1,R2);
    printf("Accept_count = %d \n", accept_count);
    
    // Save in csv:s
    save_matrix_to_csv( R1, N, 3, filename_R1);
    save_matrix_to_csv( R2, N, 3, filename_R2);

    bool open_with_write = true;
    save_vector_to_csv(E_L, N, filename_energy, open_with_write);
    save_vector_to_csv(x_distribution, N, filename_xdist, open_with_write);

    // Destroy and free arrays
    destroy_2D_array(R1,N); destroy_2D_array(R2,N);
    gsl_rng_free(r); free(E_L), free(x_distribution);

    return 0;
}
