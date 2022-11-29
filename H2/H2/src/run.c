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

int
run(
    int argc,
    char *argv[]
   )
{
    int N = 10000;
    double alpha = 0.1;
    double d = 1.0;
    int acceptance =0;
    double** R1 = create_2D_array(N,3);
    double** R2 = create_2D_array(N,3);
    double R1_test[3];
    double R2_test[3];
    double random_number1 =0;
    double random_number2 =0;
    char filename_R1[] = {"R1.csv"};
    char filename_R2[] = {"R2.csv"};

    printf("first \n");

    gsl_rng * r;
    r = init_random_num_generator();


    for (int kx=0; kx<3; ++kx)
    {
        random_number1 = gsl_ran_flat(r,-0.5,0.5);
        random_number2 = gsl_ran_flat(r,-0.5,0.5);

        R1[0][kx] = d * random_number1; //(((double) rand()% 1000)- 500.0)/1000.0;
        R2[0][kx] = d * random_number2; //(((double) rand()% 1000)- 500.0)/1000.0;
        printf("initvalue: %f,   %f\n", R1[0][kx], R2[0][kx]);
    }

    for(int i=0; i<N-1; ++i)
    {


        for (int kx=0; kx<3; ++kx)
        {
            random_number1 = gsl_ran_flat(r,-0.5,0.5);
            random_number2 = gsl_ran_flat(r,-0.5,0.5);

            R1_test[kx] = R1[i][kx] + d * random_number1;//(((double) rand()% 1000)- 500.0)/1000.0;
            R2_test[kx] = R2[i][kx] + d * random_number2;//((double) rand(1000)- 500.0)/1000.0;
            //printf("testvals: %f,   %f\n", R1_test[kx], R2_test[kx]);
        }

        double prob_test = distribution(R1_test, R2_test, alpha);
        double prob_old = distribution(R1[i], R2[i], alpha);
        //printf("probabilities: %f,   %f\n", prob_test, prob_old);

        if(prob_test>prob_old|| prob_test/prob_old > gsl_ran_flat(r,0.0,1.0))
        {
            //printf("accept\n");
            acceptance = acceptance +1;

            //printf("%d \n", acceptance);
            for (int kx=0; kx<3; ++kx)
            {
                //printf("hej");
                R1[i+1][kx] = R1_test[kx];
                R2[i+1][kx] = R2_test[kx];
                //printf("new value: %f,   %f\n", R1[i+1][kx], R2[i+1][kx]);
            }
        }
        else{
            for (int kx=0; kx<3; ++kx)
            {
                //printf("samma");
                R1[i+1][kx] = R1[i][kx];
                R2[i+1][kx] = R2[i][kx];
                //printf("new value: %f,   %f\n", R1[i+1][kx], R2[i+1][kx]);
            }
        }
    }

    save_matrix_to_csv(R1, N, 3, filename_R1);
    save_matrix_to_csv(R2, N, 3, filename_R2);
    printf("acceptance =%d \n", acceptance);

    destroy_2D_array(R1,N);
    destroy_2D_array(R2,N);
    gsl_rng_free(r);

    return 0;
}
