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
    int N = 1000;
    double alpha = 0.1;
    double d = 1.0;
    int acceptance =0;
    double** R1 = create_2D_array(N,3);
    double** R2 = create_2D_array(N,3);
    double R1_test[3];
    double R2_test[3];
    double random_number1 =0;
    double random_number2 =0;

    gsl_rng * r;
    r = init_random_num_generator();


    for (int k; k<3; ++k)
    {
        random_number1 = gsl_ran_flat(r,-0.5,0.5);
        random_number2 = gsl_ran_flat(r,-0.5,0.5);

        R1[0][k] = d * random_number1; //(((double) rand()% 1000)- 500.0)/1000.0;
        R2[0][k] = d * random_number2; //(((double) rand()% 1000)- 500.0)/1000.0;
    }

    for(int i=1; i<N; ++i)
    {

        for (int k; k<3; ++k)
        {
            random_number1 = gsl_ran_flat(r,-0.5,0.5);
            random_number2 = gsl_ran_flat(r,-0.5,0.5);

            R1_test[k] = R1[i][k] + d * random_number1;//(((double) rand()% 1000)- 500.0)/1000.0;
            R2_test[k] = R2[i][k] + d * random_number2;//((double) rand(1000)- 500.0)/1000.0;
        }

        double prob_test = distribution(R1_test, R2_test, alpha);
        double prob_old = distribution(R1[0], R2[0], alpha);

        if(prob_test>prob_old|| prob_test/prob_old > gsl_ran_flat(r,0.0,1.0))
        {
            acceptance +=1;
            for (int k; k<3; ++k)
            {
                printf("hej");
                R1[i+1][k] = R1_test[k];
                R2[i+1][k] = R2_test[k];
            }
        }
        else{
            for (int k; k<3; ++k)
            {
                R1[i+1][k] = R1[i][k];
                R2[i+1][k] = R2[i][k];
            }

        }
    }
    return 0;
}
