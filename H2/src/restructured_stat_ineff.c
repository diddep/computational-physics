
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "tools.h"

double phi_lag(double *E_local_chain, int N_steps, int Lag)
{
    double phi_k = 0, average_E_local = 0, average_squared_E_local=0, lagged_average=0;
    int lower_buffer =Lag, buffer_upper = N_steps-Lag;

    for(int step=0; step<N_steps; ++step)
    {
        double E_sample = E_local_chain[step];
        average_E_local += E_sample;
        average_squared_E_local += E_sample*E_sample;
    }

    average_E_local /=N_steps; average_squared_E_local /=N_steps;

    for(int step=0; step< buffer_upper; ++step)
    {
        lagged_average += E_local_chain[step]*E_local_chain[step+Lag];
    }
    lagged_average/=(N_steps-Lag);

    phi_k = (lagged_average -average_E_local*average_E_local)/(average_squared_E_local-average_E_local*average_E_local);

    return phi_k;
}

double statistical_ineff_from_BLAV(double *E_local_vec, int N_steps, int Block_size)
{
    int number_of_blocks = N_steps/Block_size;
    double average_i = 0, variance_block=0, variance_E_local=0, statistical_inneficiency;
    double average_tot =0;
    // array to save averages for each block
    double *block_vec = malloc(sizeof(double)*number_of_blocks);

    for(int block=0; block<number_of_blocks; ++block)
    {
        //calculate average inside a block
        //for(int step =block*Block_size; step<(block+1)*Block_size; ++step)
        for(int step =0; step< Block_size; ++step)
        {
            block_vec[block] += E_local_vec[Block_size*(block) +step]; 
            //block_vec[block] += E_local_vec[step];
        }
        block_vec[block] /=Block_size;
    }
    variance_block = variance(block_vec,number_of_blocks);
    variance_E_local = variance(E_local_vec,N_steps);

    statistical_inneficiency = (variance_block*Block_size)/variance_E_local;
    free(block_vec);
    return statistical_inneficiency;
}


// double restructured_block_average(double *E_local_vec, int N_steps, int number_of_blocks)
// {
//     int block_size = N_steps/number_of_blocks;
//     double average_i = 0, variance_block=0, variance_E_local=0, statistical_inneficiency, error_bar=0;
//     double average_tot =0;
    
//     for(int step =0; step<N_steps; ++step)
//     {
//         average_tot += E_local_vec[step];
//     }
//     average_tot /=N_steps;

    

// }

