
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

