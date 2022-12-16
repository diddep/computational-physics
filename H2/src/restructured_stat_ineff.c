
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "tools.h"


/*Function that takes in MCMC-chain for energy and calculates correlation function for a lagtime
    args:
        double *E_local_chain = 1D-array with length N_steps, is the mcmc-chain of the energy
        int N_steps = Number for steps in MCMC-chain
        int Lag = lag time for calculating correlation function.

    returns:
        phi_k the evaluated correlation function
*/
double phi_lag(double *E_local_chain, int N_steps, int Lag)
{
    //Initializing variables used in calculations
    double phi_k = 0, average_E_local = 0, average_squared_E_local=0, lagged_average=0;
    int lower_buffer =Lag, buffer_upper = N_steps-Lag;

    //calculating average energy in chain
    for(int step=0; step<N_steps; ++step)
    {
        double E_sample = E_local_chain[step];
        average_E_local += E_sample;
        average_squared_E_local += E_sample*E_sample;
    }

    //normalizing average
    average_E_local /=N_steps; average_squared_E_local /=N_steps;

    //calculating lagged average for a specified lag time
    for(int step=0; step< buffer_upper; ++step)
    {
        lagged_average += E_local_chain[step]*E_local_chain[step+Lag];
    }
    //normalizing lagged average
    lagged_average/=(N_steps-abs(Lag));

    //calculates correlation function phi_k
    phi_k = (lagged_average -average_E_local*average_E_local)/(average_squared_E_local-average_E_local*average_E_local);

    return phi_k;
}

/* Function that calculates statistical inneficieny from block averaging for a given block size

    args:
        double *E_local_chain = 1D-array with length N_steps, is the mcmc-chain of the energy
        int N_steps = Number for steps in MCMC-chain
        int Block_size = block size for calculating statistical inneficiency.

    returns:
        statistical_inneficiency= statistical inneficiency for a given block size
*/
double statistical_ineff_from_BLAV(double *E_local_vec, int N_steps, int Block_size)
{
    //initializing variables used in calculations
    int number_of_blocks = N_steps/Block_size;
    double average_i = 0, variance_block=0, variance_E_local=0, statistical_inneficiency;
    double average_tot =0;
    // array to save averages for each block
    double *block_vec = malloc(sizeof(double)*number_of_blocks);


    //looping through the blocks
    for(int block=0; block<number_of_blocks; ++block)
    {
        //calculate average inside a block
        for(int step =0; step< Block_size; ++step)
        {
            block_vec[block] += E_local_vec[Block_size*(block) +step]; 
            
        }
        //normalizes the average inside block
        block_vec[block] /=Block_size;
    }
    //calculates variance of block averages aswell as total variance of mcmc-chain
    variance_block = variance(block_vec,number_of_blocks);
    variance_E_local = variance(E_local_vec,N_steps);

    //calculates statistical inneficiency
    statistical_inneficiency = (variance_block*Block_size)/variance_E_local;
    free(block_vec);
    return statistical_inneficiency;
}
