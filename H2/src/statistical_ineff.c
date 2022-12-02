
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


#include "tools.h"

// len(Phi_k_vec) = 2*M_c +1?
double correlation_function(double *Phi_k_vec ,double *E_local_vec, int N_steps, int M_c)
{
    double phi_k = 0, average_E_local = 0, average_squared_E_local=0, statistial_inefficiency=0;
    int lower_buffer =M_c, buffer_upper = N_steps-M_c;

    //Calculating the average over the chain, should maybe change to only bufferzone
    for(int step=0; step < N_steps;++step)
    {
        average_E_local += E_local_vec[step] /N_steps;
        
        average_squared_E_local += pow(E_local_vec[step],2) /N_steps;
    }
    //printf("avg energy = %f\n" , (average_E_local));

    //Looping over k for phi_k between -M_c<k<M_c
    for(int kx =-M_c; kx<M_c; ++kx)
    {
        //calculating averages
        for(int step=lower_buffer; step<buffer_upper; ++step)
        {
            phi_k = E_local_vec[step]*E_local_vec[step + kx] / (N_steps-2*M_c);
            //phi_k = (phi_k -pow(average_E_local,2)/(average_squared_E_local- pow(average_E_local,2)));
            //statistial_inefficiency += phi_k;
            //printf("phi_k= %f\n=",phi_k);
        }

        phi_k = (phi_k -pow(average_E_local,2)/(average_squared_E_local- pow(average_E_local,2)));
        statistial_inefficiency += phi_k;
        //tricky indexing, should be array with 2M_c elements
        Phi_k_vec[M_c+kx] = phi_k;
    }

    return statistial_inefficiency;
}


double block_average(double *block_average_vec, double *E_local_vec, int N_steps, int number_of_blocks)
{
    int block_size = N_steps/number_of_blocks;
    double average_i = 0, variance_block=0, variance_E_local=0, statistical_inneficiency, error_bar=0; 

    for(int block=0; block< number_of_blocks; ++block)
    {
        for(int step=0; step<block_size; ++step)
        {
            //should maybe be [block*(block_size-1) +step]
            average_i += E_local_vec[block*(block_size) +step] /(block_size);
            //printf("avava =%f\n", average_i);
        }
        block_average_vec[block] = average_i;
    }

    variance_block = variance(block_average_vec, number_of_blocks);
    printf("block variance =%f\n",variance_block);
    variance_E_local = variance(E_local_vec, N_steps);
    printf("total variance =%f\n",variance_E_local);
    statistical_inneficiency = (double) number_of_blocks* variance_block/variance_E_local;
    printf("stateff = %f\n", statistical_inneficiency);
    error_bar= sqrt(statistical_inneficiency*variance_E_local/ N_steps);
    printf("error calculated from block avg =+- %f\n", error_bar);

    return statistical_inneficiency;
}