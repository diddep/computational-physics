#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "tools.h"
#include "integrals.h"


void get_points(double *array_of_points, int number_of_points, double lower_bound, double upper_bound, gsl_rng *r)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[ix] = gsl_ran_flat(r, lower_bound, upper_bound);
    }
}

void pdf_transformed_points_task2(double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[0][ix] = acos(-2.*array_of_points[0][ix] +1.) / M_PI;
    }
}

// TODO: Update transform to relevant function
void pdf_transformed_points_task3(double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < 3; ix++)
    {
        for(int jx = 0; jx < number_of_points; jx++)
        {
            // TODO: eta_1 and eta_2 should be two different numbers, see page 10 in MC lecture notes
            double eta_1 = array_of_points[ix][jx], eta_2 = array_of_points[ix][jx];
            //array_of_points[ix][jx] = sqrt(-2. * eta_1) * cos(2. * M_PI * eta_2);
            array_of_points[ix][jx] *= 1;
        }
    }
}

void function_value_task1(double *function_value, double **array_of_points, int number_of_points)
{
    double x_coord[number_of_points];
    for(int ix = 0; ix < number_of_points; ix++)
    {
        double x_coord = array_of_points[0][ix];
        function_value[ix] = x_coord * (1 - x_coord);
    }
}

void function_value_task3_unweighted(double *function_value, double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        double x_coord = array_of_points[0][ix], y_coord = array_of_points[1][ix], z_coord = array_of_points[2][ix];

        double factor_1 = pow(M_PI,(double) -3./2.);
        double factor_2 = pow(x_coord, 2.) + pow(x_coord * y_coord, 2.)  + pow(x_coord * y_coord * z_coord, 2.);
        double factor_3 = exp(-(pow(x_coord, 2.) + pow(y_coord, 2.) +pow(z_coord, 2.)));
        function_value[ix] = factor_1 * factor_2 * factor_3;
    }
}

void function_value_task3_weighted(double *function_value, double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        double x_coord = array_of_points[0][ix], y_coord = array_of_points[1][ix], z_coord = array_of_points[2][ix];

        double factor_1 = pow(M_PI,(double) -3./2.);
        double factor_2 = pow(x_coord, 2.) + pow(x_coord * y_coord, 2.)  + pow(x_coord * y_coord * z_coord, 2.);
        double factor_3 = exp(-(pow(x_coord, 2.) + pow(y_coord, 2.) +pow(z_coord, 2.)));
        function_value[ix] = factor_2;
    }
}

void weight_function_task2(double *weight_vector, double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        double x_coord = array_of_points[0][ix];
        weight_vector[ix] = M_PI / 2. * sin(M_PI * x_coord);
    }
}

void weight_function_task3(double *weight_vector,double **array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        double x_coord = array_of_points[0][ix], y_coord = array_of_points[1][ix], z_coord = array_of_points[2][ix];
        //FIXME: Somewhere here is a bug that disrupts the weight factor causing nans
        double factor_1 = pow(M_PI, (double) -3./2.);
        double term_1 = pow(x_coord, 2.), term_2 = pow(y_coord, 2.), term_3 = pow(z_coord, 2.);

        double factor_3 = exp(- 1. * (term_1 + term_2 + term_3));
        printf("factor_1 = %f\n", factor_1);
        printf("factor_3 = %f\n", factor_3);
        weight_vector[ix] = factor_1 * factor_3;
        printf("weigt[%d] = %f\n", ix, weight_vector[ix]);
    }
}

void evaluate_integral(double **array_of_points, int number_of_points, double *integral_results, bool is_integral_1, bool is_weighted)
{
    double integral_value = 0, integral_std = 0;
    double variance = 0; int ndim;
    double *function_value = calloc(number_of_points, sizeof(double));
    double *weight_function = calloc(number_of_points, sizeof(double));

    if(is_integral_1)
    {
        ndim = 1;
        function_value_task1(function_value, array_of_points, number_of_points);
        if(is_weighted)
        {
            weight_function_task2(weight_function, array_of_points, number_of_points);
            for(int ix = 0; ix < number_of_points; ix++)
            {
                function_value[ix] /= weight_function[ix];
            }
        }
    } else {
        ndim = 3;
        function_value_task3_unweighted(function_value, array_of_points, number_of_points);
        
        if(is_weighted)
        {
            //function_value_task3_weighted(function_value, array_of_points, number_of_points);
        
            weight_function_task3(weight_function, array_of_points, number_of_points);
            for(int ix = 0; ix < number_of_points; ix++)
            {
                function_value[ix] /= weight_function[ix];
            }
        }
    }
    
    for(int ix = 0; ix < number_of_points; ix++)
    {
        integral_value += function_value[ix] / number_of_points;
        variance += pow(function_value[ix], 2.) / number_of_points \
                    - pow(function_value[ix] / number_of_points, 2.); 

    }

    integral_std = sqrt(variance / number_of_points);

    integral_results[0] = integral_value; 
    integral_results[1] = integral_std; 
}

void MC_integration(bool is_integral_1, bool is_weighted, int number_of_calculations, int number_of_results)
{
    double **calculation_results = create_2D_array(number_of_calculations, number_of_results);
    // TODO: When switching from task 1 to task3, change file-filder in name below
    char filename_calculation_results[] = {"../task3/calculation_results.csv"};
    char filename_sampled_points[] = {"../task3/sampled_points.csv"};
    gsl_rng *r;
    r = init_random_num_generator();

    for(int ix = 0; ix < number_of_calculations; ix++)
    {   
        double **array_of_points; int number_of_points, ndim;
        if(is_integral_1)
        {
            ndim = 1;
            number_of_points = pow((int) 10., (int) ix + 1);
            
            array_of_points = create_2D_array(ndim, number_of_points);
            double lower_bound = 0, upper_bound = 1;
            get_points(array_of_points[0], number_of_points, lower_bound, upper_bound, r);

            if(is_weighted)
            {
                pdf_transformed_points_task2(array_of_points, number_of_points);
            }
        } else {
            ndim = 3;
            number_of_points = 1e4;
            array_of_points = create_2D_array(number_of_points, ndim);

            // TODO: Update lower and upper bounds
            double lower_bound = -1e1, upper_bound = 1e1;
            for(int jx = 0; jx < ndim; jx++)
            {
                get_points(array_of_points[jx], number_of_points, lower_bound, upper_bound, r);
            } 
            if(is_weighted)
                {
                    pdf_transformed_points_task3(array_of_points, number_of_points);
                }
        }
       
        //double integral_value = 0, integral_std = 0;
        //integral_results[0] = integral_value;  integral_results[1] = integral_std;
        double integral_results[2];
        evaluate_integral(array_of_points, number_of_points, integral_results, is_integral_1, is_weighted);

        for(int jx = 0; jx < number_of_results; jx++)
        {
            calculation_results[ix][jx] = integral_results[jx];
        }
        save_matrix_to_csv(array_of_points, ndim, number_of_points, filename_sampled_points);
    }
    save_matrix_to_csv(calculation_results, number_of_calculations, number_of_results, filename_calculation_results);
    destroy_2D_array(calculation_results, number_of_calculations);
    gsl_rng_free(r);
}

double probability_task_3(double *coordinate_vec)
{
    double probability =0;
    double x= coordinate_vec[0];
    double y = coordinate_vec[1];
    double z = coordinate_vec[2];
    double arg = -(x*x + y*y + z*z);

    probability = exp(arg);
    return probability;
}

double integrand_task_3(double** mcmc_chain, double *integrand_chain ,int N_samples)
{
    double integrand = 0, integrand_average=0;
    double x=0, y=0,z=0; 
    char filename_integrand[] = {"integrand_chain_task3.csv"};
    double denominator = pow(M_PI, (double) 3/2);

    
    for(int step=0; step< N_samples; ++step)
    {   
        x = mcmc_chain[step][0], y = mcmc_chain[step][1], z =mcmc_chain[step][2];
        integrand = (x*x +x*x*y*y+x*x*y*y*z*z)*exp(-(x*x + y*y + z*z));
        integrand /= denominator; 
        integrand_chain[step]= integrand;
        integrand_average += integrand/N_samples;
    }
    bool open_with_write = true;
    save_vector_to_csv(integrand_chain, N_samples, filename_integrand, open_with_write);

    return integrand_average;
}

void markov_chain_monte_carlo(double **mcmc_chain, int N_samples, double step_size)
{
    gsl_rng *r;
    r = init_random_num_generator();
    int dimension=3, accepts=0;
    double test_position[dimension], old_position[dimension];
    double random_number = 0, probability_new=0, probability_old=0;
    char filename_mcmc[] = {"mcmc_chain_task3.csv"};
    for(int coordinate=0; coordinate<dimension; ++coordinate)
    {
        random_number = gsl_ran_flat(r, -0.5,0.5);
        mcmc_chain[0][coordinate] = step_size*random_number;
    }

    for(int step=0; step< N_samples-1; ++step)
    {
        for(int coordinate=0; coordinate< dimension; ++coordinate)
        {
            old_position[coordinate] = mcmc_chain[step][coordinate];
        }
        for(int coordinate=0; coordinate< dimension; ++coordinate)
        {
            random_number = gsl_ran_flat(r, -0.5,0.5);
            test_position[coordinate] = mcmc_chain[step][coordinate] + step_size*random_number;
        }
        probability_new =probability_task_3(test_position);
        probability_old = probability_task_3(old_position);

        if(probability_new> probability_old||probability_new/probability_old>gsl_ran_flat(r, 0.0, 1.0))
        {
            accepts+=1;

            for(int coordinate=0; coordinate< dimension; ++coordinate)
            {
                mcmc_chain[step+1][coordinate] = test_position[coordinate];
            }
        }
        else
        {
            for(int coordinate=0; coordinate< dimension; ++coordinate)
            {
                mcmc_chain[step+1][coordinate] = old_position[coordinate];
            }
        }
    }
    save_matrix_to_csv(mcmc_chain, N_samples, dimension, filename_mcmc);

    gsl_rng_free(r);
}

void task1()
{
    bool is_integral_1 = true, is_weighted = true;
    int number_of_calculations = 4, number_of_results = 2;
    MC_integration(is_integral_1, is_weighted, number_of_calculations ,number_of_results);

}

void task3()
{
    bool is_integral_1 = false, is_weighted = true;
    int number_of_calculations = 1, number_of_results = 2;
    MC_integration(is_integral_1, is_weighted, number_of_calculations, number_of_results);
}

void task3_correct()
{
    int N_samples = 100000, dimension=3; 
    double step_size = 0.5;
    double **mcmc_chain = create_2D_array(N_samples, dimension);
    double *integrand_chain= malloc(sizeof(double)*N_samples);
    markov_chain_monte_carlo(mcmc_chain, N_samples,step_size);
    double integral = integrand_task_3(mcmc_chain, integrand_chain, N_samples);

    destroy_2D_array(mcmc_chain, N_samples);
    free(integrand_chain);
    printf("integral value = %f\n", integral);
}

int
run(int argc, char *argv[])
{
    //task1();
    //task3();
    task3_correct();
}