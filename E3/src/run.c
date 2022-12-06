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

void MCMC_integration(bool is_integral_1, bool is_weighted, int number_of_calculations, int number_of_results)
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

void task1()
{
    bool is_integral_1 = true, is_weighted = true;
    int number_of_calculations = 4, number_of_results = 2;
    MCMC_integration(is_integral_1, is_weighted, number_of_calculations ,number_of_results);

}

void task3()
{
    bool is_integral_1 = false, is_weighted = true;
    int number_of_calculations = 1, number_of_results = 2;
    MCMC_integration(is_integral_1, is_weighted, number_of_calculations, number_of_results);
}

int
run(int argc, char *argv[])
{
    //task1();
    task3();
}