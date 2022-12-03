#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "tools.h"
#include "integrals.h"


void get_points(double *array_of_points, int number_of_points, double lower_bound, double upper_bound)
{   
    gsl_rng *r;
    r = init_random_num_generator();

    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[ix] = gsl_ran_flat(r, lower_bound, upper_bound);
    }

    gsl_rng_free(r);
}

// FIXME: Incorrect transform method (see page 9 in MC lecture notes)
// TODO: Make three dimensional?
void pdf_transformed_points_task2(double *array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[ix] = acos(array_of_points[ix])/M_PI;
    }
}

// TODO: Update transform to relevant function
//TODO: Make three dimensional
void pdf_transformed_points_task3(double *array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[ix] = acos(array_of_points[ix])/M_PI;
    }
}

// TODO: Update and make loop through all steps
// TODO: Return vector with length number_of_points
double function_value_task1(double *point_to_evaluate)
{
    double x_coord = point_to_evaluate[0], y_coord = point_to_evaluate[1], z_coord = point_to_evaluate[2];
    return x_coord * (1 - x_coord);
}

// TODO: Update and make loop through all steps
// TODO: Return vector with length number_of_points
double function_value_task3(double *point_to_evaluate)
{
    double x_coord = point_to_evaluate[0], y_coord = point_to_evaluate[1], z_coord = point_to_evaluate[2];
    
    double factor_1 = pow(M_PI,(double) -3./2.);
    double factor_2 = pow(x_coord, 2.) + pow(x_coord * y_coord, 2.)  + pow(x_coord * y_coord * z_coord, 2.);
    double factor_3 = exp(-(pow(x_coord, 2.) + pow(y_coord, 2.) +pow(z_coord, 2.)));
    
    return factor_1 * factor_2 * factor_3;
}

double weight_function_task2(double *point_to_evaluate)
{
    double x_coord = point_to_evaluate[0], y_coord = point_to_evaluate[1], z_coord = point_to_evaluate[2];

    return sin(M_PI * x_coord);
}

double weight_function_task3(double *point_to_evaluate)
{
    double x_coord = point_to_evaluate[0], y_coord = point_to_evaluate[1], z_coord = point_to_evaluate[2];

    double factor_1 = pow(M_PI,(double) -3./2.);
    double factor_3 = exp(-(pow(x_coord, 2.) + pow(y_coord, 2.) +pow(z_coord, 2.)));
    
    return factor_1 * factor_3;
}

// TODO: Make array_of_points carry x, y and z coordinate
void evaluate_integral(double **array_of_points, int number_of_points, double *integral_results, bool is_integral_1, bool is_weighted)
{
    double integral_value = integral_results[0], integral_std = integral_results[1];
    double variance = 0, weight_function = 1; int ndim;
    double *function_value = calloc(number_of_points, sizeof(double));

    if(is_integral_1)
    {
        ndim = 1;
        function_value = function_value_task1(array_of_points);
        if(is_weighted)
        {
            weight_function = weight_function_task2(array_of_points);
            for(int ix = 0; ix < number_of_points; ix++)
            {
                function_value[ix] /= weight_function;
            }
        }
    } else {
        ndim = 3;
        function_value = function_value_task3(array_of_points);
        if(is_weighted)
        {
            weight_function = weight_function_task3(array_of_points);
            for(int ix = 0; ix < number_of_points; ix++)
            {
                function_value[ix] /= weight_function;
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

//TODO: Figure out logical struction so that MCMC_integration can be used for all 4 tasks.
void MCMC_integration(bool is_integral_1, bool is_weighted, int number_of_calculations, int number_of_results)
{

    double **calculation_results = create_2D_array(number_of_calculations, number_of_results);
    char filename_calculation_results[] = {"../task1/calculation_results_task1.csv"};
    char filename_sampled_points[] = {"../task2/sampled_points.csv"};

    for(int ix = 0; ix < number_of_calculations; ix++)
    {   
        double **array_of_points; int number_of_points, ndim;
        if(is_integral_1)
        {
            ndim = 1;
            number_of_points = pow((int) 10., (int) ix + 1);
            // TODO: dimension of array of points need to be flipped, make sync with other functions
            array_of_points = create_2D_array(ndim, number_of_points);
            double lower_bound = 0, upper_bound = 1;
            get_points(array_of_points[0], number_of_points, lower_bound, upper_bound);

            if(is_weighted)
            {
                pdf_transformed_points_task2(array_of_points[0], number_of_points);
            }
        } else {
            ndim = 3;
            number_of_points = 1e4;
            array_of_points = create_2D_array(number_of_points, ndim);

            // TODO: Update lower and upper bounds
            double lower_bound = 0, upper_bound = 1;
            for(int jx = 0; jx < ndim; jx++)
            {
                get_points(array_of_points[jx], number_of_points, lower_bound, upper_bound);
            } 
            if(is_weighted)
                {
                    pdf_transformed_points_task3(array_of_points, number_of_points);
                }
        }
       
        double integral_value = 0, integral_std = 0;
        double integral_results[2]; integral_results[0] = integral_value;  integral_results[1] = integral_std;

        evaluate_integral(array_of_points, number_of_points, integral_results, is_integral_1, is_weighted);

        for(int jx = 0; jx < number_of_results; jx++)
        {
            calculation_results[ix][jx] = integral_results[jx];
        }
        save_vector_to_csv(array_of_points, number_of_points, filename_sampled_points, true);
    }
    
    save_matrix_to_csv(calculation_results, number_of_calculations, number_of_results, filename_calculation_results);
    destroy_2D_array(calculation_results, number_of_calculations);
}

void task1()
{
    bool is_integral_1 = true, is_weighted = true;
    int number_of_calculations = 4, number_of_results = 2;
    MCMC_integration(is_integral_1, is_weighted, number_of_calculations,number_of_results  );

}

// TODO: Decide on whether to have separate function or reuse task1 with different bools
void task3()
{
    bool is_integral_1 = false, is_weighted = true;
    int number_of_calculations = 1, number_of_results = 2;
    MCMC_integration(is_integral_1, is_weighted, number_of_calculations,number_of_results);
}

int
run(int argc, char *argv[])
{
    task1();
    //TODO: Uncomment when ready to run
    //task3();
}