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
void pdf_transformed_points_task2(double *array_of_points, int number_of_points)
{
    for(int ix = 0; ix < number_of_points; ix++)
    {
        array_of_points[ix] = acos(array_of_points[ix])/M_PI;
    }
}

double function_value_task1(double point_to_evaluate)
{
    return point_to_evaluate * (1 - point_to_evaluate);
}

// TODO: Update with new 3d function
double function_value_task3(double point_to_evaluate)
{
    return point_to_evaluate * (1 - point_to_evaluate);
}

double weight_function_task2(double point_to_evaluate)
{
    return sin(M_PI * point_to_evaluate);
}

// TODO: Update with new 3d function
double weight_function_task3(double point_to_evaluate)
{
    return sin(M_PI * point_to_evaluate);
}

void evaluate_integral(double *array_of_points, int number_of_points, double *integral_results, bool is_integral_1, bool is_weighted)
{
    double integral_value = integral_results[0], integral_std = integral_results[1];
    double variance = 0, function_value = 0, weight_function = 1; 

    for(int ix = 0; ix < number_of_points; ix++)
    {

        if(is_integral_1)
        {
            function_value = function_value_task1(array_of_points[ix]);
            if(is_weighted)
            {
                weight_function = weight_function_task2(array_of_points[ix]);
                function_value /= weight_function;
            }
        } 
        else {
            function_value = function_value_task3(array_of_points[ix]);
            if(is_weighted)
            {
                weight_function = weight_function_task3(array_of_points[ix]);
                function_value /= weight_function;
            }
        }

        integral_value += function_value / number_of_points;
        variance += pow(function_value, 2.) / number_of_points \
                    - pow(function_value / number_of_points, 2.); 
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
        int number_of_points = pow((int) 10., (int) ix + 1);
        double *array_of_points = calloc(number_of_points, sizeof(double));
        double lower_bound = 0, upper_bound = 1;

        get_points(array_of_points, number_of_points, lower_bound, upper_bound);
        if(is_integral_1 && is_weighted)
        {
            pdf_transformed_points_task2(array_of_points, number_of_points);
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