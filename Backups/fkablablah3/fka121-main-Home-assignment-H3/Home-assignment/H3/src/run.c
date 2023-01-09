#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>


#include "tools.h"

double weight_factor(double x_coordinate, double E_T, double delta_tau);
double diffusion_monte_carlo(int N_steps, int N0_walkers, double gamma, double *ET_vec, double N_walker_vec, double delta_tau);



int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code

    int N_steps = 1000, N_0_walkers=100;

    double gamma = 0.5;
    int *N_walker_vec = malloc(sizeof(int)*N_steps), *E_T_vec = malloc(sizeof(int)*N_steps);




    return 0;
}


/*
idea:
    create position array before diffusion part len=N_curr
    create array that contains number of walkers to spawn/Kill during reactive loop
     
    [m0,m1,m2,---mn]
    use this array 

*/

double diffusion_monte_carlo(int N_steps, int N0_walkers, double gamma, double *ET_vec, double N_walker_vec, double delta_tau)
{
    int Number_of_walkers= N0_walkers, max_number_walkers=1e6;
    int highest_number_walker = 1000;

    gsl_rng * r;
    r = init_random_num_generator();
    double random_number = 0;
    double E_T=0;

    //not so elegant solution to begin with but setting max number of walkers to 1e6
    //coordinate_array containts x-coords and array_of_death contains info on what to do with elements
    int *array_of_death = malloc(sizeof(int)*max_number_walkers);
    double *coordinate_array = malloc(sizeof(double)*max_number_walkers);

    //Looping through a specified number of steps for the dmc 
    for(int step =0; step< N_steps; ++step)
    {

        /*
        diffusive part, loop through N_i walkers and displace each randomly
        TODO how to keep track of positions? Different number of walkers each step
        */
       double *coordinate_handling = malloc(sizeof(double)*Number_of_walkers);
       int *array_of_death_handling = malloc(sizeof(int)*Number_of_walkers);
    

       for(int walker=0; walker< Number_of_walkers; ++walker)
       { 
        //random number, updating position and saving in handling vec
        random_number =  gsl_ran_gaussian (r, 1.0);
        double x_coord = coordinate_array[walker] + sqrt(delta_tau)*random_number;
        coordinate_handling[walker] = x_coord;
       }

        //clearing position array for saving new data
       for(int walker=0; walker<max_number_walkers; ++walker)
       {
        coordinate_array[walker] =0.0;
       }

       /*
       reactive part of DMC, looping through array and killing / spawning walkers

       */
      int new_number_of_walkers = 0;
      for(int walker =0; walker< Number_of_walkers; ++walker)
      {
        random_number = gsl_ran_flat(r, 0.0, 1.0);
        double x_coord = coordinate_handling[walker];
        double W = weight_factor(x_coord,E_T, delta_tau);

        int m_spawn= (int) W+random_number;
        array_of_death[walker] =m_spawn;

        new_number_of_walkers =+ m_spawn;
      }
      if(new_number_of_walkers>highest_number_walker)
      {
        highest_number_walker = new_number_of_walkers;
      }

      //spawning in the walkers

      for(int walker =0; walker<Number_of_walkers; ++walker)
      {
        int m_spawn = array_of_death[walker];
        double x_coord = coordinate_handling[walker];

        if(m_spawn>0)
        {
            for(int spawn=0; spawn< m_spawn; ++m_spawn)
            {
                coordinate_array[spawn+ walker];
            }
        }
      }
        Number_of_walkers= new_number_of_walkers;

        double new_E_T = ET_vec[step];
        ET_vec[step+1] = new_E_T;

    }

}

double weight_factor(double x_coordinate, double E_T, double delta_tau)
{
    double potential=0.0, weight_factor=0.0;

    potential = 0.5*(1-exp(-x_coordinate));

    weight_factor = exp(-(potential-E_T)*delta_tau);
    
    return weight_factor;
}