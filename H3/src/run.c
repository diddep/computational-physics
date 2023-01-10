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
double diffusion_monte_carlo(int N_steps, int N0_walkers, double gamma, double *ET_vec, int *N_walker_vec, double delta_tau);
double clean_DMC(int N_steps, int N_eq_steps, int N0_walkers, double gamma, double delta_tau, double *E_T_vector);



int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code

    int N_steps = 1000, N_eq_steps=10000, N_0_walkers=200;

    double gamma = 0.5, ET=0.5, delta_tau=0.02;
    int *N_walker_vec = malloc(sizeof(int)*N_steps+1);
    double  *E_T_vec = malloc(sizeof(double)*N_steps+1);
    //ET = diffusion_monte_carlo(N_steps, N_0_walkers, gamma, E_T_vec, N_walker_vec, delta_tau);
    ET = clean_DMC(N_steps, N_eq_steps, N_0_walkers, gamma, delta_tau, E_T_vec);
    printf("final ET= %f\n", ET);

    
    
    return 0;
}


/*
idea:
    create position array before diffusion part len=N_curr
    create array that contains number of walkers to spawn/Kill during reactive loop
     
    [m0,m1,m2,---mn]
    use this array 

*/

double diffusion_monte_carlo(int N_steps, int N0_walkers, double gamma, double *ET_vec, int *N_walker_vec, double delta_tau)
{
    int Number_of_walkers= N0_walkers, max_number_walkers=1e6;
    int highest_number_walker = 1000;
    int equilibration_steps=100; 

    gsl_rng * r;
    r = init_random_num_generator();
    double random_number = 0;
    double E_T=0.5;
    double x_start = -5.0, x_end=5.0;
    double x_placement_increment = (x_end-x_start)/((double) N0_walkers);

    //not so elegant solution to begin with but setting max number of walkers to 1e6
    //coordinate_array containts x-coords and array_of_death contains info on what to do with elements
    int *array_of_death = malloc(sizeof(int)*max_number_walkers);
    double *coordinate_array = malloc(sizeof(double)*max_number_walkers);

    //initial placement for walkers
    for(int walker=0; walker< N0_walkers; ++walker)
    {
        double x_coord = x_start + x_placement_increment*((double) walker);
        coordinate_array[walker] = x_coord;
        //printf("xcord %f\n", x_coord);
    }

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
        //printf("rand gauss = %f\n", random_number);
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
        E_T = ET_vec[step];
        double W = weight_factor(x_coord,E_T, delta_tau);

        int m_spawn= (int) W+random_number;
        //printf("mspawn= %d\n", m_spawn);
        array_of_death[walker] =m_spawn;

        new_number_of_walkers = new_number_of_walkers+ m_spawn;
      }
      printf("new number of walkers= %d\n", new_number_of_walkers);
      if(new_number_of_walkers>highest_number_walker)
      {
        highest_number_walker = new_number_of_walkers;
      }

      //spawning in the walkers

      for(int walker =0; walker<Number_of_walkers; ++walker)
      {
        int m_spawn = array_of_death[walker];

        //clearing memory to save in next iteration
        array_of_death[walker]=0;
        double x_coord = coordinate_handling[walker];
        //printf("walker= %d\n", walker);
        //printf("numbwalker= %d\n", Number_of_walkers);
        //printf("    step= %d\n", step);
        

        if(m_spawn>0)
        {
            printf("spawning %d\n", m_spawn);
            for(int spawn=0; spawn< m_spawn; ++spawn)
            {
                coordinate_array[spawn+ walker] =x_coord;
            }
        }
      }
        Number_of_walkers= new_number_of_walkers;
        double New_ET = ET_vec[step] - gamma *log((double) Number_of_walkers/((double) N0_walkers));
        ET_vec[step+1] = New_ET;
        printf("new ET = %f\n", New_ET);
    }

    return ET_vec[N_steps+1];
}

double clean_DMC(int N_steps, int N_eq_steps, int N0_walkers, double gamma, double delta_tau, double *E_T_vector)
{
    int Number_of_walkers= N0_walkers, max_number_walkers=1e6;
    int highest_number_walker = 10000;
    int equilibration_steps=100;

    gsl_rng * r;
    r = init_random_num_generator();
    double random_number = 0;
    double E_T= 0.5;
    E_T_vector[0]=0.5;

    double x_start = -5.0, x_end=5.0;
    double x_placement_increment = (x_end-x_start)/((double) N0_walkers);

    //int *array_of_death = malloc(sizeof(int)*max_number_walkers);
    double *coordinate_array = malloc(sizeof(double)*max_number_walkers);


    //initial placement for walkers
    for(int walker=0; walker< N0_walkers; ++walker)
    {
        double x_coord = x_start + x_placement_increment*((double) walker);
        coordinate_array[walker] = x_coord;
        //printf("xcord %f\n", x_coord);
    }

    for(int time_step=0; time_step<N_steps; ++time_step)
    {
        printf("numbwalk %d\n", Number_of_walkers);
        double *coordinate_handling = malloc(sizeof(double)*Number_of_walkers);
        int new_number_of_walkers=0;

        //saving positions from previous iteration and clearing coordinate_array to save results from this iteration

        for(int walker=0; walker < Number_of_walkers; ++walker)
        {
            coordinate_handling[walker] = coordinate_array[walker];
            //printf("coord= %f\n", coordinate_array[walker]);
            //printf("    number of walkers= %d\n", Number_of_walkers);
            coordinate_array[walker] = 0.0;
        }

        for(int old_walker=0; old_walker<Number_of_walkers; ++old_walker)
        {
            
            random_number =  gsl_ran_gaussian (r, 1.0);
        
            double new_coordinate = coordinate_handling[old_walker] + sqrt(delta_tau)*random_number;
            //printf("new coord = %f\n", new_coordinate);

            E_T= E_T_vector[time_step];
            //printf("ET %f\n", E_T);
            double W = weight_factor(new_coordinate, E_T, delta_tau);
            
            random_number = gsl_ran_flat(r, 0.0, 1.0);
            //printf("    new coord %f\n", new_coordinate);
            //printf("weight %f\n", W);
            int m_spawn = (int) W+ random_number;

            printf("    spawn= %d\n", m_spawn);
            new_number_of_walkers =new_number_of_walkers+ m_spawn;

            //spawn new walkers
            if(m_spawn>0)
            {
                for(int new_walker=0; new_walker<m_spawn; ++new_walker)
                {
                    //printf("Old walker %d \n", old_walker);
                    //printf("New walker %d\n", new_walker);
                    coordinate_array[old_walker+new_walker]= new_coordinate;
                }
            }
        }

        Number_of_walkers = new_number_of_walkers;
        //printf("    number of walkers %d\n", Number_of_walkers);
        // calculating new energy

        double E_average = 0.0;
        double E_T_new = 0.0;
        if(time_step< N_eq_steps)
        {
            for(int step=0; step<time_step; ++step)
            {
                E_average = E_average+ E_T_vector[step];
            }
            E_average /= (double) time_step+1;
            //printf("E_average eq= %f\n", E_average);
        }
        else
        {
            for(int step= N_eq_steps-1; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average /= (double) time_step +1 - N_eq_steps;
            //printf("E_average prod = %f\n", E_average);
        }
        E_T_new = E_average - gamma * log((double) Number_of_walkers/N0_walkers);
        E_T_vector[time_step+1] = E_T_new;
        //printf( "E_T new %f\n", E_T_new);
        free(coordinate_handling);
    }

    free(coordinate_array);
    return E_T_vector[N_steps];
}


double weight_factor(double x_coordinate, double E_T, double delta_tau)
{
    double potential=0.0, weight_factor=0.0, exp_term=0.0;

    //TODO abs coord?
    exp_term = (1- exp(-x_coordinate));
    potential = 0.5* exp_term*exp_term;

    weight_factor = exp(-(potential-E_T)*delta_tau);
    printf("w= %f \n", weight_factor);
    
    return weight_factor;
}

//double calculate_energy(int current_number_of_walkers, int inital_number_of_walker, double gamma)
