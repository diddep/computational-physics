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

void update_coordinates(double *x_coordinates, int number_walkers, gsl_rng * r, double delta_tau);
int spawn_kill(double *x_coordinates, int * array_of_death, int number_walkers, gsl_rng * r, double delta_tau, double ET);
double restructured_DMC(int N_steps, int N_eq_steps,int save_cutoff, int N0_walkers, double gamma, double delta_tau, double *E_T_vector);


int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code

    int N_steps = 5*5000, N_eq_steps=1000, N_0_walkers=200, save_cutoff=N_steps-2500;

    double gamma = 0.5, ET=0.5, delta_tau=0.02;
    int *N_walker_vec = malloc(sizeof(int)*N_steps+1);
    double  *E_T_vec = malloc(sizeof(double)*N_steps+1);
    ET = restructured_DMC(N_steps, N_eq_steps,save_cutoff, N_0_walkers, gamma, delta_tau, E_T_vec);

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


double restructured_DMC(int N_steps, int N_eq_steps,int save_cutoff, int N0_walkers, double gamma, double delta_tau, double *E_T_vector)
{
    // 
    int Number_of_walkers= N0_walkers, max_number_walkers=1e3, saved_walkers=0;
    //int equilibration_steps=N_eq_steps;

    //what time step to save after for distribution, max number of positions to save
    int max_save =1e6;

    gsl_rng * r;
    r = init_random_num_generator();
    double random_number = 0;
    double E_T= 0.5;
    E_T_vector[0]=E_T;

    double x_start = -5.0, x_end=5.0;
    double x_placement_increment = (x_end-x_start)/((double) N0_walkers);

    //int *array_of_death = malloc(sizeof(int)*max_number_walkers);
    double *coordinate_array = calloc(sizeof(double),max_number_walkers);
    double *evolution_walker_array = calloc(sizeof(double), N_steps);

    //array to save last iterations in
    double *save_coordinates = calloc(sizeof(double), max_save);


    //initial placement for walkers
    for(int walker=0; walker< N0_walkers; ++walker)
    {
        double x_coord = x_start + x_placement_increment*((double) walker);
        coordinate_array[walker] = x_coord;
        //printf("xcord %f\n", x_coord);
    }

    for(int time_step=0; time_step<N_steps; ++time_step)
    {
        E_T = E_T_vector[time_step];

        evolution_walker_array[time_step]= (double) Number_of_walkers;

        // array for saving positions during run
        double *coordinate_handling = calloc(sizeof(double),Number_of_walkers);
        
        int *array_of_death = calloc(sizeof(int),Number_of_walkers);
        int new_number_of_walkers=0;

        //saving positions from previous iteration and clearing coordinate_array to save results from this iteration

        for(int walker=0; walker < Number_of_walkers; ++walker){coordinate_handling[walker] = coordinate_array[walker];}

        for(int walker=0; walker < max_number_walkers; ++walker){coordinate_array[walker]=0.0;}
        
        //if above save cutoff, save coordinates in large array
        if(time_step>save_cutoff)
        {
            for(int walker =0; walker< Number_of_walkers; ++walker)
            {
                save_coordinates[saved_walkers+walker]= coordinate_handling[walker];
            }
            saved_walkers +=Number_of_walkers;
        }

        //updating walker positions
        update_coordinates(coordinate_handling, Number_of_walkers, r, delta_tau);

        //calculating which walkers to kill

        new_number_of_walkers = spawn_kill(coordinate_handling, array_of_death, Number_of_walkers, r, delta_tau, E_T);

        //printf("new_number of walkers= %d\n", new_number_of_walkers);
        // spawning new walkers

        int tot_spawned =0;
        for(int old_walker =0; old_walker< Number_of_walkers; ++old_walker)
        {
            double xcoord = coordinate_handling[old_walker];
            int m_spawn = array_of_death[old_walker];
            
            if(m_spawn>0)
            {
                for(int new_walker =0; new_walker<m_spawn; ++new_walker)
                {
                    coordinate_array[tot_spawned+ new_walker] = xcoord;
                }
            }
            tot_spawned += m_spawn;
        }


        //calculate new energy

        Number_of_walkers = new_number_of_walkers;
        double E_average  =0.0, new_ET=0.0;

        if(time_step>0)
        {
            for(int step=0; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average/=time_step;
        }
        else
        {
            E_average= 0.5;
        }

        new_ET = E_average - gamma *log((double) Number_of_walkers/N0_walkers);
        E_T_vector[time_step+1]= new_ET;
        //printf("new energy = %f\n", new_ET);

    free(coordinate_handling), free(array_of_death) ;
    }

    //saving arrays

    char filename_distribution[200], filename_ET_vec[200], filename_evolution_walkers[200];
    char cwd[200], buf[200];

    int written = snprintf(buf, 200, "%s", cwd);

    snprintf(buf + written, 200 - written, "./csv/task1_distribution.csv");
    strcpy(filename_distribution, buf);

    snprintf(buf + written, 200 - written, "./csv/task1_ET_vec.csv");
    strcpy(filename_ET_vec, buf);

    snprintf(buf + written, 200 - written, "./csv/evolution_walkers.csv");
    strcpy(filename_evolution_walkers, buf);

    bool open_with_write = true;

    //save_vector_to_csv(coordinate_array, Number_of_walkers, filename_distribution, true);
    save_vector_to_csv(save_coordinates, saved_walkers, filename_distribution, true);
    save_vector_to_csv(E_T_vector,N_steps, filename_ET_vec, true);

    save_vector_to_csv(evolution_walker_array,N_steps, filename_evolution_walkers, true);
    
    printf("total number of saved walkers%d\n", saved_walkers);

    free(save_coordinates), free(coordinate_array);

   return E_T_vector[N_steps];
}

double weight_factor(double x_coordinate, double E_T, double delta_tau)
{
    double potential=0.0, weight_factor=0.0, exp_term=0.0;

    exp_term = (1.- exp(-x_coordinate));
    potential = 0.5* exp_term*exp_term;

    weight_factor = exp(-(potential-E_T)*delta_tau);
    //printf("w= %f \n", weight_factor);
    
    return weight_factor;
}

void update_coordinates(double *x_coordinates, int number_walkers, gsl_rng * r, double delta_tau)
{
    double random_number =0;
    double x_old = 0.0, x_new=0.0;

    // TODO måste vi göra ny rng?
    for(int walker =0; walker<number_walkers; ++walker)
    {
        random_number =  gsl_ran_gaussian (r, 1.0);
        x_old = x_coordinates[walker];
        x_new = x_old+ random_number* sqrt(delta_tau);
        x_coordinates[walker] = x_new;
    }

}

//returns new number of walkers, saves what walkers to copy/kill
int spawn_kill(double *x_coordinates, int * array_of_death, int number_walkers, gsl_rng * r, double delta_tau, double ET)
{
    int m_spawn=0, new_number_of_walkers=0;
    double x_coord=0.0, wf=0.0, random_number=0.0; 

    for(int walker=0; walker< number_walkers; ++walker)
    {
        random_number = gsl_ran_flat(r, 0.0, 1);
        x_coord =x_coordinates[walker];
        wf= weight_factor(x_coord,ET,delta_tau);

        //printf("wf + U = %f\n", wf+random_number);
        m_spawn = (int) floor(wf+ random_number);

        //printf("m spawn = %d\n", m_spawn);
        array_of_death[walker] = m_spawn;
        
        new_number_of_walkers += m_spawn;
    }
    return new_number_of_walkers;
}
