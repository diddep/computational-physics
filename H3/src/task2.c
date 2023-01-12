
// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <time.h>
// #include <unistd.h>
// #include <limits.h>
// #include <string.h>


// #include "tools.h"

// #define DIM3 3
// # define PI 3.14159265358979323846


// double potential(double * walker_coordinates);
// double weight_factor_6d(double *walker_coordinates, double delta_tau, double ET);
// void diffusion(double **walker_coordinates, int Number_of_walkers, gsl_rng * r, double delta_tau);
// int spawn_kill_6d(double **walker_coordinates, int *array_of_death,int Number_of_walkers, gsl_rng *r, double delta_tau, double ET);
// void death_row(double **old_walker_coordinates,double **new_walker_coordinates,int *array_of_death,int old_number_of_walkers);
// void initialize_positions(double **walker_coordinates, int number_of_walkers, gsl_rng *r);
// double dmc_6dim(int Number_of_steps,int N_eq_steps, int N0_walkers, double ET, double delta_tau, double gamma);


// int
// run(
//     int argc,
//     char *argv[]
//    )
// {
//     // Write your code here
//     // This makes it possible to test
//     // 100% of you code

//     int N_steps = 1000, N_eq_steps=100, N_0_walkers=1000, save_cutoff=N_steps-3500;

//     double gamma = 0.5, ET=-3, delta_tau=0.01;
//     int *N_walker_vec = malloc(sizeof(int)*N_steps+1);
//     double  *E_T_vec = malloc(sizeof(double)*N_steps+1);

//     ET = dmc_6dim(N_steps, N_eq_steps,N_0_walkers, ET, delta_tau,gamma); 
//     printf("final ET= %f\n", ET);

//     return 0;
// }





// /*
// plan:
//     walkers i 2D array #walkers X 6

//     takes in 2D arrays

//     functions needed:
//         placement
//             -initial placement
//         displacement
//             takes 2d array and displaces walkers
//         weight
//             takes 2d array and calculates killing / weights
//         kill/spawn
//             takes in 2 2d arrays, one with positions and one to save new walkers in

//         this will make it easier if we want to implement drift later
// */

// double dmc_6dim(int Number_of_steps,int N_eq_steps, int N0_walkers, double ET_0, double delta_tau, double gamma)
// {
//     gsl_rng * r;
//     int Number_of_walkers= N0_walkers, max_number_walkers=1e4, saved_walkers=0;
//     //int equilibration_steps=N_eq_steps;


//     //what time step to save after for distribution, max number of positions to save
//     int max_save =1e6;

//     double *E_T_vector = calloc(sizeof(double), Number_of_steps);
//     double *evolution_walker_vec = calloc(sizeof(double), Number_of_steps);

//     r = init_random_num_generator();
//     double random_number = 0;
//     double E_T= ET_0;
//     E_T_vector[0]=E_T;

    

//     double **big_array = create_2D_array(max_number_walkers, 6);


//     for(int time_step=0; time_step< Number_of_steps; ++time_step)
//     {
//         double **coordinate_handling = create_2D_array(Number_of_walkers, 6);
//         int * array_of_death = calloc(sizeof(int), Number_of_walkers);
//         double new_number_of_walkers = 0;
//         E_T = E_T_vector[time_step];

//         //saving walkers in coordinate handling and cleaning big array
//         //maybe break out this part
//         for(int walker=0; walker<Number_of_walkers; ++walker)
//         {
//             for(int coordinate=0; coordinate<6; ++coordinate)
//             {
//                 coordinate_handling[walker][coordinate]= big_array[walker][coordinate];
//                 big_array[walker][coordinate]=0;
//             }
//         }

//         //displace walkers
//         diffusion(coordinate_handling, Number_of_walkers, r,delta_tau);

//         //calculate which to kill
//         new_number_of_walkers= spawn_kill_6d(coordinate_handling, array_of_death,Number_of_walkers,r,delta_tau,ET);

//         //kill and copy
//         death_row(coordinate_handling, big_array, array_of_death, Number_of_walkers);

//         //calculate new energy

//         Number_of_walkers = new_number_of_walkers;
//         evolution_walker_vec[time_step]= (double) Number_of_walkers;

//         double E_average  =0.0, new_ET=0.0;

//         if(time_step<1)
//         {
//             E_average= 0.5;
//         }
//         else if(time_step <N_eq_steps|| time_step==N_eq_steps)
//         {
//             for(int step=0; step<time_step; ++step)
//             {
//                 E_average += E_T_vector[step];
//             }
//             E_average/=time_step;
//         }
//         else
//         {
//             for(int step=N_eq_steps; step<time_step; ++step)
//             {
//                 E_average += E_T_vector[step];
//             }
//             E_average/=(time_step+N_eq_steps);
//         }

//         new_ET = E_average - gamma *log((double) Number_of_walkers/N0_walkers);
//         E_T_vector[time_step+1]= new_ET;
//         //printf("new energy = %f\n", new_ET);

//         free(coordinate_handling), free(array_of_death);
//     }

//     char filename_distribution[200], filename_ET_vec[200], filename_evolution_walkers[200];
//     char cwd[200], buf[200];

//     int written = snprintf(buf, 200, "%s", cwd);

//     snprintf(buf + written, 200 - written, "./csv/task2_distribution.csv");
//     strcpy(filename_distribution, buf);

//     snprintf(buf + written, 200 - written, "./csv/task2_ET_vec.csv");
//     strcpy(filename_ET_vec, buf);

//     snprintf(buf + written, 200 - written, "./csv/task2_evolution_walkers.csv");
//     strcpy(filename_evolution_walkers, buf);

//     bool open_with_write = true;

//     //save_vector_to_csv(coordinate_array, Number_of_walkers, filename_distribution, true);
//     //save_vector_to_csv(save_coordinates, saved_walkers, filename_distribution, true);
//     save_vector_to_csv(E_T_vector,Number_of_steps, filename_ET_vec, true);

//     save_vector_to_csv(evolution_walker_vec, Number_of_steps, filename_evolution_walkers, true);
    
//     //printf("total number of saved walkers%d\n", saved_walkers);

//     //free(save_coordinates), free(coordinate_array);

//    return E_T_vector[Number_of_steps];


// }


// void initialize_positions(double **walker_coordinates, int number_of_walkers, gsl_rng *r)
// {

//     double random_number = 0.0; 
//     double *rand_vec = calloc(sizeof(double), 6);
//     double theta1 =0, phi1 =0,r1=0;
//     double theta2 =0, phi2 =0,r2=0;
//     double x1=0, y1=0, z1=0, x2=0, y2=0, z2=0;
    

//     for(int walker = 0; walker < number_of_walkers; ++walker)
//     {
//         for(int coord=0; coord< 6; ++coord)
//         {
//             random_number = gsl_ran_flat(r, 0.0, 1);
//             rand_vec[coord]= random_number;
//         }

//         r1 = 0.7 + rand_vec[0], theta1 = acos(2*rand_vec[1] -1), phi1 = 2.0*PI *rand_vec[2];
//         r2 = 0.7 + rand_vec[3], theta2 = acos(2*rand_vec[4] -1), phi2 = 2.0*PI *rand_vec[5];

//         x1 = r1 *sin(theta1)*cos(phi1), y1=r1*sin(theta1)*sin(phi1), z1= r1*cos(theta1);
//         x2 = r2 *sin(theta2)*cos(phi2), y2=r2*sin(theta2)*sin(phi2), z2= r2*cos(theta2);

//         walker_coordinates[walker][0]=x1, walker_coordinates[walker][1]=y1; walker_coordinates[walker][2]=z1;
//         walker_coordinates[walker][3]=x2, walker_coordinates[walker][4]=y2; walker_coordinates[walker][5]=z2;
        
//     }

//     free(rand_vec);
// }


// int spawn_kill_6d(double **walker_coordinates, int *array_of_death,int Number_of_walkers, gsl_rng *r, double delta_tau, double ET)
// {
//     double random_number =0.0; double wf=0.0; 
//     int new_number_walkers =0, m_spawn =0;

//     for(int walker=0; walker<Number_of_walkers; ++walker)
//     {
//         wf = weight_factor_6d(walker_coordinates[walker], delta_tau, ET);
//         random_number = gsl_ran_flat(r, 0.0, 1);
//         m_spawn  = (int) floor(wf+ random_number);
//         array_of_death[walker]= m_spawn;
//         new_number_walkers+= m_spawn;
//     }
// }

// //function that actually spawns the walkers
// /*
// old_walkers_coordinates = old number of walkers with updated number of walkers [old_number_of_walkers][6]
// new_walker_coordinates = array for new copies to be saved into [new_number_of_walkers][6]
// array_of_death = array containing what walkers to spawn/kill [old_number_of_walkers]

// */
// void death_row(
//     double **old_walker_coordinates,
//     double **new_walker_coordinates,
//     int *array_of_death,
//     int old_number_of_walkers
//     )
//     {
//         int m_spawn=0, total_spawned=0;
//         double x1=0,y1=0,z1=0,x2=0,y2=0,z2=0;

//         //looping over old walkers
//         for(int old_walker=0; old_number_of_walkers; ++old_walker)
//         {
//             //if we should spawn new copies
//             if(m_spawn>0)
//             {
//                 x1 = old_walker_coordinates[old_walker][0];
//                 y1 = old_walker_coordinates[old_walker][1];
//                 z1 = old_walker_coordinates[old_walker][2];
//                 x2 = old_walker_coordinates[old_walker][3];
//                 y2 = old_walker_coordinates[old_walker][4];
//                 z2 = old_walker_coordinates[old_walker][5];
                
                
//                 //looping over instances that should be spawned
//                 for(int new_walker=0; new_walker< m_spawn; ++new_walker)
//                 {
//                     new_walker_coordinates[new_walker+ total_spawned][0] = x1;
//                     new_walker_coordinates[new_walker+ total_spawned][1] = y1;
//                     new_walker_coordinates[new_walker+ total_spawned][2] = z1;
//                     new_walker_coordinates[new_walker+ total_spawned][3] = x2;
//                     new_walker_coordinates[new_walker+ total_spawned][4] = y2;
//                     new_walker_coordinates[new_walker+ total_spawned][5] = z2;
//                 }
//                 total_spawned+=m_spawn;
//             }

//         }
//     }

// void diffusion(double **walker_coordinates, int Number_of_walkers, gsl_rng * r, double delta_tau)
// {
//     double random_number = 0.0;

//     for(int walker=0; walker < Number_of_walkers; ++walker)
//     {
//         for(int coordinate=0; coordinate<6; ++coordinate)
//         {
//             random_number = gsl_ran_gaussian (r, 1.0);
//             //updating coodinate of walker
//             walker_coordinates[walker][coordinate] += sqrt(delta_tau)*random_number;
//         }
//     }
// }



// double weight_factor_6d(double *walker_coordinates, double delta_tau, double ET)
// {
//     double W =0.0, V=0.0;

//     V = potential(walker_coordinates);

//     W =  exp(-(V-ET)*delta_tau);

//     return W;
// }

// double potential(double * walker_coordinates)
// {
//     //variables
//     double potential =0., r1=0, r2=0, r12=0;

//     //arrays for calculating r1 and r2 lengths
//     double *R1 = calloc(sizeof(double), 3),  *R2 = calloc(sizeof(double), 3);

//     //splitting up walker coordinates into two different arrays
//     for(int coordinate =0; coordinate<3; ++coordinate)
//     {
//         R1[coordinate] = walker_coordinates[coordinate];
//         R2[coordinate] = walker_coordinates[3+coordinate];
//     }
//     //lengths for potential
//     r12 = distance_between_vectors(R1, R2, 3);
//     r1 = vector_norm(R1, 3), r2 = vector_norm(R2,3);

//     potential = -2./r1 - 2./r2 + 1./r12;
//     free(R1), free(R2);

//     return potential;
// }