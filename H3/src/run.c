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

#define DIM3 3
# define PI 3.14159265358979323846


double potential(double * walker_coordinates);
double weight_factor_6d(double *walker_coordinates, double delta_tau, double ET);
void diffusion(double **walker_coordinates, int Number_of_walkers, gsl_rng * r, double delta_tau);
int spawn_kill_6d(double **walker_coordinates, int *array_of_death,int Number_of_walkers, gsl_rng *r, double delta_tau, double ET);
void death_row(double **old_walker_coordinates,double **new_walker_coordinates,int *array_of_death,int old_number_of_walkers);
void initialize_positions(double **walker_coordinates, int number_of_walkers, gsl_rng *r);
double dmc_6dim(int Number_of_steps,int N_eq_steps, int N0_walkers, double ET, double delta_tau, double gamma);



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

    int task=2;
    
    if(task ==1)
    {
        int N_steps = 10000, N_eq_steps=2000, N_0_walkers=200, save_cutoff=N_steps-5000;

        double gamma = 0.5, ET=0.5, delta_tau=0.02;
        int *N_walker_vec = malloc(sizeof(int)*N_steps+1);
        double  *E_T_vec = malloc(sizeof(double)*N_steps+1);
        ET = restructured_DMC(N_steps, N_eq_steps,save_cutoff, N_0_walkers, gamma, delta_tau, E_T_vec);

        printf("final ET= %f\n", ET);
    }

    if(task ==2)
    {
        //N_steps =10k stabil teori: fel när E average nollställs
        int N_steps = 2*1e5, N_eq_steps=10000, N_0_walkers=1000, save_cutoff=N_steps-3500;

        double gamma = 0.5, ET=-3, delta_tau=0.01;
        int *N_walker_vec = malloc(sizeof(int)*N_steps+1);
        double  *E_T_vec = malloc(sizeof(double)*N_steps+1);

        ET = dmc_6dim(N_steps, N_eq_steps,N_0_walkers, ET, delta_tau,gamma); 
        printf("final ET= %f\n", ET);
    }

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

        if(time_step<1)
        {
            E_average= 0.5;
        }
        else if(time_step <N_eq_steps|| time_step==N_eq_steps)
        {
            for(int step=0; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average/=time_step;
        }
        else
        {
            for(int step=N_eq_steps; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average/=(time_step-N_eq_steps);
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




//task 2






double dmc_6dim(int Number_of_steps,int N_eq_steps, int N0_walkers, double ET_0, double delta_tau, double gamma)
{
    gsl_rng * r;
    int Number_of_walkers= N0_walkers, max_number_walkers=1e6, saved_walkers=0;
    //int equilibration_steps=N_eq_steps;


    //what time step to save after for distribution, max number of positions to save
    int max_save =1e6;

    double *E_T_vector = calloc(sizeof(double), Number_of_steps);
    double *evolution_walker_vec = calloc(sizeof(double), Number_of_steps);

    r = init_random_num_generator();
    double random_number = 0;
    double E_T= ET_0;
    E_T_vector[0]=E_T;
    double **big_array = create_2D_array(max_number_walkers, 6);
    
    initialize_positions(big_array, Number_of_walkers,r);


    for(int time_step=0; time_step< Number_of_steps; ++time_step)
    {
        //printf("timestep = %d\n", time_step);
        double **coordinate_handling = create_2D_array(Number_of_walkers, 6);
        int * array_of_death = calloc(sizeof(int), Number_of_walkers);
        double new_number_of_walkers = 0;
        E_T = E_T_vector[time_step];

        //saving walkers in coordinate handling and cleaning big array
        //maybe break out this part
        for(int walker=0; walker<Number_of_walkers; ++walker)
        {
            for(int coordinate=0; coordinate<6; ++coordinate)
            {
                //if(time_step<1){printf("coord = %f\n", big_array[walker][coordinate]);}
                
                //printf("coord = %f\n",big_array[walker][coordinate]);
                coordinate_handling[walker][coordinate]= big_array[walker][coordinate];
                big_array[walker][coordinate]=0;
            }
        }

        //displace walkers
        //if(time_step<1){printf("old coord %f\n",coordinate_handling[10][3]);}
        
        diffusion(coordinate_handling, Number_of_walkers, r,delta_tau);
        //if(time_step<1){printf("new coord %f\n",coordinate_handling[10][3]);}
        
        //calculate which to kill
        new_number_of_walkers= spawn_kill_6d(coordinate_handling, array_of_death,Number_of_walkers,r,delta_tau,E_T);

        //kill and copy
        death_row(coordinate_handling, big_array, array_of_death, Number_of_walkers);

        //calculate new energy

        destroy_2D_array(coordinate_handling, Number_of_walkers);
        Number_of_walkers = new_number_of_walkers;
        evolution_walker_vec[time_step]= (double) Number_of_walkers;
        printf("new number of walkers = %d\n", Number_of_walkers);

        double E_average  =0.0, new_ET=0.0;

        if(time_step<1)
        {
            E_average= -3;
        }
        else if(time_step <N_eq_steps|| time_step==N_eq_steps)
        {
            for(int step=0; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average/=time_step;
        }
        else
        {
            for(int step=N_eq_steps; step<time_step; ++step)
            {
                E_average += E_T_vector[step];
            }
            E_average/=(time_step-N_eq_steps);
        }

        new_ET = E_average - gamma *log((double) Number_of_walkers/N0_walkers);
        E_T_vector[time_step+1]= new_ET;
        printf("new energy = %f\n", new_ET);
        printf("time step %d\n", time_step);

        free(array_of_death);
    }

    char filename_ET_vec[] = {"./csv/task2_ET_vec.csv"};
    char filename_evolution_walkers[] = {"./csv/task2_evolution_walkers.csv"};

    bool open_with_write = true;

    save_vector_to_csv(E_T_vector, Number_of_steps, filename_ET_vec, true);
    save_vector_to_csv(evolution_walker_vec, Number_of_steps, filename_evolution_walkers, true);

    
    
    double final_energy= E_T_vector[Number_of_steps];
    
    free(E_T_vector);
    return final_energy;


}


void initialize_positions(double **walker_coordinates, int number_of_walkers, gsl_rng *r)
{

    double random_number = 0.0; 
    double *rand_vec = calloc(sizeof(double), 6);
    double theta1 =0, phi1 =0,r1=0;
    double theta2 =0, phi2 =0,r2=0;
    double x1=0, y1=0, z1=0, x2=0, y2=0, z2=0;
    

    for(int walker = 0; walker < number_of_walkers; ++walker)
    {
        for(int coord=0; coord< 6; ++coord)
        {
            random_number = gsl_ran_flat(r, 0.0, 1);
            rand_vec[coord]= random_number;
        }

        r1 = 0.7 + rand_vec[0], theta1 = acos(2*rand_vec[1] -1), phi1 = 2.0*PI *rand_vec[2];
        r2 = 0.7 + rand_vec[3], theta2 = acos(2*rand_vec[4] -1), phi2 = 2.0*PI *rand_vec[5];

        x1 = r1 *sin(theta1)*cos(phi1), y1=r1*sin(theta1)*sin(phi1), z1= r1*cos(theta1);
        x2 = r2 *sin(theta2)*cos(phi2), y2=r2*sin(theta2)*sin(phi2), z2= r2*cos(theta2);
        //printf("x1 = %f\n", x1);
        //printf("y1= %f\n", y1);
        //printf("z1 = %f\n", z1);


        walker_coordinates[walker][0]=x1, walker_coordinates[walker][1]=y1; walker_coordinates[walker][2]=z1;
        walker_coordinates[walker][3]=x2, walker_coordinates[walker][4]=y2; walker_coordinates[walker][5]=z2;
        
    }

    free(rand_vec);
}


//TODO something fucked up
int spawn_kill_6d(double **walker_coordinates, int *array_of_death,int Number_of_walkers, gsl_rng *r, double delta_tau, double ET)
{
    double random_number =0.0; double wf=0.0; 
    int new_number_walkers =0, m_spawn =0;
    //printf("ET%f\n", ET);
    

    for(int walker=0; walker<Number_of_walkers; ++walker)
    {
        wf = weight_factor_6d(walker_coordinates[walker], delta_tau, ET);
        random_number = gsl_ran_flat(r, 0.0, 1);
        m_spawn  = (int) floor(wf+ random_number);
        if(m_spawn>5){m_spawn=5;}
        //printf("spawn %d\n", m_spawn);
        array_of_death[walker] = m_spawn;
        new_number_walkers+= m_spawn;
    }
    //printf("new number of walkers %d\n", new_number_walkers);
    return new_number_walkers;
}

//function that actually spawns the walkers
/*
old_walkers_coordinates = old number of walkers with updated number of walkers [old_number_of_walkers][6]
new_walker_coordinates = array for new copies to be saved into [new_number_of_walkers][6]
array_of_death = array containing what walkers to spawn/kill [old_number_of_walkers]

*/
void death_row(
    double **old_walker_coordinates,
    double **new_walker_coordinates,
    int *array_of_death,
    int old_number_of_walkers
    )
    {
        int m_spawn=0, total_spawned=0;
        double x1=0,y1=0,z1=0,x2=0,y2=0,z2=0;

        //looping over old walkers
        for(int old_walker=0; old_walker<old_number_of_walkers; ++old_walker)
        {
            
            m_spawn= array_of_death[old_walker];
            //printf("spawn = %d\n",m_spawn );
            //printf("old walker = %d\n", old_walker);
            //printf("old number of walker = %d\n", old_number_of_walkers);
            
            //if we should spawn new copies
            if(m_spawn>0)
            {
                x1 = old_walker_coordinates[old_walker][0];
                y1 = old_walker_coordinates[old_walker][1];
                z1 = old_walker_coordinates[old_walker][2];
                x2 = old_walker_coordinates[old_walker][3];
                y2 = old_walker_coordinates[old_walker][4];
                z2 = old_walker_coordinates[old_walker][5];
                
                
                //looping over instances that should be spawned
                for(int new_walker=0; new_walker< m_spawn; ++new_walker)
                {
                    new_walker_coordinates[new_walker+ total_spawned][0] = x1;
                    new_walker_coordinates[new_walker+ total_spawned][1] = y1;
                    new_walker_coordinates[new_walker+ total_spawned][2] = z1;
                    new_walker_coordinates[new_walker+ total_spawned][3] = x2;
                    new_walker_coordinates[new_walker+ total_spawned][4] = y2;
                    new_walker_coordinates[new_walker+ total_spawned][5] = z2;
                }
                total_spawned+=m_spawn;
            }

        }
    }

void diffusion(double **walker_coordinates, int Number_of_walkers, gsl_rng * r, double delta_tau)
{
    double random_number = 0.0;

    for(int walker=0; walker < Number_of_walkers; ++walker)
    {
        for(int coordinate=0; coordinate<6; ++coordinate)
        {
            random_number = gsl_ran_gaussian (r, 1.0);
            //updating coodinate of walker

            //printf("    old pos = %f\n",walker_coordinates[walker][coordinate]);
            walker_coordinates[walker][coordinate] += sqrt(delta_tau)*random_number;
            //printf("updated pos = %f\n",walker_coordinates[walker][coordinate]);
        }
    }
}



double weight_factor_6d(double *walker_coordinates, double delta_tau, double ET)
{
    double W =0.0, V=0.0;

    V = potential(walker_coordinates);

    W =  exp(-(V-ET)*delta_tau);

    return W;
}

double potential(double * walker_coordinates)
{
    //variables
    double potential =0., r1=0, r2=0, r12=0;

    //arrays for calculating r1 and r2 lengths
    double *R1 = calloc(sizeof(double), 3),  *R2 = calloc(sizeof(double), 3);

    //splitting up walker coordinates into two different arrays
    for(int coordinate =0; coordinate<3; ++coordinate)
    {
        R1[coordinate] = walker_coordinates[coordinate];
        R2[coordinate] = walker_coordinates[3+coordinate];
    }
    //lengths for potential
    r12 = distance_between_vectors(R1, R2, 3);
    r1 = vector_norm(R1, 3), r2 = vector_norm(R2,3);

    potential = -2./r1 - 2./r2 + 1./r12;
    //printf("potential %f\n", potential);
    free(R1), free(R2);

    return potential;
}