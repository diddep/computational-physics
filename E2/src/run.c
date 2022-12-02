/* 
 * Help routines for E2
 *
 * E2code.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N_PARTICLES 32
#define PI 3.141592653589

#include <stdbool.h>
#include "tools.h"


/*
 * trans_matrix[N_PARTICLES][N_PARTICLES]: empty allocated array which
 * will be filled with sine transformation matrix
 * N_PARTICLES: number of particles in system
 */
void construct_transformation_matrix(
    double trans_matrix[N_PARTICLES][N_PARTICLES], int n_particles)
{
    double factor = 1 / ((double)n_particles + 1);
    for(int i = 0; i < n_particles; i++){
	for(int j = 0; j < n_particles; j++){
	    trans_matrix[i][j] = sqrt(2 * factor)
				 * sin((j + 1) * (i + 1) * PI * factor);
	}
    }
}

/*
 * Transformation matrix constucted in above function
 * q cartesian coordinate of particles
 * Q output normal modes coordinate
 * N_PARTICLES is number of particles in system
 */
void transform_to_normal_modes(double trans_matrix[N_PARTICLES][N_PARTICLES],
			       int n_particles,
			       double *q, double *Q)
{
    for(int i = 0; i < n_particles; i++){
	double sum = 0;
	for(int j = 0; j < n_particles; j++){
	    sum += q[j] * trans_matrix[i][j];
	}
	Q[i] = sum;
    }
}
/* Calculating the acceleration for array of positions*/
void calc_acc(double *a, double *u, double *m, double kappa, double alpha , int size_of_u)
{
  
    /* Calculating the acceleration on the boundaries */

    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[0] *= 1 + alpha * u[1];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    a[size_of_u - 1] *= 1 + alpha * ( - u[size_of_u - 2]);
    
    /* Calculating the acceleration of the inner points */
    for (int i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
        a[i] *= (1 + alpha * (u[i + 1] - u[i-1]));
    }
}

/* Set initial conditions of input arrays, particles start at zero, velocities so all energy is in mode 1. All masses is 1.*/
void set_initial_condition(double *v, double *q, double *m)
{   
    double E0 = N_PARTICLES;
    double Q0=0;
    double P0 = sqrt(2*E0); 

    for(int ix = 0; ix < N_PARTICLES; ix++)
    {
        q[ix] = 0;
        m[ix] = 1;
        v[ix] = sqrt((double) 2 / (N_PARTICLES + 1)) * sqrt(m[ix]) * P0 * sin( (ix + 1) * PI / (N_PARTICLES + 1)) / m[ix];
       
    }
}
/* Take a velocity verlet timestep such that states go from t -> t + dt*/
void velocity_verlet_timestep(double dt, double *v, double *q, double *a, double *m, double kappa, double alpha, int n_particles)
{
    // v(t+dt/2)
    for(int ix = 0; ix < n_particles; ix++)
    {
        v[ix] += 0.5 * dt * a[ix];
    }

    // q(t+dt)
    for(int ix = 0; ix < n_particles; ix++)
    {
        q[ix] += dt * v[ix];
    }

    // a(t+dt)
    calc_acc((double*) a, (double*) q, (double*) m, (double) kappa, (double) alpha , (int) n_particles);
    
    // v(t+dt)
    for(int ix = 0; ix < n_particles; ix++)
    {
        v[ix] += 0.5 * dt * a[ix];
    }
}

/* Function to run full velocity verlet simulation */
void velocity_verlet(int n_timesteps, int timestep_interval, double dt, double *v, double *q, \
 double *m, double kappa, double alpha)
{
    // Initiate file names
    char filename_result[] = {"verlet_results.csv"};
    char filename_positions[] = {"verlet_positions.csv"};
    char filename_velocities[] = {"verlet_velocities.csv"};
    char filename_param[] = {"verlet_params.csv"};
    char filename_energies[] = {"verlet_energies.csv"};
    char filename_energy_average[] = {"verlet_energy_averages.csv"};
    bool is_write;
     
    // Initiate arrays
    double T = 1;
    double a[N_PARTICLES];
    double E[N_PARTICLES];
    double Q[N_PARTICLES];
    double P[N_PARTICLES];
    double omega[N_PARTICLES];
    double E_tot_dt[N_PARTICLES];
    double E_time_average[N_PARTICLES];
    
    // Set constants, time will increase when the current timestep modulo timestep_interval is zero
    int n_time_rows = n_timesteps/timestep_interval*dt;
    int count = 0;
    
    // Setting energies to zero and eigenfrequencies for each particle according to formula
    for(int ix = 0; ix < N_PARTICLES; ix++)
    {
        omega[ix] = 2 * sqrt( kappa / m[ix] ) * sin( (ix + 1) * PI / (2 * (N_PARTICLES + 1) ) );
        E_tot_dt[ix] = 0; E_time_average[ix] = 0;
    }

    // Initiate transformation matrix
    double trans_matrix[N_PARTICLES][N_PARTICLES];
    construct_transformation_matrix((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES);

    // Get acceleration
    calc_acc((double*) a, (double*) q, (double*) m, (double) kappa, (double) alpha , (int) N_PARTICLES);
     
    for(int tx = 1; tx < n_timesteps + 1; tx++)
    {
        // t -> t + dt
        velocity_verlet_timestep((double) dt, (double*) v, (double*) q, (double*) a, (double*) m, (double) kappa, (double) alpha, (int) N_PARTICLES);

        if(tx % timestep_interval == 0 || tx == 1)
        {
            //Tranform to normal modes
            transform_to_normal_modes((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES, q, Q);
            transform_to_normal_modes((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES, v, P);

            // Calculate current and time-averaged energy 
            for(int ix = 0; ix < N_PARTICLES; ix++)
            {
                E[ix] = 0.5 * ( pow(P[ix], 2) + pow(omega[ix], 2) * pow(Q[ix], 2) );

                if( count != 0)
                {   
                    // Get current time-factor for time-averaged energy
                    T = ( count * timestep_interval * dt);

                    // Get current time-factor for time-averaged energy
                    E_tot_dt[ix] +=  E[ix] * dt;
                }
                // Set time-averaged energy
                E_time_average[ix] = E_tot_dt[ix]/T;
            }
            // Setting bool to true to open csv files with write, false to open with append
            if (tx == 1) { is_write = true; } else { is_write = false; }
            
            // Saving values to csv files.
            double result_vec[] = {tx * dt};
            save_vector_to_csv(result_vec, 1, filename_result, is_write);
            save_vector_to_csv(q, N_PARTICLES, filename_positions, is_write);
            save_vector_to_csv(v, N_PARTICLES, filename_velocities, is_write);
            save_vector_to_csv(E, N_PARTICLES, filename_energies, is_write);
            save_vector_to_csv(E_time_average, N_PARTICLES, filename_energy_average, is_write);

            // Incrementing counter
            count++;
        }
    }
    // Save simulation parameters dt and alpha in vector
    double param_vec[] = {dt, alpha};
    save_vector_to_csv(param_vec, 2, filename_param, true);
}

/* "Main function" of program*/
int run()
{
    // Initiating arrays
    double q[N_PARTICLES];
    double v[N_PARTICLES];
    double m[N_PARTICLES];

    // Setting constants
    double kappa = 1, alpha = 0.1;

    // Set initial conditions for v, q and m (mass is just one)
    set_initial_condition((double*) v, (double*) q, (double*) m);
    
    // Simulation parameters. timestep_interval governs interval spacing between calculating energies and saving to csv.
    int end_time = 1e6; double dt = 1e-1;
    int n_timesteps = end_time / dt;
    int timestep_interval = 1e3;

    // Running the velocity verlet function
    velocity_verlet((int) n_timesteps, (int) timestep_interval, (double) dt, (double*) v, (double*) q, \
                    (double*) m, (double) kappa, (double) alpha);

    return 0;
}
