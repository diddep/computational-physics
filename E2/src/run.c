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

void velocity_verlet(int n_timesteps, int timestep_interval, double dt, double *v, double *q, \
 double *m, double kappa, double alpha)
{
    char filename_result[] = {"verlet_results.csv"};
    char filename_positions[] = {"verlet_positions.csv"};
    char filename_velocities[] = {"verlet_velocities.csv"};
    char filename_energies[] = {"verlet_energies.csv"};
    char filename_param[] = {"verlet_params.csv"};
    char filename_energy_average[] = {"verlet_energy_averages.csv"};
    bool is_write;
     
    
    double trans_matrix[N_PARTICLES][N_PARTICLES];
    double E[N_PARTICLES];
    double Q[N_PARTICLES];
    double P[N_PARTICLES];
    double omega[N_PARTICLES];
    
    
    int n_time_rows = n_timesteps/timestep_interval*dt;
    double **E_time_average = create_2D_array(n_time_rows, N_PARTICLES);
    printf("no rows: %d\n", n_time_rows);
    int count = 0;
    

    for(int ix = 0; ix < N_PARTICLES; ix++)
    {
        omega[ix] = 2 * sqrt( kappa / m[ix] ) * sin( (ix + 1) * PI / (2 * (N_PARTICLES + 1) ) );
        for(int txx = 0; txx < n_time_rows; txx++)
        {
            E_time_average[txx][ix] = 0;
        }
    }

    construct_transformation_matrix((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES);

    double a[N_PARTICLES];
    calc_acc((double*) a, (double*) q, (double*) m, (double) kappa, (double) alpha , (int) N_PARTICLES);
     
    for(int tx = 1; tx < n_timesteps + 1; tx++)
    {
        velocity_verlet_timestep((double) dt, (double*) v, (double*) q, (double*) a, (double*) m, (double) kappa, (double) alpha, (int) N_PARTICLES);

        if(tx % timestep_interval == 0 || tx == 1)
        {
            transform_to_normal_modes((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES, q, Q);
            transform_to_normal_modes((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES, v, P);
        
            for(int ix = 0; ix < N_PARTICLES; ix++)
            {
                E[ix] = 0.5 * ( pow(P[ix], 2) + pow(omega[ix], 2) * pow(Q[ix], 2) );
            }

            if(tx == 1)
            {
                is_write = true;
            } else {
                is_write = false;
            }

            if(count != 0)
            {   
                double T = ( count * timestep_interval * dt);
                for(int ixx = 0; ixx < N_PARTICLES; ixx++)
                {
                    if(count == 1)
                    {
                        E_time_average[count][ixx] = E[ixx] * dt;
                    } else 
                    {
                        for(int txx = 0; txx < count; txx++)
                        {
                            E_time_average[count][ixx] += E_time_average[txx][ixx] * dt;;
                        }

                        E_time_average[count][ixx] += E[ixx] * dt;
                        E_time_average[count][ixx] *= 1/T;
                    }                   
                }
            }
            //printf("count: %d\n", count);
            printf("E[%d][0]: %f\n",count, E_time_average[count][0]);
            count++;
            //printf("tx/timestep_interval: %d\n", tx/timestep_interval);

            double result_vec[] = {tx * dt};

            save_vector_to_csv(result_vec, 2, filename_result, is_write);
            save_vector_to_csv(q, N_PARTICLES, filename_positions, is_write);
            save_vector_to_csv(v, N_PARTICLES, filename_velocities, is_write);
            save_vector_to_csv(E, N_PARTICLES, filename_energies, is_write);
            
        }
    }
    double param_vec[] = {dt, alpha};
    save_vector_to_csv(param_vec, 2, filename_param, true);
    save_matrix_to_csv(E_time_average, n_time_rows, N_PARTICLES, filename_energy_average);
    destroy_2D_array(E_time_average, n_timesteps/timestep_interval);
}

int run()
{
    
    double q[N_PARTICLES];
    
    double v[N_PARTICLES];
    double m[N_PARTICLES];
    double kappa = 1, alpha = 0.1;

    set_initial_condition((double*) v, (double*) q, (double*) m);
    
    int end_time = 1e6; double dt = 1e-1;
    int n_timesteps = end_time / dt;
    int timestep_interval = 1e3;

    velocity_verlet((int) n_timesteps, (int) timestep_interval, (double) dt, (double*) v, (double*) q, \
    (double*) m, (double) kappa, (double) alpha);


    // Evolove system in time 

    
    return 0;
}
