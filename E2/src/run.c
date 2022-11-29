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
 * q cartesian coordinate of paricles
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
        v[ix] = sqrt((double) 2 / (N_PARTICLES + 1)) * sqrt(m[ix]) * P0 * sin( ix * PI / (N_PARTICLES + 1)) / m[ix];
       
    }
}

void velocity_verlet(int n_timesteps, double dt, double *v, double *q, \
 double *m, double kappa, double alpha)
{
    char filename_result[] = {"velocity_verlet.csv"};
    char filename_positions[] = {"verlet_pos.csv"};
    bool is_write;
     
    double a[N_PARTICLES];

    calc_acc((double*) a, (double*) q, (double*) m, (double) kappa, (double) alpha , (int) N_PARTICLES);
     
    for(int tx = 1; tx < n_timesteps + 1; tx++)
    {
        // v(t+dt/2)
        for(int ix = 0; ix < N_PARTICLES; ix++)
        {
            v[ix] += 0.5 * dt * a[ix];
        }

        // q(t+dt)
        for(int ix = 0; ix < N_PARTICLES; ix++)
        {
            q[ix] += dt * v[ix];
        }

        // a(t+dt)
        calc_acc((double*) a, (double*) q, (double*) m, (double) kappa, (double) alpha , (int) N_PARTICLES);

        // v(t+dt)
        for(int ix = 0; ix < N_PARTICLES; ix++)
        {
            v[ix] += 0.5 * dt * a[ix];
        }

        double result_vec[] = {tx * dt};
        if(tx == 0)
        {
            is_write = true;
        } else {
            is_write = false;
        }

        save_vector_to_csv(result_vec, 1, filename_result, is_write);
        save_vector_to_csv(q, N_PARTICLES, filename_positions, is_write);
        
    }


}

int run()
{
    
    double trans_matrix[N_PARTICLES][N_PARTICLES];
    double q[N_PARTICLES];
    double Q[N_PARTICLES];
    
    double v[N_PARTICLES];
    double m[N_PARTICLES];
    double kappa = 1, alpha = 1;

    set_initial_condition((double*) v, (double*) q, (double*) m);
    
    int end_time = 250; double dt =  1e-1;
    int n_timesteps = end_time / dt;

    velocity_verlet((int) n_timesteps, (double) dt, (double*) v, (double*) q, \
    (double*) m, (double) kappa, (double) alpha);

    construct_transformation_matrix((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES);

    // Evolove system in time 

    transform_to_normal_modes((double (*)[N_PARTICLES]) trans_matrix, (int) N_PARTICLES, (double*) q, (double*) Q);

    return 0;
}
