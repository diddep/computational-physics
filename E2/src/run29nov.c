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

void velocity_verlet(int n_timesteps, double dt, double **v, double **q, \
 double *m, double kappa, double alpha)
{
    double m = 1;
    double kappa = 1, alpha = 1
    
    

#define kappa 1
}

int run()
{
    
    double trans_matrix[N_PARTICLES][N_PARTICLES];
    double q[N_PARTICLES];
    double Q[N_PARTICLES];

    construct_transformation_matrix(trans_matrix, N_PARTICLES);

    // Evolove system in time 

    transform_to_normal_modes(trans_matrix, N_PARTICLES, q, Q);

    return 0;
}
