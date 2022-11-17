/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"

void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm 
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa)
{

    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];

    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

    }
    }

int main(int argc, char **argv)
{
    // Setting standard variables
    //int t_meas_duration = 0.25; // (ps)
    double dt = 1e-3;
    int n_timesteps = 250;
    int n_particles = 3; double kappa = 62.46; 
    
    // Retrieving mass vector with carbon mass in asu with sizeof(n_particles)
    double carbon_amu = 12.01; double m_asu = 9649; // m_asu=eV*(ps)^2/Å^2 
    double carbon_asu = carbon_amu/m_asu;
    double *m = calloc(sizeof(double), n_particles);
    
    for(int ix = 0; ix < n_particles; ix++)
    {
        m[ix] = carbon_asu;
    }

    // (Empty allocated) velocity array: sizeof(v) = n_particles
    double *v = calloc(sizeof(double), n_particles);
    

    // (Empty allocated) q arrays: sizeof(q) = n_timesteps +1
    double *q_1 = calloc(sizeof(double), n_timesteps+1);
    double *q_2 = calloc(sizeof(double), n_timesteps+1);
    double *q_3 = calloc(sizeof(double), n_timesteps+1);

    // Setting initial conditions (Å)
    q_1[0] = 0.01; q_2[0] = 0; q_3[0] = 0; v[0] = 0; 

    // Performing velocity verlet algorithm
    velocity_verlet((int) n_timesteps, (int) n_particles, (double *) v, \
                    (double *) q_1, (double *) q_2, (double *) q_3, \
                    (double) dt, (double *) m, (double) kappa);

    // Saving trajectory vectors to file
    double **q_matrix = NULL;
    //create_2D_array(&q_matrix, n_particles, n_timesteps);
    create_2D_array(&q_matrix, n_timesteps, n_particles);
    for(int ix = 0; ix < n_timesteps; ix++)
    {
        q_matrix[ix][0] = q_1[ix];
        q_matrix[ix][1] = q_2[ix];
        q_matrix[ix][2] = q_3[ix];
    }
    
    char filename_position[] = {"Saved_trajectory.csv"};
    save_matrix_to_csv(q_matrix, n_particles, n_timesteps, filename_position);

    double **v_matrix = NULL;
    create_2D_array(&v_matrix, n_timesteps - 1, n_particles);

    for(int ix = 0; ix < n_timesteps - 1; ix++)
    {
        for(int jx = 0; jx < n_particles; jx++)
        {
            v_matrix[ix][jx] = (q_matrix[ix+1][jx]-q_matrix[ix][jx])/dt;
        }
    }
    print_vector(v, n_particles);
    char filename_velocity[] = {"Saved_velocity.csv"};
    save_matrix_to_csv(v_matrix, n_particles, n_timesteps-1, filename_velocity);
    
    destroy_2D_array(q_matrix); destroy_2D_array(v_matrix);

    return 0;
}