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

void calc_acc_carb(double *a, double *u, double *m, double kappa, int size_of_u)
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

void calc_acc_co2(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    //a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    //a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    a[0] = kappa*(u[1] - u[0])/m[0];
    a[size_of_u - 1] = - kappa*(u[size_of_u - 1]-u[size_of_u - 2])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        //a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
        //a[i] = - kappa*(u[1] - u[0])/m[0] + kappa*(u[size_of_u - 1]-u[size_of_u - 2])/m[size_of_u - 2];
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
void velocity_verlet_carb(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa, double **v_matrix)
{

    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];

    calc_acc_carb(a, q, m, kappa, n_particles);
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
        calc_acc_carb(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
            v_matrix[i][j] = v[j];
        }
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

    }
}

void velocity_verlet_co2(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa, double **v_matrix)
{
    printf("size of u -1: %i\n", n_particles-1);
    printf("m[0]= %f: \n", m[0]);
    printf("m[n_particles-1]= %f: \n", m[n_particles-1]);
    printf("m[n_particles-2]= %f: \n", m[n_particles-2]);
    
    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];

    calc_acc_co2(a, q, m, kappa, n_particles);
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
        calc_acc_co2(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
            v_matrix[i][j] = v[j];
        }
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

    }
}

void carbon()
{
    // Setting standard variables
    //int t_meas_duration = 0.25; // (ps)7
    //convertions factor for N/m to asu is 1/16.0218
    double dt = 1e-5; // /2;
    double force_conv = 1.0/16.0218;
    int n_timesteps = 250*100; //4;
    int n_particles = 3; double kappa = 1000 * force_conv;
    printf("kappa: %f \n", kappa);
    
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
    q_1[0] = 0.01; q_2[0] = 0.005; q_3[0] = -0.005; v[0] = 0;

    //Creating matrix for velocity
    double **v_matrix = NULL;
    create_2D_array(&v_matrix, n_timesteps+1, n_particles);

    // Performing velocity verlet algorithm
    velocity_verlet_carb((int) n_timesteps, (int) n_particles, (double *) v, \
                    (double *) q_1, (double *) q_2, (double *) q_3, \
                    (double) dt, (double *) m, (double) kappa, (double **) v_matrix);

    // Saving trajectory vectors to file
    double **q_matrix = NULL;

    // Create_2D_array(&q_matrix, n_particles, n_timesteps);
    create_2D_array(&q_matrix, n_timesteps, n_particles);
    for(int ix = 0; ix < n_timesteps; ix++)
    {
        q_matrix[ix][0] = q_1[ix];
        q_matrix[ix][1] = q_2[ix];
        q_matrix[ix][2] = q_3[ix];
    }
    
    char filename_position[] = {"Saved_trajectory_carb_small_dt.csv"};
    save_matrix_to_csv(q_matrix, n_particles, n_timesteps, filename_position);

    char filename_position_p1[] = {"Saved_trajectory_particle1_carb_small_dt.csv"};
    save_vector_to_csv(q_1, n_timesteps, filename_position_p1);
    char filename_position_p2[] = {"Saved_trajectory_particle2_carb_small_dt.csv"};
    save_vector_to_csv(q_2, n_timesteps, filename_position_p2);
    char filename_position_p3[] = {"Saved_trajectory_particle3_carb_small_dt.csv"};
    save_vector_to_csv(q_3, n_timesteps, filename_position_p3);

    //print_vector(v, n_particles);
    
    char filename_velocity[] = {"Saved_velocity_carb_small_dt.csv"};
    save_matrix_to_csv(v_matrix, n_particles, n_timesteps+1, filename_velocity);
    // save_matrix_to_csv(v_matrix, n_particles, n_timesteps-1, filename_velocity);
    
    destroy_2D_array(q_matrix); destroy_2D_array(v_matrix);
}

void carbon_dioxide()
{
    // Setting standard variables
    //int t_meas_duration = 0.25; // (ps)7
    //convertions factor for N/m to asu is 1/16.0218
    double dt = 1e-5; // /2;
    double force_conv = 1.0/16.0218;
    int n_timesteps = 250*100; //4;
    int n_particles = 3; double kappa = 1600 * force_conv;
    printf("kappa: %f \n", kappa);
    
    // Retrieving mass vector with carbon mass in asu with sizeof(n_particles)
    double carbon_amu = 15.9994; double m_asu = 9649; // m_asu=eV*(ps)^2/Å^2 
    double carbon_asu = carbon_amu/m_asu;
    double *m = calloc(sizeof(double), n_particles);

    double oxy_amu = 12.01; // m_asu=eV*(ps)^2/Å^2 
    double oxy_asu = oxy_amu/m_asu;
    
    for(int ix = 0; ix < n_particles; ix++)
    {
        m[ix] = oxy_asu;
        if(ix == 1)
        {
            m[ix] = carbon_asu;
        }
    }
    print_vector(m, n_particles);

    // (Empty allocated) velocity array: sizeof(v) = n_particles
    double *v = calloc(sizeof(double), n_particles);
    

    // (Empty allocated) q arrays: sizeof(q) = n_timesteps +1
    double *q_1 = calloc(sizeof(double), n_timesteps+1);
    double *q_2 = calloc(sizeof(double), n_timesteps+1);
    double *q_3 = calloc(sizeof(double), n_timesteps+1);

    // Setting initial conditions (Å)
    q_1[0] = 0.01; q_2[0] = 0.005; q_3[0] = -0.005; v[0] = 0;

    //Creating matrix for velocity
    double **v_matrix = NULL;
    create_2D_array(&v_matrix, n_timesteps+1, n_particles);

    // Performing velocity verlet algorithm
    velocity_verlet_co2((int) n_timesteps, (int) n_particles, (double *) v, \
                    (double *) q_1, (double *) q_2, (double *) q_3, \
                    (double) dt, (double *) m, (double) kappa, (double **) v_matrix);

    // Saving trajectory vectors to file
    double **q_matrix = NULL;

    // Create_2D_array(&q_matrix, n_particles, n_timesteps);
    create_2D_array(&q_matrix, n_timesteps, n_particles);
    for(int ix = 0; ix < n_timesteps; ix++)
    {
        q_matrix[ix][0] = q_1[ix];
        q_matrix[ix][1] = q_2[ix];
        q_matrix[ix][2] = q_3[ix];
    }
    
    char filename_position[] = {"Saved_trajectory_co2_small_dt.csv"};
    save_matrix_to_csv(q_matrix, n_particles, n_timesteps, filename_position);

    char filename_position_p1[] = {"Saved_trajectory_particle1_co2_small_dt.csv"};
    save_vector_to_csv(q_1, n_timesteps, filename_position_p1);
    char filename_position_p2[] = {"Saved_trajectory_particle2_co2_small_dt.csv"};
    save_vector_to_csv(q_2, n_timesteps, filename_position_p2);
    char filename_position_p3[] = {"Saved_trajectory_particle3_co2_small_dt.csv"};
    save_vector_to_csv(q_3, n_timesteps, filename_position_p3);

    //print_vector(v, n_particles);
    
    char filename_velocity[] = {"Saved_velocity_co2_small_dt.csv"};
    save_matrix_to_csv(v_matrix, n_particles, n_timesteps+1, filename_velocity);
    // save_matrix_to_csv(v_matrix, n_particles, n_timesteps-1, filename_velocity);
    
    destroy_2D_array(q_matrix); destroy_2D_array(v_matrix);
}

int main(int argc, char **argv)
{
    carbon();
    carbon_dioxide();
    return 0;
}