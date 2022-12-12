
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "try_lattice_constants.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

double
velocity_verlet(double positions[][3], double v[][3], double lattice_param, double cell_length, int end_time, \
double dt, int n_cols, int nbr_atoms, bool temp_scaling, bool press_scaling, double temp_eq, double press_eq, bool write_not_append, \
double tau_P, double tau_T)
{
    // Initialize variables
    double cell_volume = pow(cell_length, 3.);
    double lattice_volume = pow(lattice_param, 3.);
    double aluminium_amu = 26.98; double m_asu = 9649;
    double aluminium_asu = aluminium_amu/m_asu; //Mass in atomic simulation units

    double E_kinetic = 0; double E_kinetic_per_unitcell = 0;
    double E_potential = 0; double E_potential_per_unitcell = 0;
    double E_total_per_unitcell = 0;

    double virial = 0; double virial_per_unitcell = 0;
    double temp_inst = 0; double temp_inst_per_unitcell = 0;
    double press_inst = 0; double press_inst_per_unitcell = 0;
    double kB = 8.61733 * 1e-5; // Boltzmann constant in eV/K
    
    //char *filename_result, *filename_pos, *filename_param;
    char filename_result[] = {"../csv/vel_verlet_eq.csv"};
    char filename_pos[] = {"../csv/position_track_eq.csv"};
    char filename_param[] = {"../csv/parameters_eq.csv"};

    char filename_result_prod[] = {"../csv/vel_verlet_prod.csv"};
    char filename_pos_prod[] = {"../csv/position_track_prod.csv"};
    char filename_param_prod[] = {"../csv/parameters_prod.csv"};
    char filename_radial_distribution[] = {"../csv/radial_distribution.csv"};
    bool is_write;
    
    //Creating empty arrays
    //double v[nbr_atoms][n_cols];
    double f[nbr_atoms][n_cols];

    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            //v[ix][jx] = 0;
            f[ix][jx] = 0;
        }
    }
    
    // Variables for duration of measurement
    int n_timesteps = end_time / dt; 

    // Declaring scaling variables, if temp/press_scaling = false scaling is turned off and alpha_T/P just remains 1
    double alpha_T = 1; double alpha_P = 1; 

    //Velocity verlet algorithm as in E1
    get_forces_AL((double (*)[3]) f, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);

    
    for(int tx = 1; tx < n_timesteps + 1; tx++)
    {
        /* v(t+dt/2) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                v[ix][jx] += 0.5 * dt * f[ix][jx] / aluminium_asu;
                //printf("v[%i][%i] %f\n",ix, jx, v[ix][jx]);
            }
        }
        
        /* q(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                // alpha_P for scaling. If press_scaling==false this is 1.
                positions[ix][jx] += dt * v[ix][jx];
                positions[ix][jx] *= cbrt(alpha_P);
                //positions[ix][jx] *= pow(alpha_P, (double) 1./3.);
            }
        }

        // After scaling positions lattice_parameter is scaled to change pressure
        cell_length *= cbrt(alpha_P); cell_volume = cell_length * cell_length * cell_length;
        

        /* a(t+dt) */
        get_forces_AL((double (*)[3]) f, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);

        // Resetting kinetic energy for every timestep to ensure summing over particles and not over timesteps
        E_kinetic = 0;

        /* v(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {   
                // alpha_T for scaling. If temp_scaling==false this is 1.
                v[ix][jx] += 0.5 * dt * f[ix][jx] / aluminium_asu;
                v[ix][jx] *= sqrt(alpha_T);
                E_kinetic += 0.5 * aluminium_asu * v[ix][jx] * v[ix][jx];
            }
        }
        
        // Calculate potential, kinetic and total energy per unitcell
        E_potential = get_energy_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        E_potential_per_unitcell = E_potential / pow(4.,3.);
        E_kinetic_per_unitcell = E_kinetic / pow(4.,3.);
        E_total_per_unitcell = E_potential_per_unitcell + E_kinetic_per_unitcell;

        // Calculate virial to use for pressure calculation
        virial = get_virial_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        virial_per_unitcell = virial / pow(4.,3.);

        double N = 256;
        double N1 = 4;
        temp_inst_per_unitcell = 2 / (3 * N * kB) * E_kinetic;
        //temp_inst_per_unitcell = 2 / (3 * N1 * kB) * E_kinetic_per_unitcell;

        // Change scaling parameter for temperature
        if(temp_scaling == true){

            alpha_T = 1 + 2 * dt / tau_T * (temp_eq - temp_inst_per_unitcell)/temp_inst_per_unitcell;


         // Calculate pressure in eV/Å^3
        press_inst_per_unitcell = ((N * kB * temp_inst_per_unitcell) + virial) / cell_volume; //lattice_volume;
        
        //press_inst_per_unitcell = ((N * kB * temp_inst_per_unitcell) + virial_per_unitcell) / cell_volume; //lattice_volume;

        press_inst_per_unitcell *= 1.60219*1e2*1e4; // 1eV/Å^3 = 160.2 Gpa, 1 GPa = 10^4 Bar 
        
        // Changing scaling parameter for pressure
        if(press_scaling == true){

            double kappa_T = 0.01385*1e-4; // Neg eller pos?

            // Unsure if plus or minus. Scaling with "correct" sign seams to change pressure in wrong direction
            alpha_P = 1 - kappa_T * dt / tau_P * (press_eq - press_inst_per_unitcell);
        }
        
        // Creating vectors so to save results in csv files. Can be plotted with python files plot_energy.py and plot_position_track.py
        double parameter_vec[] = {end_time, dt, lattice_param, temp_scaling, press_scaling, temp_eq, press_eq, tau_T, tau_P};
        double result_vec[] = {tx*dt, cell_length, E_potential_per_unitcell, E_kinetic_per_unitcell, E_total_per_unitcell, temp_inst_per_unitcell, press_inst_per_unitcell, alpha_T, alpha_P};
        double position_track_vec[] = {tx*dt,positions[23][0], positions[23][1], positions[23][2],\
                                             positions[135][0], positions[135][1], positions[135][2],\
                                             positions[189][0], positions[189][1], positions[189][2],\
                                             temp_inst_per_unitcell};

        // Saving results to csv files
        if(tx == 1){is_write = write_not_append;} else {is_write = false;};
        if(temp_scaling == true || press_scaling == true)
        {   save_vector_to_csv(result_vec, 9, filename_result, is_write); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos, is_write); // false -> fopen with "a"
            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param, is_write);
            }
            printf("Calibration: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        } else {
            save_vector_to_csv(result_vec, 9, filename_result_prod, is_write); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos_prod, is_write); // false -> fopen with "a"

            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param_prod, is_write);
            }
            printf("Production: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        }
        // Printing temperature for each timestep to keep track during longer measurements
        //printf("Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
    }
    return cell_length;
}


//

double
velocity_verlet_with_radial(double positions[][3], double v[][3], double lattice_param, double cell_length, int end_time, \
double dt, int n_cols, int nbr_atoms, bool temp_scaling, bool press_scaling, double temp_eq, double press_eq, bool write_not_append, \
double tau_P, double tau_T, double *radial_histogram_vector, int number_of_bins)
{
    // Initialize variables
    double cell_volume = pow(cell_length, 3);
    double lattice_volume = pow(lattice_param, 3);
    double aluminium_amu = 26.98; double m_asu = 9649;
    double aluminium_asu = aluminium_amu/m_asu; //Mass in atomic simulation units

    double E_kinetic = 0; double E_kinetic_per_unitcell = 0;
    double E_potential = 0; double E_potential_per_unitcell = 0;
    double E_total_per_unitcell = 0;

    double virial = 0; double virial_per_unitcell = 0;
    double temp_inst = 0; double temp_inst_per_unitcell = 0;
    double press_inst = 0; double press_inst_per_unitcell = 0;
    double kB = 8.61733 * 1e-5; // Boltzmann constant in eV/K
    
    //char *filename_result, *filename_pos, *filename_param;
    char filename_result[] = {"../csv/vel_verlet_eq.csv"};
    char filename_pos[] = {"../csv/position_track_eq.csv"};
    char filename_param[] = {"../csv/parameters_eq.csv"};

    char filename_result_prod[] = {"../csv/vel_verlet_prod.csv"};
    char filename_pos_prod[] = {"../csv/position_track_prod.csv"};
    char filename_param_prod[] = {"../csv/parameters_prod.csv"};
    char filename_radial_distribution[] = {"../csv/radial_distribution.csv"};
    bool is_write;
    
    //Creating empty arrays
    //double v[nbr_atoms][n_cols];
    double f[nbr_atoms][n_cols];

    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            //v[ix][jx] = 0;
            f[ix][jx] = 0;
        }
    }
    
    // Variables for duration of measurement
    int n_timesteps = end_time / dt; 

    // Declaring scaling variables, if temp/press_scaling = false scaling is turned off and alpha_T/P just remains 1
    double alpha_T = 1; double alpha_P = 1; 

    //Velocity verlet algorithm as in E1
    get_forces_AL((double (*)[3]) f, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);

    
    get_radial_dist_AL(number_of_bins, radial_histogram_vector, (double (*)[3])positions, cell_length, nbr_atoms);
    
    for(int tx = 1; tx < n_timesteps + 1; tx++)
    {
        /* v(t+dt/2) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                v[ix][jx] += 0.5 * dt * f[ix][jx] / aluminium_asu;
                //printf("v[%i][%i] %f\n",ix, jx, v[ix][jx]);
            }
        }
        
        /* q(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                // alpha_P for scaling. If press_scaling==false this is 1.
                positions[ix][jx] += dt * v[ix][jx];
                positions[ix][jx] *= cbrt(alpha_P);
                //positions[ix][jx] *= pow(alpha_P, (double) 1/3);
            }
        }

        // After scaling positions lattice_parameter is scaled to change pressure
        cell_length *= cbrt(alpha_P); cell_volume = cell_length*cell_length*cell_length;
        

        /* a(t+dt) */
        get_forces_AL((double (*)[3]) f, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        get_radial_dist_AL(number_of_bins, radial_histogram_vector, (double (*)[3])positions, cell_length, nbr_atoms);
    
        // Resetting kinetic energy for every timestep to ensure summing over particles and not over timesteps
        E_kinetic = 0;

        /* v(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {   
                // alpha_T for scaling. If temp_scaling==false this is 1.
                v[ix][jx] += 0.5 * dt * f[ix][jx] / aluminium_asu;
                v[ix][jx] *= sqrt(alpha_T);
                E_kinetic += 0.5 * aluminium_asu * v[ix][jx] * v[ix][jx];
            }
        }
        
        // Calculate potential, kinetic and total energy per unitcell
        E_potential = get_energy_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        E_potential_per_unitcell = E_potential / 64;// pow(4.,3.);
        E_kinetic_per_unitcell = E_kinetic / 64; //pow(4.,3.);
        E_total_per_unitcell = E_potential_per_unitcell + E_kinetic_per_unitcell;

        // Calculate virial to use for pressure calculation
        virial = get_virial_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        virial_per_unitcell = virial / 64;//pow(4.,3.);
    
        double N = 256;

        // Calculate temperature (4 atoms in unit cell, kB in eV/K)
        temp_inst_per_unitcell = 2 / (3 * N * kB) * E_kinetic;
        //temp_inst_per_unitcell = 2 / (3 * N1 * kB) * E_kinetic_per_unitcell;

        // Change scaling parameter for temperature
        if(temp_scaling == true){

            alpha_T = 1 + 2 * dt / tau_T * (temp_eq - temp_inst_per_unitcell)/temp_inst_per_unitcell;
        }
        press_inst_per_unitcell = ((N * kB * temp_inst_per_unitcell) + virial) / cell_volume; //lattice_volume;
        press_inst_per_unitcell *= 1.60219*1e2*1e4; // 1eV/Å^3 = 160.2 gPA
        
        // Changing scaling parameter for pressure
        if(press_scaling == true){

            // Isothermal compressability for aluminium in Gpa^-1
            double kappa_T = 0.01385*1e-4;

            alpha_P = 1 - kappa_T * dt / tau_P * (press_eq - press_inst_per_unitcell);
        }
        
        // Creating vectors so to save results in csv files. Can be plotted with python files plot_energy.py and plot_position_track.py
        double parameter_vec[] = {end_time, dt, lattice_param, temp_scaling, press_scaling, temp_eq, press_eq, tau_T, tau_P};
        double result_vec[] = {tx*dt, cell_length, E_potential_per_unitcell, E_kinetic_per_unitcell, E_total_per_unitcell, temp_inst_per_unitcell, press_inst_per_unitcell, alpha_T, alpha_P};
        double position_track_vec[] = {tx*dt,positions[23][0], positions[23][1], positions[23][2],\
                                             positions[135][0], positions[135][1], positions[135][2],\
                                             positions[189][0], positions[189][1], positions[189][2],\
                                             temp_inst_per_unitcell};

        // Saving results to csv files
        if(tx == 1){is_write = write_not_append;} else {is_write = false;};
        if(temp_scaling == true || press_scaling == true)
        {   save_vector_to_csv(result_vec, 9, filename_result, is_write); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos, is_write); // false -> fopen with "a"
            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param, is_write);
            }
            printf("Calibration: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        } else {
            save_vector_to_csv(result_vec, 9, filename_result_prod, is_write); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos_prod, is_write); // false -> fopen with "a"

            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param_prod, is_write);
            }
            printf("Production: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        }
       }
    return cell_length;
}

/* Returns the forces */
void get_radial_dist_AL(int number_of_bins, double *radial_histogram_vector, double positions[][3], double cell_length, int nbr_atoms)
{
  int i, j;
  double cell_length_inv, cell_length_sq, bin_length;
  double rcut, rcut_sq;
  //double densityi, dens, drho_dr, force;
  //double dUpair_dr;
  double sxi, syi, szi, sxij, syij, szij, rij,  rij_sq;
  
  double *sx = malloc(nbr_atoms * sizeof (double));
  double *sy = malloc(nbr_atoms * sizeof (double));
  double *sz = malloc(nbr_atoms * sizeof (double));
  
  rcut = 6.06; // Embedded atom method potential. 
  rcut_sq = rcut * rcut;

  
  cell_length_inv = 1 / cell_length;
  cell_length_sq = cell_length * cell_length;
  bin_length= cell_length/(double)number_of_bins;
  
  for (i = 0; i < nbr_atoms; i++){
    sx[i] = positions[i][0] * cell_length_inv;
    sy[i] = positions[i][1] * cell_length_inv;
    sz[i] = positions[i][2] * cell_length_inv;
  }

  /* Compute radial distribution on atoms. */
	/* Loop over atoms again :-(. */
  
  for (i = 0; i < nbr_atoms; i++) {
    /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) 
    {
      /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
      
      /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
      
      /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
    
      /*Add position into bin depending on radial size*/
      rij = sqrt( rij_sq );
      
      int bin = (int) floor(  rij/bin_length + 0.5);

      radial_histogram_vector[bin] +=1;

    }
  }
  
  free(sx); free(sy); free(sz); sx = NULL; sy = NULL; sz = NULL;
}