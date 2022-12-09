#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

double
try_lattice_constants(int N_atoms, int n_rows, int n_cols){
    // Initializing values needed for task
    double E_pot = 0; int n_lattice_params = 8;
    char filename_result[] = {"try_lattice_constants.csv"};
    double lattice_params[n_lattice_params]; double lattice_param_init = 4.05; // 4.0478; // Lattice_param, denoted a0 in document. Should be 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    
    // Set print_q=true if we want to print the lattice_parameters code is looping through.
    // Lattice params is an array around lattice_param_init with lattice_param_spacing
    bool print_q = false; double lattice_param_spacing = 0.05;
    if(print_q = true){
        printf("n/2-ix, lattice_param \n");
        for(int ix = 0; ix < n_lattice_params; ix++){
            printf("%i, ", n_lattice_params/2-ix);
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*lattice_param_spacing;
            printf("%f\n", lattice_params[ix]);
        }
        printf("\n");
    } else{
        for(int ix = 0; ix < n_lattice_params; ix++){
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*lattice_param_spacing;
        }
        printf("\n");
    }
    
    double smallest_E_potperunitcell = 0; double smallest_lattice_param = 0; 
    // For loop to run over lattice parameters
    for(int ix = 0; ix < n_lattice_params; ix++){
        // Initialize matric and update values
        double pos_matrix[256][3];//double pos_matrix[4*256][3];
        double lattice_param = lattice_params[ix];
        double lattice_volume = pow(lattice_param, 3);
        
        // Retrieve fcc with already made function
        init_fcc((double (*)[3]) pos_matrix, (int) 4, (double) lattice_param); // 4 unit cells in each direction
        
        // Retrieve potential energy with ready made function
        double cell_length = 4 * lattice_param; // 4-unit cells 
        E_pot = get_energy_AL((double (*)[3]) pos_matrix, (double) cell_length, (int) N_atoms);
        
        // Scaling result with number of unit cells and initializing result vector
        double E_pot_per_unitcell = E_pot/pow(4,3);
        double result_vec[] = {ix, lattice_param, lattice_volume, E_pot, E_pot_per_unitcell};

        //Saving smallest lattice param if it results in lowest potential energy per unit cell
        if(E_pot_per_unitcell<smallest_E_potperunitcell){
            smallest_E_potperunitcell = E_pot_per_unitcell;
            smallest_lattice_param = lattice_param;
        }

        // Saving results to csv file "try_lattice_constants.csv"
        if(ix==0){
            save_vector_to_csv(result_vec, 5, filename_result, true); // true -> fopen with "w"
        } else {
            save_vector_to_csv(result_vec, 5, filename_result, false); // false -> fopen with "a"
        }

        //Printing results in terminal if print_book is set to true
        bool print_bool = true;
        if(print_bool == true) { 
            printf("Lattice_volume: %f\n", lattice_volume);
            printf("E_pot: %f\n", E_pot);
            printf("E_pot_perunitcell: %f\n", E_pot_per_unitcell);
            printf("\n");
        }
    }
    
    return smallest_lattice_param;
}

gsl_rng *
init_random_num_generator()
{
    int seed = 42;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed); 
    return r;
}


void
displace_fcc(double positions[][3], int N, double lattice_param)
{
    int i, j, k;
    int xor_value;
    double abs_displacement = 0.065; 

    //Initializing random number generator
    gsl_rng * r;
    r = init_random_num_generator();
    
    for (i = 0; i < 2 * N; i++){
        for (j = 0; j < 2 * N; j++){
            for (k = 0; k < N; k++){
                for(int lx = 0; lx < 3; lx++){
                    double random_number = gsl_ran_flat(r,-abs_displacement, abs_displacement);
                    positions[i * N * 2 * N + j * N + k][lx] += random_number*lattice_param;
                }
            }
        }
    };
    gsl_rng_free(r);
}

double
velocity_verlet(double positions[][3], double v[][3], double lattice_param, double cell_length, int end_time, \
double dt, int n_cols, int nbr_atoms, bool temp_scaling, bool press_scaling, double temp_eq, double press_eq, bool write_not_append)
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
    double tau_T, tau_P;
    
    //char *filename_result, *filename_pos, *filename_param;
    char filename_result[] = {"vel_verlet_eq.csv"};
    char filename_pos[] = {"position_track_eq.csv"};
    char filename_param[] = {"parameters_eq.csv"};

    char filename_result_prod[] = {"vel_verlet_prod.csv"};
    char filename_pos_prod[] = {"position_track_prod.csv"};
    char filename_param_prod[] = {"parameters_prod.csv"};
    
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
                positions[ix][jx] *= pow(alpha_P, (double) 1/3);
            }
        }

        // After scaling positions lattice_parameter is scaled to change pressure
        cell_length *= pow(alpha_P, (double) 1/3); cell_volume = pow(cell_length, 3);
        

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
        E_potential_per_unitcell = E_potential / pow(4,3);
        E_kinetic_per_unitcell = E_kinetic / pow(4,3);
        E_total_per_unitcell = E_potential_per_unitcell + E_kinetic_per_unitcell;

        // Calculate virial to use for pressure calculation
        virial = get_virial_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        virial_per_unitcell = virial / pow(4,3);

        // Calculate temperature (4 atoms in unit cell, kB in eV/K)
        temp_inst_per_unitcell = 2 / (3 * 4 * kB) * E_kinetic_per_unitcell;

        // Change scaling parameter for temperature
        if(temp_scaling == true){

            // Value of tau_T is uncertain, and unsure if we should mult. by dt. Have heard approx 50.
            tau_T = 200*dt;

            // Unsure if plus or minus. Scaling with "correct" sign seams to change temp in wrong direction
            alpha_T = 1 + 2 * dt / tau_T * (temp_eq - temp_inst_per_unitcell)/temp_inst_per_unitcell;
            //alpha_T = 1 - 2 * dt / tau_T * (temp_eq - temp_inst_per_unitcell)/temp_inst_per_unitcell;
        }

         // Calculate pressure in eV/Å^3
        press_inst_per_unitcell = ((4 * kB * temp_inst_per_unitcell) + virial_per_unitcell) / cell_volume; //lattice_volume;
        press_inst_per_unitcell *= 1.60219*1e2*1e4; // 1eV/Å^3 = 160.2 gPA
        
        // Changing scaling parameter for pressure
        if(press_scaling == true){

            // Value of tau_P is uncertain, and unsure if we should mult. by dt. Have heard slightly more than 50.
            tau_P = 25 * dt;

            // isothermal compressability for aluminium in Gpa^-1
            double kappa_T = 0.01385*1e-4; // Neg eller pos?

            // Unsure if plus or minus. Scaling with "correct" sign seams to change pressure in wrong direction
            alpha_P = 1 - kappa_T * dt / tau_P * (press_eq - press_inst_per_unitcell);
            //alpha_P = 1 + kappa_T * dt / tau_P * (press_eq - press_inst_per_unitcell);
        }
        
        // Creating vectors so to save results in csv files. Can be plotted with python files plot_energy.py and plot_position_track.py
        double parameter_vec[] = {end_time, dt, lattice_param, temp_scaling, press_scaling, temp_eq, press_eq, tau_T, tau_P};
        double result_vec[] = {tx*dt, cell_length, E_potential_per_unitcell, E_kinetic_per_unitcell, E_total_per_unitcell, temp_inst_per_unitcell, press_inst_per_unitcell, alpha_T, alpha_P};
        double position_track_vec[] = {tx*dt,positions[23][0], positions[23][1], positions[23][2],\
                                             positions[135][0], positions[135][1], positions[135][2],\
                                             positions[189][0], positions[189][1], positions[189][2],\
                                             temp_inst_per_unitcell};

        // Saving results to csv files
        if(temp_scaling == true || press_scaling == true)
        {   save_vector_to_csv(result_vec, 9, filename_result, write_not_append); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos, write_not_append); // false -> fopen with "a"
            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param, write_not_append);
            }
            printf("Calibration: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        } else {
            save_vector_to_csv(result_vec, 9, filename_result_prod, write_not_append); // true -> fopen with "w"
            save_vector_to_csv(position_track_vec, 10, filename_pos_prod, write_not_append); // false -> fopen with "a"

            if(tx == n_timesteps){
                save_vector_to_csv(parameter_vec, 9, filename_param_prod, write_not_append);
            }
            printf("Production: Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
        }
        // Printing temperature for each timestep to keep track during longer measurements
        //printf("Inst. Temp, Press at t = [%i]: %f,   %f\n", tx, temp_inst_per_unitcell, press_inst_per_unitcell);
    }
    return cell_length;
}
/*
void 
H1_task1()
{
     int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;
    
    // Testing lattice parameters
    bool get_small_lattice_param = true;
    if(get_small_lattice_param == true){
    double smallest_lattice_param = 0;
    smallest_lattice_param = try_lattice_constants((int) nbr_atoms, (int) n_rows, (int) n_cols);
    printf("%f\n", smallest_lattice_param);
    }
}

void 
H1_task2()
{
     int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;

    // Initialising position and velocity arrrays
    double position[nbr_atoms][n_cols];
    double velocity[nbr_atoms][n_cols];
    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            velocity[ix][jx] = 0;
        }
    }

    // Choosing lattice param
    double lattice_param = 4.05; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling;
    double temp_eq, press_eq;

    // Production run
    end_time = 10; dt = 1e-2;
    temp_scaling = false, press_scaling = false;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq);
}

void 
H1_task3()
{
    int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;

    // Initialising position and velocity arrrays
    double position[nbr_atoms][n_cols];
    double velocity[nbr_atoms][n_cols];
    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            velocity[ix][jx] = 0;
        }
    }

    // Choosing lattice param
    double lattice_param = 4.05; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling;
    double temp_eq, press_eq;


    // Equalibration run
    end_time = 10; dt = 1e-2;
    temp_scaling = true; press_scaling = true;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq);


    // Production run
    end_time = 10; dt = 1e-2;
    temp_scaling = false, press_scaling = false;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq);

}
*/
void 
H1_task4()
{
    int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;

    // Initialising position and velocity arrrays
    double position[nbr_atoms][n_cols];
    double velocity[nbr_atoms][n_cols];
    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            velocity[ix][jx] = 0;
        }
    }

    // Choosing lattice param
    double lattice_param = 4.05; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling, write_not_append;
    double temp_eq, press_eq;


    // Melting run
    end_time = 10; dt = 1e-2;
    temp_scaling = true; press_scaling = true;
    temp_eq = 3000; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append);

    // Cooling run
    end_time = 10; dt = 1e-2;
    temp_scaling = true; press_scaling = true;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = false;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append);


    // Production run
    end_time = 10; dt = 1e-2;
    temp_scaling = false, press_scaling = false;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = false;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append);

}

int
run(
    int argc,
    char *argv[]
   )
{
    // H1_task1();
    //H1_task2();
    //H1_task3();
    //H1_task3();
    H1_task4();

    return 0;
}
