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
                    //printf("step: %i\n", i * N * 2 * N + j * N + k);
                    //printf("position value b4: %f\n", positions[i * N * 2 * N + j * N + k][lx]);
                    
                    positions[i * N * 2 * N + j * N + k][lx] += random_number*lattice_param;
                    //printf("position value after: %f\n", positions[i * N * 2 * N + j * N + k][lx]);
                    //printf("random number: %f\n", random_number*lattice_param);

                }
            }
        }
    };
    gsl_rng_free(r);
}

void velocity_verlet(int n_timesteps, int nbr_atoms, int n_cols, double q[][3], double v[][3], double dt, double m, double lattice_param)
{
    double cell_length = 4 * lattice_param;
    
    //double **f = create_2D_array((int) nbr_atoms, n_cols);

    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            //f[ix][jx] = 0;
            //printf("v[%i][%i]: %f\n", ix, jx, v[ix][jx]);
        }
    }
    
   
    for(int tx = 1; tx < 1+1; tx++)
    {   
        // OBS!!!! Should get forces be row below or outside of tx for-loop?
        // Calculate forces based on atom positions
        double f[nbr_atoms][n_cols];
        printf("v %f\n", v[213][0]);
        get_forces_AL((double (*)[3]) f, (double (*)[3]) q, (double) cell_length, (int) nbr_atoms);
        /* v(t+dt/2) */
        printf("tx: %i\n", tx);
        printf("v %f\n", v[213][0]);
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                //printf("f[%i][%i] %f\n",ix, jx, f[ix][jx]);
                //printf("b4 v[%i][%i] %f\n", ix, jx, v[ix][jx]);
                v[ix][jx] += 0.5 * dt * f[ix][jx] / m;
                //printf("aft v[%i][%i] %f\n", ix, jx, v[ix][jx]);
                //printf("f[%i][%i] %f\n",ix, jx, f[ix][jx]);
                //printf("aft v[%i][%i] %f\n",ix, jx, v[ix][jx]);
            }
        }

        /* q(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                //q[ix][jx] += dt * v[ix][jx];
                //printf("v[%i][%i] %f\n",ix, jx, v[ix][jx]);
                //f[ix][jx] = 0;
            }
        }
        
        // Calculate forces based on (new) atom positions
        //get_forces_AL((double (*)[3]) f, (double (*)[3]) q, (double) cell_length, (int) nbr_atoms);

        /* v(t+dt) */
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                //v[ix][jx] += 0.5 * dt * f[ix][jx] / m;
                //printf("Hej");
            }
        }
        printf("Hejsan");
        // save positions - get interesting values for timestep / periodix boundary cond?
        // get potenital kinetic and total energies
        // Get temperature and pressure
    }
    //destroy_2D_array(f, nbr_atoms);
}

void
sol_eq_of_motion(double positions[][3], double lattice_param, int n_rows, int n_cols, int nbr_atoms)
{
    // Initialize variables
    double cell_length = 4 * lattice_param; double lattice_volume = pow(lattice_param, 3);
    double aluminium_amu = 26.98; double m_asu = 9649;
    double aluminium_asu = aluminium_amu/m_asu;

    double E_kinetic = 0; double E_kinetic_per_unitcell = 0;
    double E_potential = 0; double E_potential_per_unitcell = 0;
    double E_total_per_unitcell = 0;

    double virial = 0; double virial_per_unitcell = 0;
    double temp_inst = 0; double temp_inst_per_unitcell = 0;
    double press_inst = 0; double press_inst_per_unitcell = 0;
    double kB = 8.61733 * 1e-5; // Boltzzmann constant in eV/K
    
    char filename_result[] = {"eq_of_motion.csv"};
    
    //double q[nbr_atoms][n_cols];
    double v[nbr_atoms][n_cols];
    double f[nbr_atoms][n_cols];

    
    for(int ix = 0; ix < nbr_atoms; ix++){
        for(int jx = 0; jx < n_cols; jx++){
            //q[ix][jx] = positions[ix][jx];
            v[ix][jx] = 0;
            f[ix][jx] = 0;
        }
    }
    double end_time = 10; double dt = 1e-2;
    int n_timesteps = end_time / dt; 

    double alpha_T = 1; double alpha_P = 1;
    double temp_eq = 500; double press_eq = 1*1e-4; // 1 bar = 1e-4 GPa
    bool temp_scaling = true; bool press_scaling = true;

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
                positions[ix][jx] += dt * v[ix][jx] * pow(alpha_P, 1/3);
                //printf("v[%i][%i] %f\n",ix, jx, v[ix][jx]);
            }
        }

        /* a(t+dt) */
        get_forces_AL((double (*)[3]) f, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);

        /* v(t+dt) */
        E_kinetic = 0;
        for(int ix = 0; ix < nbr_atoms; ix++)
        {
            for(int jx = 0; jx < n_cols; jx++)
            {
                v[ix][jx] += 0.5 * dt * f[ix][jx] / aluminium_asu * sqrt(alpha_T); // alpha_T for scaling. If temp_scaling==false this is 1.
                E_kinetic += 0.5 * aluminium_asu * v[ix][jx] * v[ix][jx];
                //printf("v[%i][%i] %f\n",ix, jx, v[ix][jx]);
                
            }
        }
        // save positions - get interesting values for timestep / periodix boundary cond?
        // get potenital kinetic and total energies
        
        // Get temperature and pressure
        
        E_potential = get_energy_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        E_potential_per_unitcell = E_potential / pow(4,3);
        E_kinetic_per_unitcell = E_kinetic / pow(4,3);
        E_total_per_unitcell = E_potential_per_unitcell + E_kinetic_per_unitcell;

        virial = get_virial_AL((double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
        virial_per_unitcell = virial / pow(4,3);

        temp_inst_per_unitcell = 2 / (3 * 4 * kB) * E_kinetic_per_unitcell; // 4 atoms in an Fcc unit cell and kB in eV/K
        if(temp_scaling == true){
            double tau_T = 50; // What should this be // Åtminstone några hundra tidssteg, typ 50/Viktor
            alpha_T = 1 + 2*dt/tau_T * (temp_eq - temp_inst_per_unitcell)/temp_inst_per_unitcell;
        }
         // ev/Å^3 to SI gger konstant för asu to SI, Pascale är Joule per kubikmeter
        press_inst_per_unitcell = ((4 * kB * temp_inst_per_unitcell) + virial_per_unitcell) / lattice_volume;
        if(press_scaling == true){
            double tau_P = 50; // What should this be // Åtminstone några hundra tidssteg, typ något mer än 100 /Viktor
            double kappa_T = 0.01385; // isothermal compressability for aluminium in Gpa^-1
            alpha_P = 1 - kappa_T * dt / tau_P * (press_eq - press_inst_per_unitcell);
        }
        
        double result_vec[] = {tx*dt, lattice_param, E_potential_per_unitcell, E_kinetic_per_unitcell, E_total_per_unitcell, temp_inst_per_unitcell, press_inst_per_unitcell, alpha_T, alpha_P};

        // Saving results to csv file "try_lattice_constants.csv"
        if(tx == 1){
            save_vector_to_csv(result_vec, 9, filename_result, true); // true -> fopen with "w"
        } else {
            save_vector_to_csv(result_vec, 9, filename_result, false); // false -> fopen with "a"
        }

        printf("Inst. Temp at t = [%i]: %f\n", tx, temp_inst_per_unitcell);
    }
}

int
run(
    int argc,
    char *argv[]
   )
{
    int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;
    
    bool get_small_lattice_param = false;
    if(get_small_lattice_param == true){
    double smallest_lattice_param = 0;
    smallest_lattice_param = try_lattice_constants((int) nbr_atoms, (int) n_rows, (int) n_cols);
    printf("%f\n", smallest_lattice_param);
    }

    
    double lattice_param = 4.05; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)

    double position[nbr_atoms][n_cols];
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    sol_eq_of_motion((double (*)[3]) position, (double) lattice_param, (int) n_rows, (int) n_cols, (int) nbr_atoms);
    
    return 0;
}
