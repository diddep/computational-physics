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
    bool print_q=false; double lattice_param_spacing = 0.05;
    if(print_q=true){
        printf("n/2-ix, lattice_param \n");
        for(int ix = 0; ix < n_lattice_params; ix++){
            printf("%i, ", n_lattice_params/2-ix);
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*lattice_param_spacing;
            printf("%f\n", lattice_params[ix]);
        }
        printf("\n");
    } else {
        for(int ix = 0; ix < n_lattice_params; ix++){
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*lattice_param_spacing;
        }
        printf("\n");
    }
    
    double smallest_E_potperunitcell = 0; double smallest_lattice_param = 0; 
    // For loop to run over lattice parameters
    for(int ix = 0; ix < n_lattice_params; ix++){
        // Initialize matric and update values
        double pos_matrix[4*256][3];
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
    }
    gsl_rng_free(r);
}


void
sol_eq_of_motion(double positions[][3], double lattice_param, int n_rows, int n_cols, int nbr_atoms)
{
    double cell_length = 4 * lattice_param; 
    double forces[4*256][3];
    get_forces_AL((double (*)[3]) forces, (double (*)[3]) positions, (double) cell_length, (int) nbr_atoms);
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

    double position[4*256][3];
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    sol_eq_of_motion((double (*)[3]) position, (double) lattice_param, (int) n_rows, (int) n_cols, (int) nbr_atoms);
    
    return 0;
}
