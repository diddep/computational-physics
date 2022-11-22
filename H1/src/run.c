#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"

double
try_lattice_constants(int N_atoms, int n_rows, int n_cols){
    // Initializing and creating matrix to store positions
    double **pos_matrix = NULL;  double E_pot = 0;
    char filename_result[] = {"try_lattice_constants.csv"};
    create_2D_array(&pos_matrix, n_rows, n_cols); // Has to be size [4*N*N*N][3]

    int n_lattice_params = 8;
    double lattice_params[n_lattice_params]; double lattice_param_init = 4.05; // 4.0478; // Lattice_param, denoted a0 in document. Should be 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    bool print_q=false;
    if(print_q=true){
        printf("n/2-ix, lattice_param \n");
        for(int ix = 0; ix < n_lattice_params; ix++){
            printf("%i, ", n_lattice_params/2-ix);
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*0.05;
            printf("%f\n", lattice_params[ix]);
        }
        printf("\n");
    } else {
        for(int ix = 0; ix < n_lattice_params; ix++){
            lattice_params[ix] = lattice_param_init - (n_lattice_params/2-ix)*0.05;
        }
        printf("\n");
    }
    
    double smallest_E_potperunitcell = 0; double smallest_lattice_param = 0; 
    // For loop to run over lattice parameters
    for(int ix = 0; ix < n_lattice_params; ix++){
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
    destroy_2D_array(pos_matrix);
    return smallest_lattice_param;
}


void
sol_eq_of_motion()
{
    // 1. Introduce small displacement in inital lattice positions

}
int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
    
    int N_atoms = 256; int n_rows = 4 * pow(N_atoms, 3); int n_cols = 3; // Antal unit cells in each direction
    
    bool get_small_lattice_param = false;
    if(get_small_lattice_param == true){
    double smallest_lattice_param = 0;
    smallest_lattice_param = try_lattice_constants((int) N_atoms, (int) n_rows, (int) n_cols);
    printf("%f\n", smallest_lattice_param);
    }

    
    double lattice_param = 4.05; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double **pos_matrix = NULL;
    create_2D_array(&pos_matrix, n_rows, n_cols); // Has to be size [4*N*N*N][3]
    init_fcc((double (*)[3]) pos_matrix, (int) 4, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) pos_matrix, (int) 4, (double) lattice_param);
    
    destroy_2D_array(pos_matrix);
    
    return 0;
}
