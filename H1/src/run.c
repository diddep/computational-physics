#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

void
try_lattice_constants(){
    double **pos_matrix = NULL;
    int N_atoms = 256; int n_rows = 4 * pow(N_atoms, 3); int n_cols = 3;
    
    
    double lattice_param = 4; // 4.0478; // Lattice_param, denoted a0 in document. Should be 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double lattice_volume = pow(lattice_param, 3);
    create_2D_array(&pos_matrix, n_rows, n_cols); // Has to be size [4*N*N*N][3]

    init_fcc((double (*)[3]) pos_matrix, (int) N_atoms, (double) lattice_param);
    
    double cell_length = n_rows * lattice_param;
    double E_pot = 0;
    E_pot = get_energy_AL((double (*)[3]) pos_matrix, (double) cell_length, (int) N_atoms);
    
    double E_pot_per_unitcell = E_pot/lattice_param;


    printf("Lattice_volume: %f\n", lattice_volume);
    printf("E_pot: %f\n", E_pot);
    printf("E_pot_perunitcell: %f\n", E_pot_per_unitcell);

    destroy_2D_array(pos_matrix);
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
    
    try_lattice_constants();

    return 0;
}
