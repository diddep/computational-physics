#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "try_lattice_constants.h"
#include "velocity_verlet.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

void H1_task1(), H1_task2(), H1_task3(), H1_task4(), H1_task6();

int
run(
    int argc,
    char *argv[]
   )
{
    //H1_task1();
    //H1_task2();
    //H1_task3();
    //H1_task4();
    H1_task6();

    return 0;
}

void 
H1_task1()
{
    int nbr_atoms = 256; int n_rows = nbr_atoms; int n_cols = 3; int n_unitcells = 4;
    
    // Testing lattice parametersdhwudh
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
    double lattice_param = 4.03; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);

    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling, write_not_append;
    double temp_eq, press_eq, tau_T, tau_P;

    // Production run
    end_time = 10; dt = 1e-2;
    temp_scaling = false, press_scaling = false;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    tau_T = 1 * dt; tau_P = 1*dt;
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);
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
    double lattice_param = 4.03; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling, write_not_append;
    double temp_eq, press_eq, tau_T, tau_P;

    // end_time = 25; dt = 1e-3; 
    // tau_T = 100*dt; tau_P = 50*dt;

    // Equalibration run
    end_time = 25; dt = 1e-2;
    temp_scaling = true; press_scaling = true;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    tau_T = 100*dt; tau_P = 300*dt; //*dt; 50,5
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);


    // Production run
    end_time = 5; dt = 1e-2;
    temp_scaling = false, press_scaling = false;
    temp_eq = 773.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = true;
    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);
}

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
    double lattice_param = 4.03; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling, write_not_append;
    double temp_eq, press_eq, tau_T, tau_P;


    // Melting run
    end_time = 20; dt = 1e-3;
    temp_scaling = true; press_scaling = true;
    temp_eq = 5000; press_eq = 1; //773.15 K och 1 Bar
    tau_T = 100*dt; tau_P = 300*dt; //50*dt;
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);

    // Cooling run
    end_time = 20; dt = 1e-3;
    temp_scaling = true; press_scaling = true;
    temp_eq = 973.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = false;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);


    // Production run
    end_time = 20; dt = 1e-3;
    temp_scaling = false, press_scaling = false;
    temp_eq = 973.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);

}


void
H1_task6()
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
    double lattice_param = 4.03; // True is around 4.0478 Å (Masahiko Morinaga, https://bit.ly/3ERRFt3)
    double cell_length = 4 * lattice_param;
    
    // Initialice and displace fcc
    init_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param); // 4 unit cells in each direction
    
    displace_fcc((double (*)[3]) position, (int) n_unitcells, (double) lattice_param);
    
    // Declaring parameters for velocity verlet. 
    // If temp/press_scaling = false scaling is turned off and scaling factors remains = 1
    int end_time; double dt;
    bool temp_scaling, press_scaling, write_not_append;
    double temp_eq, press_eq, tau_T, tau_P;


    // Melting run
    end_time = 20; dt = 1e-3;
    temp_scaling = true; press_scaling = true;
    temp_eq = 5000; press_eq = 1; //773.15 K och 1 Bar
    tau_T = 100*dt; tau_P = 300*dt; //*dt; 50,5
    write_not_append = true;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);

    // Cooling run
    end_time = 20; dt = 1e-3;
    temp_scaling = true; press_scaling = true;
    temp_eq = 973.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = false;

    cell_length = velocity_verlet((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T);


    // Production run
    end_time = 20; dt = 1e-3;
    temp_scaling = false, press_scaling = false;
    temp_eq = 973.15; press_eq = 1; //773.15 K och 1 Bar
    write_not_append = true;

    int number_of_bins = 150;
    double *radial_distribution_vector = calloc(sizeof(double),  number_of_bins);
    char filename_radial_dist[] = {"../csv/radial_distribution.csv"};
    double normalisation_factor_radial_dist = (double) nbr_atoms*end_time/(dt);


    cell_length = velocity_verlet_deluxe((double (*)[3]) position, (double (*)[3]) velocity, (double) lattice_param, (double) cell_length, (int) end_time, (double) dt, (int) n_cols, (int) nbr_atoms, \
                    (bool) temp_scaling, (bool) press_scaling, (double) temp_eq, (double) press_eq, (bool) write_not_append, (double) tau_P, (double) tau_T, radial_distribution_vector, number_of_bins);

    
    for(int bin = 0; bin<number_of_bins; ++bin)
    {   

        //also need to divide by number of time steps
        radial_distribution_vector[bin] /=normalisation_factor_radial_dist; 
    }

    bool is_empty = true;
    save_vector_to_csv(radial_distribution_vector, number_of_bins, filename_radial_dist, is_empty);
    free(radial_distribution_vector);
}


