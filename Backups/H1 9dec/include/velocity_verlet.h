

double
velocity_verlet(double positions[][3], double v[][3], double lattice_param, double cell_length, int end_time, double dt, int n_cols, int nbr_atoms, bool temp_scaling, bool press_scaling, double temp_eq, double press_eq, bool write_not_append, double tau_P, double tau_T);

double
velocity_verlet_deluxe(double positions[][3], double v[][3], double lattice_param, double cell_length, int end_time, \
double dt, int n_cols, int nbr_atoms, bool temp_scaling, bool press_scaling, double temp_eq, double press_eq, bool write_not_append, \
double tau_P, double tau_T, double *radial_histogram_vector, int number_of_bins);
