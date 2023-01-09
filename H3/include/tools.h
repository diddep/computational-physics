/* **********************************************
 *
 * Add v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
#include <stdbool.h>
void
elementwise_addition(
        double *res,
        double *v1,
        double *v2,
        unsigned int len
);

/* **********************************************
 * Carl: Own addition
 * Subtract v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/

void
elementwise_subtraction(
        double *res,
        double *v1,
        double *v2,
        unsigned int len
);

/* **********************************************
 *
 * Multiply v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
void
elementwise_multiplication(
        double *res,
        double *v1,
        double *v2,
        unsigned int len
);

/* **********************************************
 *
 * Calculate the dot product between
 * v1 and v2
 *
 * the result is returned as a double
 *
 * **********************************************/
double
dot_product(
        double *v1,
        double *v2,
        unsigned int len
);

/* **********************************************
 *
 * Allocate the memory to a 2D array
 *
 * **********************************************/

double**
create_2D_array(
        unsigned int row_size,
        unsigned int column_size
);

void
destroy_2D_array_pointers(
        double **array
);

void
destroy_2D_array(
        double **array,
        int n_rows
);

void
print_2D_array(
        double **array,
        unsigned int row_size,
        unsigned int column_size
);


/* **********************************************
 *
 * Calculate the matrix product between
 * v1 of dim(m x n) and v2 of dim(n x p)
 *
 * the result is dim(m x p)-matrix with elements of type double
 *
 * **********************************************/
void
matrix_multiplication(
        double **result,
        double **v1,
        double **v2,
        unsigned int m,
        unsigned int n,
        unsigned int p
);

double
vector_norm(
        double *v1,
        unsigned int len
);

void
normalize_vector(
        double *v1,
        unsigned int len
);

double
average(
        double *v1,
        unsigned int len
);

double
standard_deviation(
        double *v1,
        unsigned int len
);

double
distance_between_vectors(
        double *v1,
        double *v2,
        unsigned int len
);

void
print_vector(
        double *vec,
        unsigned int ndims
);

int
save_vector_to_csv(
        double *vec,
        unsigned int ndims,
        char *filename,
        bool is_empty
);

int save_transposedvector_to_csv(
        double *vec,  // Vector to save
        unsigned int ndims,  // Number of dimensions
        char *filename, // filename
        bool is_empty
);

int
save_matrix_to_csv(
        double **matrix,
        unsigned int nrows,
        unsigned int ncols,
        char *filename
);

gsl_rng *
init_random_num_generator();

double variance(
        double *quantity_vec,
        int number_of_elements
        );
