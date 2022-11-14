/* **********************************************
 *
 * Add v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
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
void
create_2D_array(
		double ***array,
		unsigned int column_size,
		unsigned int row_size
	       );

void
destroy_2D_array(
		 double **array
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
