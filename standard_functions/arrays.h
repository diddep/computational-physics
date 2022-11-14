
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
 *
 * Multiply v1 by constant a
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
void
constant_multiplication(
        double *res,
        double *v1,
        double a,
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

double
L2_norm(
        double *v1,
        unsigned int len
);