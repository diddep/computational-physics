#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"
    
void
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
	for(int ix=0; ix < len; ix++)
	{
		res[ix] = v1[ix] + v2[ix];
	}
}

void
elementwise_subtraction(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
	for(int ix=0; ix < len; ix++)
	{
		res[ix] = v1[ix] - v2[ix];
	}
}

void
elementwise_multiplication(
			   double *res,
			   double *v1,
			   double *v2,
			   unsigned int len
	                  )
{
	for(int ix = 0; ix < len; ix++)
	{
		res[ix] = v1[ix] * v2[ix];
	}
}

double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   )
{
    double result = 0;
	for(int ix = 0; ix < len; ix++)
	{
		result += v1[ix] * v2[ix];
	}
	return result;
}

void
create_2D_array(
		double ***array,
		unsigned int row_size,
		unsigned int column_size // Carl: Changed place between row_size and column_size
	       )
{
	double *asentries = (double*) malloc(row_size * column_size * sizeof(double)); // Use calloc ,sizeof(double)
	*array = (double **) malloc(row_size * sizeof(double *));
	
	for(int row = 0; row < row_size; row++)
	{
		(*array)[row] = asentries + column_size * row;
	}
	
	for(int row = 0; row < row_size; row++)
	{
		for(int col = 0; col < column_size; col++)
			{
				(*array)[row][col] = 0;
			}
	}
}

void
destroy_2D_array(
		 double **array
		){
	free(*array);
	free(array);
}

/*
void
matrix_multiplication(
		      double **result,
		      double **v1,
		      double **v2,
		      unsigned int m,
		      unsigned int n
		     )
{
	printf("V1: \n");
	for(int ix = 0; ix < m; ++ix){
	for(int jx = 0; jx < n; ++jx){
	    v1[ix][jx] = ix;
        printf("%f ",v1[ix][jx]);
	}
    printf("\n");
    }
    printf("\n");
	printf("V2: \n");
    for(int ix = 0; ix < n; ++ix){
	for(int jx = 0; jx < m; ++jx){
	    v2[ix][jx] = jx;
        printf("%f ",v2[ix][jx]);
	}
    printf("\n");
    }
	printf("result: \n");
	for(int ix = 0; ix < m; ++ix){
	for(int jx = 0; jx < m; ++jx){ //Change later to allow general dimensions mxn X nxp
	for(int kx = 0; kx < m-1; ++kx){
	    result[ix][jx] += v1[ix][kx] * v2[kx][jx];
	}
	printf("%f ",result[ix][jx]);
    }
	printf("\n");
	}
}
*/
void
matrix_multiplication(
		      double **result,
		      double **v1,
		      double **v2,
		      unsigned int m,
		      unsigned int n,
			  unsigned int p
		     )
{
	for(int ix = 0; ix < m; ix++)
	{
		for(int jx = 0; jx < p; jx++)
		{
			for(int kx = 0; kx < n; kx++)
			{
	    		result[ix][jx] += v1[ix][kx] * v2[kx][jx];
			}
    	}
	}
}

double
vector_norm(
	    double *v1,
	    unsigned int len
	   )
{
	double result = 0;
	// Euclidian L2 norm
	for(int ix =0; ix < len; ix++)
	{
		result += v1[ix]*v1[ix];
	}
	result = sqrt(result);

    return result;
}


void
normalize_vector(
		 double *v1,
		 unsigned int len
		)
{
	double norm = vector_norm(v1, len);
	for(int ix =0; ix < len; ix++)
	{
		v1[ix] = v1[ix]/norm;
	}
}

double
average(
	double *v1,
	unsigned int len
       )
{
	double result = 0;
	for(int i =0; i < len; i++)
	{
		result += v1[i];
	}
	result /= len;

    return result;
}


double
standard_deviation(
	    	   double *v1,
	    	   unsigned int len
	          )
{
	double ave = average(v1, len);
	double diff[len];
	
	for(int ix = 0; ix < len; ix++)
	{
		diff[ix] = v1[ix] - ave;
	}

	double std = 0;
	for(int ixx = 0; ixx < len; ixx++)
	{
		std += diff[ixx]*diff[ixx];
	}
	std /= len;
	std = sqrt(std);
    return std;
}

double
distance_between_vectors(
			 double *v1,
			 double *v2,
			 unsigned int len
	                )
{
	double *distance = calloc(sizeof(double), len);
	elementwise_subtraction(distance, v1, v2, len);
	double result = vector_norm(distance, len);
	free(distance);

    return result;
}
