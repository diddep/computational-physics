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
	double *asentries = (double*) calloc(row_size * column_size, sizeof(double)); // Use calloc ,sizeof(double)
	*array = (double **) calloc(row_size, sizeof(double *));
	
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

void
print_2D_array(
		double **array,
		unsigned int row_size,
		unsigned int column_size // Carl: Changed place between row_size and column_size
	       )
{
	for(int ix = 0; ix < row_size; ix++){
			for(int jx = 0; jx < column_size; jx++){
				if( ix % 1000 == 0){
					printf("%f ", array[ix][jx]);
				}
			}
			printf("\n");
		}
}

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

void
print_vector(
        double *vec,  // Vector to print
        unsigned int ndims  // Number of dimensions
        )
{
    for (int i = 0; i < ndims; i++) {
        printf("%10.5f ", vec[i]);
    }
}

int save_vector_to_csv(
        double *vec,  // Vector to save
        unsigned int ndims,  // Number of dimensions
        char *filename // filename
        )
{
    FILE *fp1;
    fp1 = fopen(filename, "w"); // Create a file
    if (fp1 == NULL)
    {
        printf("Error while opening the file.\n");
        return 1;
    }
    //fprintf(fp1, "%10.5f,%10.5f,%10.5f", vec[0], vec[1], vec[2]);
    for(int i =0; i<ndims; ++i){
        fprintf(fp1, "%10.5f", vec[i]);
        fprintf(fp1, "\n");
    }

    fprintf(fp1,"\n");

    fclose(fp1);
    return 0;
}

int save_matrix_to_csv(
        double **matrix,  // Matrix to save
        unsigned int nrows,  // Number of dimensions
		unsigned int ncols,  // Number of dimensions
        char *filename // filename
        )
{
    FILE *fp1;
    fp1 = fopen(filename, "w"); // Create a file
    if (fp1 == NULL)
    {
        printf("Error while opening the file.\n");
        return 1;
    }

	for(int ix = 0; ix<nrows; ix++)
	{
		for(int jx = 0; jx<ncols; jx++)
		{
			if(jx == ncols - 1)
			{
				fprintf(fp1, "%10.5f", matrix[ix][jx]);
			} else {
				fprintf(fp1, "%10.5f, ", matrix[ix][jx]);
			}
			
		}
	fprintf(fp1,"\n");
	}

    fclose(fp1);
    return 0;
}
