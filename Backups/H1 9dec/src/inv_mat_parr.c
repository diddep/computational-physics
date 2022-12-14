#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> // include the OpenMP library

#define N 3 // size of the matrix

// Function to calculate the inverse of a matrix
void inverse(double matrix[N][N])
{
    // Find the determinant of the matrix
    double det = 0;
    for (int i = 0; i < N; i++)
        det += matrix[0][i] * (matrix[1][(i+1)%N] * matrix[2][(i+2)%N] - matrix[1][(i+2)%N] * matrix[2][(i+1)%N]);

    // Create a new matrix to hold the inverse
    double inverse[N][N];

    // Find the inverse of the matrix using parallelism with OpenMP
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inverse[i][j] = ((matrix[(j+1)%N][(i+1)%N] * matrix[(j+2)%N][(i+2)%N]) - (matrix[(j+1)%N][(i+2)%N] * matrix[(j+2)%N][(i+1)%N])) / det;

    // Print the inverse of the matrix
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%lf ", inverse[i][j]);
        printf("\n");
    }
}

int main()
{
    // Example matrix
    double matrix[N][N] = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};

    inverse(matrix);
    return 0;
}