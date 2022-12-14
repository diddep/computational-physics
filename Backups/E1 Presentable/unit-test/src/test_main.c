#include <stdlib.h>
#include <string.h>

#include "test_main.h"
#include "tools.h"

/* *************************************
 * The test framework forces you
 * to work with global variables which
 * is a bad pattern.
 *
 * If you are interested in digging
 * deeper into testing, we recommend
 * cmocka whcih is a better designed
 * testing framework. However, it requiers
 * a somewhat more advanced understanding
 * of C.
 * ************************************/


/* ************************************
 * Gloabl variables makes the accesible
 * int the testing scope
 * ************************************/
#define SIZE 4
#define ARRAY_SIZE_M 4 //Carl: Switched from 3 to 4
#define ARRAY_SIZE_N 3 //Carl: Switched from 4 to 3

/* ************************************
 * Here follows the test we want
 * to run
 * ***********************************/
START_TEST(test_elementwise_addition)
{
    // assert re and correct answer is within 0.001f //
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    double *result_add = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v1[i] = 0.1*i;
	v2[i] = 0.01*i;
    }
    elementwise_addition(result_add, v1, v2, SIZE);
    ck_assert_double_eq_tol(result_add[0], 0, 1e-6);
    ck_assert_double_eq_tol(result_add[1], 0.11, 1e-6);
    ck_assert_double_eq_tol(result_add[2], 0.22, 1e-6);
    ck_assert_double_eq_tol(result_add[3], 0.33, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
    free(result_add); result_add = NULL;
}

START_TEST(test_elementwise_multiplication)
{
    // assert re and correct answer is within 0.001f //
    double * v1 = malloc(sizeof(double) * SIZE);
    double * v2 = malloc(sizeof(double) * SIZE);
    double *result_mul = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v1[i] = 0.1*i;
	v2[i] = 0.01*i;
    }
    elementwise_multiplication(result_mul, v1, v2, SIZE);
    ck_assert_double_eq_tol(result_mul[0], 0., 1e-6);
    ck_assert_double_eq_tol(result_mul[1], 0.001, 1e-6);
    ck_assert_double_eq_tol(result_mul[2], 0.004, 1e-6);
    ck_assert_double_eq_tol(result_mul[3], 0.009, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
    free(result_mul); result_mul = NULL;
}

#include <stdio.h>
START_TEST(test_dot_product)
{
    // assert re and correct answer is within 0.001f //
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v1[i] = 0.1*i;
	v2[i] = 0.01*i;
    }
    double result_dot_prd = dot_product(v1, v2, SIZE);
    ck_assert_double_eq_tol(result_dot_prd, 0.014, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
}



START_TEST(test_matrix_multiplication)
{
    double **result_array = NULL;
    double **matrix_A = NULL;
    double **matrix_B = NULL;
    create_2D_array(&result_array,
		    ARRAY_SIZE_M, ARRAY_SIZE_M); // Carl: Switched from N to M to create 3x3 result matrix (rel. to row-major)
    create_2D_array(&matrix_A,
		    ARRAY_SIZE_M, ARRAY_SIZE_N);
    create_2D_array(&matrix_B,
		    ARRAY_SIZE_N, ARRAY_SIZE_M);
    for(int i = 0; i < ARRAY_SIZE_M; ++i)
    { // Carl: Switced place between N and M so test is row major
	    for(int j = 0; j < ARRAY_SIZE_N; ++j)
        {
	        matrix_A[i][j] = (double) i;
	    }
    }

    for(int i = 0; i < ARRAY_SIZE_N; ++i)
    { // Carl: Switced place between N and M so test is row major
	    for(int j = 0; j < ARRAY_SIZE_M; ++j)
        {
	        matrix_B[i][j] = (double) j;
	    }
    }
   
    // Carl: Changed matrix_multiplication to accept two matrices of dim(m x n) and dim(n x p), and added input for p. Result is a dim(m x p) matrix.
    // now m,n,p where first vector is 
    matrix_multiplication(result_array, matrix_A, matrix_B, ARRAY_SIZE_M, ARRAY_SIZE_N, ARRAY_SIZE_M);  

    /* // LEAVING HERE, REMOVING THIS COMMENT BLOCK PRINTS RESULT ARRAY
    for(int i = 0; i < ARRAY_SIZE_M; ++i){}
	for(int j = 0; j < ARRAY_SIZE_M; ++j){
        printf("%f ",result_array[i][j]);
	}
    printf("\n");
    }
    printf("\n");
    */
    
    // assert re and correct answer is within 0.001f //
    ck_assert_double_eq_tol(result_array[0][0], 0, 1e-6);
    ck_assert_double_eq_tol(result_array[0][1], 0, 1e-6);
    ck_assert_double_eq_tol(result_array[0][2], 0, 1e-6);
    ck_assert_double_eq_tol(result_array[0][3], 0, 1e-6);
    
    ck_assert_double_eq_tol(result_array[1][0], 0.f, 1e-6);
    ck_assert_double_eq_tol(result_array[1][1], 3.f, 1e-6);
    ck_assert_double_eq_tol(result_array[1][2], 6.f, 1e-6);
    ck_assert_double_eq_tol(result_array[1][3], 9.f, 1e-6);
    
    ck_assert_double_eq_tol(result_array[2][0], 0.f, 1e-6);
    ck_assert_double_eq_tol(result_array[2][1], 6.f, 1e-6);
    ck_assert_double_eq_tol(result_array[2][2], 12.f, 1e-6);
    ck_assert_double_eq_tol(result_array[2][3], 18.f, 1e-6);
    
    ck_assert_double_eq_tol(result_array[3][0], 0.f, 1e-6);
    ck_assert_double_eq_tol(result_array[3][1], 9.f, 1e-6);
    ck_assert_double_eq_tol(result_array[3][2], 18.f, 1e-6);
    ck_assert_double_eq_tol(result_array[3][3], 27.f, 1e-6);
    
    destroy_2D_array(result_array); result_array = NULL;
    destroy_2D_array(matrix_A); matrix_A = NULL;
    destroy_2D_array(matrix_B); matrix_B = NULL;
}

START_TEST(test_vector_length)
{
    // assert re and correct answer is within 0.001f //
    double *v = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v[i] = 0.1*i;
    }
    double v_norm = vector_norm(v, SIZE);
    ck_assert_double_eq_tol(v_norm, 0.3741657, 1e-6);
    free(v); v = NULL;
}

START_TEST(test_vector_normed)
{
    // assert re and correct answer is within 0.001f //
    double *v_normed = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v_normed[i] = 0.1*i;
    }
    normalize_vector(v_normed, SIZE);
    ck_assert_double_eq_tol(v_normed[0], 0.f, 1e-6);
    ck_assert_double_eq_tol(v_normed[1], 0.26726124f, 1e-6);
    ck_assert_double_eq_tol(v_normed[2], 0.53452248f, 1e-6);
    ck_assert_double_eq_tol(v_normed[3], 0.80178373f, 1e-6);
    free(v_normed); v_normed = NULL;
}

START_TEST(test_average_and_std)
{
    // assert re and correct answer is within 0.001f //
    double *v = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v[i] = 0.01*i;
    }
    double v_average = average(v, SIZE);
    double v_std = standard_deviation(v, SIZE);
    ck_assert_double_eq_tol(v_std, 0.011180339f, 1e-6);
    ck_assert_double_eq_tol(v_average, 0.015f, 1e-6);
    free(v); v = NULL;
}

START_TEST(test_distance_between_vectors)
{
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
	v1[i] = 0.1*i;
	v2[i] = 0.01*i;
    }
    double v1_v2_distance = distance_between_vectors(v1, v2, SIZE);
    ck_assert_double_eq_tol(v1_v2_distance, 0.33674916, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
}


int
main()
{
    test_setup("testing", "core");
    
    // Tests
    add_test(test_elementwise_addition);
    add_test(test_elementwise_multiplication);
    add_test(test_dot_product);
    add_test(test_matrix_multiplication);
    add_test(test_vector_length);
    add_test(test_vector_normed);
    add_test(test_average_and_std);
    add_test(test_distance_between_vectors);
    
    test_teardown();
    return 0;
}
