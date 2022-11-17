/******************************************************************************
 * E1code3
 ******************************************************************************
 * reads h(t) = a*cos(2*pi*f*t + phi)
 * runs the fft of h(t)
 *
 *
 * Compile me as:
 * clang -c fft.c -o fft.o -lgsl -lgslcblas
 * clang E1code3 fft.o -o <execuutable name> -lgsl -lgslcblas
 */

/************************************************************
 * Includes
 ************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../include/fft.h"  // interface to fft routine
#include "../src/fft.c"

/*************************************************************
 * Macro defines
 *************************************************************/
#define N_POINTS 250 //original 250
#define dt 0.1

/******************************************************************************
 * Helper functions
 *****************************************************************************/
/*
 * reads time_array and signal data from file
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
*/
void read_data(char *fname, double *time_array, double *signal)
{
    FILE *fp = fopen(fname, "r");

    /* if file no found
     * error out and exit code 1
     */
    if(fp == NULL){
	perror("error:");
	exit(1);
    }

    /* skip header */
    fseek(fp, strlen("time, signal\n"), SEEK_SET);
    char line[128] = {0};
    char *token;
    int i = 0;
    while(fgets(line, sizeof(line), fp) != NULL){
	token = strtok(line, ",");
	time_array[i] = strtod(token, NULL);
	token = strtok(NULL, ",");
	signal[i] = strtod(token, NULL);
	i++;
	memset(line, 0, sizeof(line));
	token = NULL;
    }
    fclose(fp);
}

/*
 * constructs time array
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
 * @n_points - number of points
*/
void write_to_file(char *fname, double *frequencies,
		   double *spectrum, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%e,%e\n", frequencies[i], spectrum[i]);
    }
    fclose(fp);
}



/**************************************************************
 * Main routine
 **************************************************************/
int main(int argc, char **argv)
{
    double time_array[N_POINTS];
    double signal[N_POINTS];
    read_data("signal.csv", time_array, signal);
    
    /*
     * Construct array with frequencies
     */ 
    double frequencies[N_POINTS];
    for(int i = 0; i < N_POINTS; i++){
	    frequencies[i] = i / (dt * N_POINTS);
    }

    /*
     * Do the fft
     */
    double fftd_data[N_POINTS];
    powerspectrum(signal, fftd_data, N_POINTS);
    powerspectrum_shift(fftd_data, N_POINTS);
    fft_freq_shift(frequencies, dt, N_POINTS);

    /*
     * Dump fft and frequencies to file
     */
    write_to_file("powerspectrum_shift.csv", fftd_data, frequencies, N_POINTS);
    return 0;
}
