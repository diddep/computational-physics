#pragma once

/* ****************************************************************************
 * Function that calculates the forces on all atoms in units of [eV/Å]. the 
 * forces are stored in f which should be a matrix of size N x 3, where N is
 * the number of atoms and column 1,2 and 3 correspond to the x,y and z
 * component of the force respectively. pos should be a matrix containing the
 * positions of  all the atoms, L is the length of the supercell and N is the
 * number of atoms.
 *
 * Example usage:
 *
 * >>> get_forces_AL(f,pos, L, N);
 *
 * ****************************************************************************/
void get_forces_AL(double[][3] , double[][3], double, int);


/* ****************************************************************************
 * Function that calculates the potential energy in units of [eV]. pos should be
 * a matrix containing the positions of all the atoms, L is the length of the 
 * supercell and N is the number of atoms.
 *
 * Example usage:
 *
 * >>> double energy = get_energy_AL(pos, L, N);
 *
 * ****************************************************************************/
double get_energy_AL(double[][3], double, int);



/* ****************************************************************************
 * Function that calculates the virial in units of [eV]. pos should be a matrix
 * containing the positions of all the atoms, L is the length of the supercell 
 * and N is the number of atoms.
 *
 * Example usage:
 *
 * >>> double virial = get_virial_AL(pos, L, N);
 *
 * ****************************************************************************/
double get_virial_AL(double[][3], double, int);
