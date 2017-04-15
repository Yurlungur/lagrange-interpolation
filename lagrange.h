// lagrange.h
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2017-04-15 12:46:27 (jmiller)>

// This is a library for performing Lagrange interpolation of
// polynomials on (not necessarily uniform) Cartesian product grids.
// Note we do not perform a bounds check. If you decide to use
// these polynomials to extrapolate, that's on you.
// ----------------------------------------------------------------------

#ifndef LAGRANGE_H
#define LAGRANGE_H

// max and min because the ones in cstdlib can often be borked
inline int min(int a, int b) {return a<b ? a : b;}
inline int max(int a, int b) {return a<b ? b : a;}

// Gets the spacing of an ordered evenly spaced array of length at
// least 2
inline double get_dx(const double x[]) { return x[1] - x[0]; }

// get index of x in a sorted, evenly spaced array a of length n.
// If x is not in the array, the returned index will be either 0
// or n-1, depending.
// Optionally, adds a buffer b so that x lands within the indexes
// [b, n-b-1].
int find_xb(double x, const double a[], int n, int b);
int find_x(double x, const double a[], int n);

// Maps an d-dimensional index to a 1D index
// array sizes gives the size of the tensor in each direction
int flatten_index(const int index[], const int array_sizes[], int d);

// Maps a 2D index (i,j) to a 1D index. ni,nj are size of array in i,j
// directions
int index_2D_to_1D(int i, int ni,
		   int j, int nj);
// Maps a 3D index (i,j,k) to a 1D index. ni,nj,nk are size of array
// in i,j,k directions.
int index_3D_to_1D(int i, int ni,
		   int j, int nj,
		   int k, int nk);

// Evaluates the j'th Lagrange basis function l_j at point x. The
// basis function is constructed from the num_x interpolation points
// and is thus of order k = num_x-1.
double lagrange_basis(int j, double x, const double x_values[],
		      int num_x);

// Evaluates the Lagrange interpolant of y at point x for the data
// (x_values[i],y_values[i]), where 0 <= i < num_x
double lagrange_interp_1D(double x,
			  const double x_values[], int num_x,
			  const double y_values[]);

// Evaluates the Lagrange interpolant of y at point x for the data
// (x_values[i],y_values[i])
// Now assumes x_values and y_values are ordered and large
// and that the user wants an interpolant of a given order.
// requires order is even
double lagrange_interp_1Dfo(double x, int order,
			    const double x_values[], int num_x,
			    const double y_values[]);

// Evaluates the Lagrange interpolant of z at point (x,y) for the data
// (x_values[i],y_values[j],z_values[i][j]), where 0 <= i < num_x
// and 0 <= j < num_y.
// Assumes z_values is a 1D array with column-major ordering. i.e.,
// z_values[i][j] is shorthand for z_values[index_2D_to_1D(i,j)]
// If you are using C and row-major ordering, this should be automatic
// for statically declared arrays.
double lagrange_interp_2D(double x, double y,
			  const double x_values[], int num_x,
			  const double y_values[], int num_y,
			  const double z_values[]);

// Evaluates the Lagrange interpolant of z at point (x,y) for the data
// (x_values[i],y_values[j],z_values[i][j]), where 0 <= i < num_x
// and 0 <= j < num_y.
// Assumes z_values is a 1D array with column-major ordering. i.e.,
// z_values[i][j] is shorthand for z_values[index_2D_to_1D(i,j)]
// Now assumes that x_values[] and y_values[] are very large
// and interpolates with order_x in x and order_y in y.
// Works only for even order.
double lagrange_interp_2Dfo(double x, double y,
			    int order_x, int order_y,
			    const double x_values[], int num_x,
			    const double y_values[], int num_y,
			    const double z_values[]);

// Evaluates the Lagrange interpolant of f at point (x,y,z) for the data
// (x_values[i],y_values[j],z_values[k],f_values[i,j,k]),
// where 0 <= i,j,k < num_x,num_y,num_z
// Assumes f_values is a 1D array with column-major ordering. i.e.,
// z_values[i][j][k] is shorthand for z_values[index_2D_to_1D(i,j,k)]
double lagrange_interp_3D(double x, double y, double z,
			  const double x_values[], int num_x,
			  const double y_values[], int num_y,
			  const double z_values[], int num_z,
			  const double f_values[]);

// Evaluates the Lagrange interpolant of f at point (x,y,z) for the
// data
// (x_values[i],y_values[j],z_values[k],f_values[i,j,k]),
// where 0 <= i,j,k < num_x,num_y,num_z
// Assumes f_values is a 1D array with column-major ordering. i.e.,
// z_values[i][j][k] is shorthand for z_values[index_2D_to_1D(i,j,k)]
// Now assumes that x_values,y_values, and z_values are very large
// and interpolates with order_x in x, order_y in y, and order_z in z.
// Works only for even order.
double lagrange_interp_3Dfo(double x, double y, double z,
			    int order_x, int order_y, int order_z,
			    const double x_values[], int num_x,
			    const double y_values[], int num_y,
			    const double z_values[], int num_z,
			    const double f_values[]);

#endif // LAGRANGE_H
