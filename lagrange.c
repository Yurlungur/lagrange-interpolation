// lagrange.c
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2017-03-30 21:07:11 (jmiller)>

#include "lagrange.h"
#include <assert.h>
#include <stdlib.h>

extern inline int min(int a, int b);
extern inline int max(int a, int b);
extern inline double get_dx(const double x[]);

int find_xb(double x, const double a[], int n, int b) {
  assert(b<n);
  double dx = get_dx(a);
  int ix = (x - a[0])/dx;
  ix = max(ix,b);
  ix = min(ix,n-b-1);
  return ix;
}

int find_x(double x, const double a[], int n) {
  return find_xb(x,a,n,0);
}

int flatten_index(const int index[], const int array_sizes[],
		  int d) {
  if (d == 0) return 0;
  if (d == 1) return index[0];
  if (d == 2) {
    return index[1] + array_sizes[1]*index[0];
  }
  int rec = flatten_index(index,array_sizes,d-1);
  return index[d-1] + array_sizes[d-1]*rec;
}

int index_2D_to_1D(int i, int ni,
		   int j, int nj) {
  return j + nj*i;
}

int index_3D_to_1D(int i, int ni,
		   int j, int nj,
		   int k, int nk) {
  return k + nk*(j + nj*i);
}

double lagrange_basis(int j, double x,
		      const double x_values[],
		      int num_x) {
  double l_j = 1;
  double r;
  assert( 0 <= j && j < num_x && "j is in bounds" );
  if (num_x  == 0 || num_x == 1) return 1.;
  for (int m = 0; m < num_x; m++) {
    if (m != j) {
      r = (x - x_values[m])/(x_values[j] - x_values[m]);
      l_j *= r;
    }
  }
  return l_j;
}

double lagrange_interp_1D(double x,
			  const double x_values[],
			  int num_x,
			  const double y_values[]) {
  double out = 0;
  for (int j = 0; j < num_x; j++) {
    out += y_values[j]*lagrange_basis(j,x,x_values,num_x);
  }
  return out;
}

double lagrange_interp_1Dfo(double x, int order,
			    const double x_values[], int num_x, 
			    const double y_values[]) {
  assert (order % 2 == 0);
  assert (num_x > order);
  double mx[order+1];
  double my[order+1];
  int i;
  int ix = find_xb(x,x_values,num_x,order/2);
  for (i = 0; i < order+1; i++) {
    mx[i] = x_values[ix + i - order/2];
    my[i] = y_values[ix + i - order/2];
  }
  return lagrange_interp_1D(x, mx, order+1, my);
}

double lagrange_interp_2D(double x, double y,
			  const double x_values[], int num_x,
			  const double y_values[], int num_y,
			  const double z_values[]) {
  double out = 0;
  double lx,ly,z;
  for (int i = 0; i < num_x; i++) {
    lx = lagrange_basis(i,x,x_values,num_x);
    for (int j = 0; j < num_y; j++) {
      ly = lagrange_basis(j,y,y_values,num_y);
      z = z_values[index_2D_to_1D(i,num_x,j,num_y)];
      out += z*lx*ly;
    }
  }
  return out;
}

double lagrange_interp_2Dfo(double x, double y,
			    int order_x, int order_y,
			    const double x_values[], int num_x,
			    const double y_values[], int num_y,
			    const double z_values[]) {
  assert(order_x % 2 == 0);
  assert(order_y % 2 == 0);
  assert(num_x > order_x);
  assert(num_y > order_y);
  double mx[order_x+1];
  double my[order_y+1];
  double mz[(order_x+1)*(order_y+1)];
  int ix = find_xb(x,x_values,num_x,order_x/2);
  int iy = find_xb(y,y_values,num_y,order_y/2);
  int iz_small,iz_big;
  for (int i = 0; i < order_x + 1; i++) {
    mx[i] = x_values[ix + i - order_x/2];
    for (int j = 0; j < order_y + 1; j++) {
      my[j] = y_values[iy + j - order_y/2];
      iz_small = index_2D_to_1D(i,order_x+1,
				j,order_y+1);
      iz_big = index_2D_to_1D(ix + i - order_x/2, num_x,
			      iy + j - order_y/2, num_y);
      mz[iz_small] = z_values[iz_big];
    }
  }
  return lagrange_interp_2D(x,y,
			    mx,order_x+1,
			    my,order_y+1,
			    mz);
}

double lagrange_interp_3D(double x, double y, double z,
			  const double x_values[], int num_x,
			  const double y_values[], int num_y,
			  const double z_values[], int num_z,
			  const double f_values[]) {
  double out = 0;
  double lx,ly,lz,f;
  for (int i = 0; i < num_x; i++) {
    lx = lagrange_basis(i,x,x_values,num_x);
    for (int j = 0; j < num_y; j++) {
      ly = lagrange_basis(i,y,y_values,num_y);
      for (int k = 0; k < num_z; k++) {
	f = f_values[index_3D_to_1D(i,num_x,j,num_y,k,num_z)];
	lz = lagrange_basis(i,z,z_values,num_z);
	out += f*lx*ly*lz;
      }
    }
  }
  return out;
}

double lagrange_interp_3Dfo(double x, double y, double z,
			    int order_x, int order_y, int order_z,
			    const double x_values[], int num_x,
			    const double y_values[], int num_y,
			    const double z_values[], int num_z,
			    const double f_values[]) {
  assert(order_x % 2 == 0);
  assert(order_y % 2 == 0);
  assert(order_z % 2 == 0);
  assert(num_x > order_x);
  assert(num_y > order_y);
  assert(num_z > order_z);
  double mx[order_x+1];
  double my[order_y+1];
  double mz[order_z+1];
  double mf[(order_x+1)*(order_y+1)*(order_z+1)];
  int ix = find_xb(x,x_values,num_x,order_x/2);
  int iy = find_xb(y,y_values,num_y,order_y/2);
  int iz = find_xb(z,z_values,num_z,order_z/2);
  int if_small,if_big;
  for (int i = 0; i < order_x + 1; i++) {
    mx[i] = x_values[ix + i - order_x/2];
    for (int j = 0; j < order_y + 1; j++) {
      my[j] = y_values[iy + j - order_y/2];
      for (int k = 0; k < order_z + 1; k++) {
	mz[j] = z_values[iz + k - order_z/2];
	if_small = index_3D_to_1D(i,order_x+1,
				  j,order_y+1,
				  k,order_z+1);
	if_big = index_3D_to_1D(ix + i - order_x/2, num_x,
				iy + j - order_y/2, num_y,
				iz + k - order_z/2, num_z);
	mf[if_small] = f_values[if_big];
      }
    }
  }
  return lagrange_interp_3D(x,y,z,
			    mx,order_x+1,
			    my,order_y+1,
			    mz,order_z+1,
			    mf);
}
