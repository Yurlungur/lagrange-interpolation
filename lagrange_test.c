// lagrange_test.c
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2017-03-27 18:19:31 (jmiller)>

// This is code for performing a simple test of the lagrange library.
// We test the 2D version of the interpolation. 

#include "lagrange.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define XMIN -1.0
#define XMAX 1.0
#define NX_COARSE 20
#define NX_FINE 100
#define DX_COARSE (XMAX-XMIN)/NX_COARSE
#define DX_FINE (XMAX-XMIN)/NX_FINE
#define FILENAME "test.out"

double xyn(double x, double y) {
  return pow(x,4.) - pow(y,4.) + 0.5*pow(x*y,3.);
}

int main(int argc, char* argv[]) {
  int iflat;
  int nz_coarse = NX_COARSE*NX_COARSE;
  double* x_coarse;
  double* z_coarse;
  double x_fine[NX_FINE];
  double z_fine[NX_FINE*NX_FINE];
  double z_interp[NX_FINE*NX_FINE];
  double error;
  double l2_error=0.;
  double l2_error_denom=0.;
  double max_error=0.;

  if (argc < 3) {
    printf("usage %s order_x order_y\n", argv[0]);
    return 1;
  }
  int order_x = atoi(argv[1]);
  int order_y = atoi(argv[2]);
  if (order_x < 0 || order_y < 0) {
    printf("Interpolation order must be at least 0\n");
    return 1;
  }
  if (order_x >= NX_COARSE || order_y >= NX_COARSE) {
    printf("Interpolation order can be no greater than %d\n",
	   NX_COARSE-1);
  }

  printf("Beginning test with:\n"
	 "\torder in x = %d\n"
	 "\torder in y = %d.\n",
	 order_x,order_y);
  printf("We compare to the analytic function\n"
	 "x^4 - y^4 + 0.5*(x*y)**4\n"
	 "on the domain [-1,1]x[-1,1].\n");

  printf("Preparing grids.\n");
  x_coarse = (double *)malloc(sizeof(double)*NX_COARSE);
  z_coarse = (double *)malloc(sizeof(double)*nz_coarse);
  printf("Writing to grids.\n");
  for (int i = 0; i < NX_COARSE; i++) {
    x_coarse[i] = XMIN + DX_COARSE*i;
  }
  for (int i = 0; i < NX_FINE; i++) {
    x_fine[i] = XMIN + DX_FINE*i;
  }
  for (int i = 0; i < NX_COARSE; i++) {
    for (int j = 0; j < NX_COARSE; j++) {
      iflat = index_2D_to_1D(i,NX_COARSE,j,NX_COARSE);
      z_coarse[iflat] = xyn(x_coarse[i],x_coarse[j]);
    }
  }
  for (int i = 0; i < NX_FINE; i++) {
    for (int j = 0; j < NX_FINE; j++) {
      iflat = index_2D_to_1D(i,NX_FINE,j,NX_FINE);
      z_fine[iflat] = xyn(x_coarse[i],x_coarse[j]);
    }
  }

  printf("Interpolating.\n");
  for (int i = 0; i < NX_FINE; i++) {
    for (int j = 0; j < NX_FINE; j++) {
      iflat = index_2D_to_1D(i,NX_FINE,j,NX_FINE);
      z_interp[iflat] = lagrange_interp_2Dfo(x_fine[i],
					     x_fine[i],
					     order_x,order_y,
					     x_coarse, NX_COARSE,
					     x_coarse, NX_COARSE,
					     z_coarse);
      error = z_interp[iflat] - z_fine[iflat];
      l2_error += error*error;
      l2_error_denom += 1.*1.;
      max_error = fmax(error,max_error);
    }
  }
  l2_error /= l2_error_denom;
  printf("\tMax error = %lf\n\tL2 error = %lf\n",
	 max_error,l2_error);

  printf("Creating output file.\n");
  FILE *fp;
  fp = fopen(FILENAME,"w");
  fprintf(fp,"#i\tj\tx\ty\tz\tz_interp\n");
  fprintf(fp,"#NX = NY = %d\n",NX_FINE);
  for (int i = 0; i < NX_FINE; i++) {
    for (int j = 0; j < NX_FINE; j++) {
      iflat = index_2D_to_1D(i,NX_FINE,j,NX_FINE);
      fprintf(fp,"%d\t%d\t%lf\t%lf\t%lf\t%lf\n",
	      i,j,x_fine[i],x_fine[j],
	      z_fine[iflat],z_interp[iflat]);
    }
  }
  fclose(fp);
  return 0;
}
