// lagrange_test.c
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2017-04-15 12:42:24 (jmiller)>

// This is code for performing a simple test of the lagrange library.
// We test the 2D version of the interpolation. 

#include "lagrange.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define XMIN -1.0
#define XMAX 1.0
#define NX_COARSE 21
#define NX_FINE 101
#define DX_COARSE (XMAX-XMIN)/(NX_COARSE-1)
#define DX_FINE (XMAX-XMIN)/(NX_FINE-1)
#define FILENAME_FULL "test_full.out"
#define FILENAME_FO "test_fo.out"

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
  double z_interp_full[NX_FINE*NX_FINE];
  double z_interp_fo[NX_FINE*NX_FINE];
  double error_fo,error_full;
  double l2_error_fo=0.;
  double l2_error_full=0.;
  double l2_error_denom=0.;
  double max_error_full=0.;
  double max_error_fo=0.;

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

  printf("We compare to the analytic function\n"
	 "x^4 - y^4 + 0.5*(x*y)**4\n"
	 "on the domain [-1,1]x[-1,1].\n"
	 "For FULL Lagrange polynomials of order:\n"
	 "\t%d in x\n\t%d in y\n"
	 "and FIXED ORDER Lagrange polynomials of order:\n"
	 "\t%d in x\n\t%d in y\n",
	 NX_COARSE,NX_COARSE,order_x,order_y);

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
      z_fine[iflat] = xyn(x_fine[i],x_fine[j]);
    }
  }

  printf("Interpolating.\n");
  for (int i = 0; i < NX_FINE; i++) {
    for (int j = 0; j < NX_FINE; j++) {
      iflat = index_2D_to_1D(i,NX_FINE,j,NX_FINE);
      z_interp_fo[iflat] = lagrange_interp_2Dfo(x_fine[i],
						x_fine[j],
						order_x,order_y,
						x_coarse, NX_COARSE,
						x_coarse, NX_COARSE,
						z_coarse);
      z_interp_full[iflat] = lagrange_interp_2D(x_fine[i],x_fine[j],
						x_coarse, NX_COARSE,
						x_coarse, NX_COARSE,
						z_coarse);
      error_fo = z_interp_fo[iflat] - z_fine[iflat];
      error_full = z_interp_full[iflat] - z_fine[iflat];
      l2_error_fo += error_fo*error_fo;
      l2_error_full += error_full*error_full;
      l2_error_denom += 1.*1.;
      max_error_fo = fmax(error_fo,max_error_fo);
      max_error_full = fmax(error_full,max_error_full);
    }
  }
  l2_error_fo = sqrt(l2_error_fo/l2_error_denom);
  l2_error_full = sqrt(l2_error_full/l2_error_denom);
  printf("\tMax error = %lf\n\tL2 error = %lf\n"
	 "\tMax error fixed order %lf\n"
	 "\tL2 error fixed order %lf\n",
	 max_error_full,l2_error_full,
	 max_error_fo,l2_error_fo);

  printf("Creating output files.\n");
  FILE *fp_full;
  FILE *fp_fo;
  fp_full = fopen(FILENAME_FULL,"w");
  fp_fo = fopen(FILENAME_FO,"w");
  fprintf(fp_full,"#i\tj\tx\ty\tz\tz_interp\n");
  fprintf(fp_fo,"#i\tj\tx\ty\tz\tz_interp\n");
  fprintf(fp_full,"#NX = NY = %d\n",NX_FINE);
  fprintf(fp_fo,"#NX = NY = %d\n",NX_FINE);
  for (int i = 0; i < NX_FINE; i++) {
    for (int j = 0; j < NX_FINE; j++) {
      iflat = index_2D_to_1D(i,NX_FINE,j,NX_FINE);
      fprintf(fp_fo,"%d\t%d\t%lf\t%lf\t%lf\t%lf\n",
	      i,j,x_fine[i],x_fine[j],
	      z_fine[iflat],z_interp_fo[iflat]);
      fprintf(fp_full,"%d\t%d\t%lf\t%lf\t%lf\t%lf\n",
	      i,j,x_fine[i],x_fine[j],
	      z_fine[iflat],z_interp_full[iflat]);
    }
  }
  fclose(fp_fo);
  fclose(fp_full);
  free(x_coarse);
  free(z_coarse);
  return 0;
}
