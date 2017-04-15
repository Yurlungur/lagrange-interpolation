// lagrance_convergence.c
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2017-04-15 14:54:36 (jmiller)>

// This performs a convergence test for the lagrange library.

#include "lagrange.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// in case the implementation lacks this
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define XMIN -1.0*M_PI
#define XMAX M_PI
#define NUM_GRIDS 4
#define NUM_RESOLUTIONS 3
#define NUM_ORDERS 3
#define FILENAME_SLICE "convergence_slice.out"
#define FILENAME_NORM "convergence_norm.out"
#define INFINITY_NORM 1
#define L2_NORM 2

double f_ana(double x, double y, int kx, int ky) {
  assert (kx >= 0 && ky >= 0);
  double out = sin(kx*x)*cos(ky*y);
  return out;
}

double get_dx_from_bounds(double xmin, double xmax, int nx) {
  return (xmax - xmin) / (nx - 1);
}

double infinity_norm(double* error, int size) {
  double out = 0;
  int i;
  for (i = 0; i < size; ++i) {
    out = fmax(fabs(error[i]),out);
  }
  return out;
}

double l2_norm(double* error,
	       double dx, double dy,
	       int size) {
  int i;
  double dv = dx*dy;
  double v = (XMAX - XMIN)*(XMAX-XMIN);
  double out;
  double integral = 0;
  for (i = 0; i < size; ++i) {
    integral += error[i]*error[i]*dv;
  }
  out = sqrt(integral) / v;
  return out;
}

void test_convergence_1d(int nx[NUM_GRIDS], int kx, int order) {
  FILE* fp;
  int i,j;
  int ky = 0;
  double y = 0;
  double dx;
  double* x[NUM_GRIDS];
  double* f[NUM_GRIDS];
  double* f_interp[NUM_RESOLUTIONS];
  double* f_error[NUM_RESOLUTIONS];
  printf("Testing 1D convergence\n");
  printf("Using frequency:\n\tk = %d\n",kx);
  printf("Using %d resolutions:\n",NUM_RESOLUTIONS);
  printf("\tnx = [ ");
  for (i = 0; i < NUM_RESOLUTIONS; i++) {
    printf("%d ",nx[i]);
  }
  printf("]\n");
  printf("And interpolating to a fine grid with:\n\tnx = %d\n",
 	 nx[NUM_GRIDS-1]);
  printf("Using convergence order:\n\tORDER = %d\n",order);
  printf("Generating grids\n");
  for (i = 0; i < NUM_GRIDS; ++i) {
    x[i] = (double *)malloc(sizeof(double)*nx[i]);
    f[i] = (double *)malloc(sizeof(double)*nx[i]);
    dx = get_dx_from_bounds(XMIN,XMAX,nx[i]);
    for (int j = 0; j < nx[i]; ++j) {
      x[i][j] = XMIN + j*dx;
      f[i][j] = f_ana(x[i][j],y,kx,ky);
    }
  }
  for (int i = 0; i < NUM_RESOLUTIONS; i++) {
    f_interp[i] = (double*)malloc(sizeof(double)*nx[NUM_GRIDS-1]);
    f_error[i] = (double*)malloc(sizeof(double)*nx[NUM_GRIDS-1]);
  }
  printf("Interpolating and calculating error\n");
  for (int i = 0; i < NUM_RESOLUTIONS; i++) {
    for (int j = 0; j < nx[NUM_GRIDS-1]; ++j) {
      f_interp[i][j] = lagrange_interp_1Dfo(x[NUM_GRIDS-1][j],
					    order,
					    x[i], nx[i],
					    f[i]);
      f_error[i][j] = f_interp[i][j] - f[NUM_GRIDS-1][j];
    }
  }
  printf("Saving to text file\n");
  fp = fopen(FILENAME_SLICE,"w");
  fprintf(fp,"# X\tE(nx =");
  for (i = 0; i < NUM_RESOLUTIONS; i++) {
    fprintf(fp,"\t%d",nx[i]);
  }
  fprintf(fp," )\n");
  for ( j = 0; j < nx[NUM_GRIDS-1]; j++) {
    fprintf(fp,"%lf",x[NUM_GRIDS-1][j]);
    for (i = 0; i < NUM_RESOLUTIONS; ++i) {
      fprintf(fp," %lf",f_error[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  printf("Cleaning up.\n");
  for (i = 0; i < NUM_GRIDS; ++i) {
    free(x[i]); free(f[i]);
  }
  for (i = 0; i < NUM_RESOLUTIONS; ++i) {
    free(f_interp[i]); free(f_error[i]);
  }
}

void test_convergence_2d(int nx[NUM_GRIDS], int ny[NUM_GRIDS],
			 int order[NUM_ORDERS],
			 int kx, int ky) {
  double norm = L2_NORM;//INFINITY_NORM;
  FILE* fp;
  int i,j,k,o,i1d;
  double dx, dy;
  double* x[NUM_GRIDS];
  double* y[NUM_GRIDS];
  double* f[NUM_GRIDS];
  double* f_interp[NUM_RESOLUTIONS];
  double* f_error[NUM_RESOLUTIONS];
  double norm_error[NUM_ORDERS][NUM_RESOLUTIONS];
  printf("Testing 2D convergence\n");
  printf("Using frequencies:\n\tkx = %d\n\tky = %d\n",kx,ky);
  printf("Using %d resolutions:\n",NUM_RESOLUTIONS);
  printf("\tnx = [ ");
  for (i = 0; i < NUM_RESOLUTIONS; i++) {
    printf("%d ",nx[i]);
  }
  printf("]\n");
  printf("\tny = [ ");
  for (i = 0; i < NUM_RESOLUTIONS; i++) {
    printf("%d ",ny[i]);
  }
  printf("]\n");
  printf("Using %d orders\n",NUM_ORDERS);
  printf("\tO = [ ");
  for( i = 0; i < NUM_ORDERS; ++i) {
    printf("%d ",order[i]);
  }
  printf("]\n");
  printf("Generating grids");
  for (i = 0; i < NUM_GRIDS; ++i) {
    x[i] = (double *)malloc(sizeof(double)*nx[i]);
    y[i] = (double *)malloc(sizeof(double)*ny[i]);
    f[i] = (double *)malloc(sizeof(double)*nx[i]*ny[i]);
    dx = get_dx_from_bounds(XMIN,XMAX,nx[i]);
    dy = get_dx_from_bounds(XMIN,XMAX,ny[i]);
    for (int j = 0; j < nx[i]; ++j) {
      x[i][j] = XMIN + j*dx;
      for (int k = 0; k < ny[i]; ++k) {
	y[i][k] = XMIN + k*dy;
	i1d = index_2D_to_1D(j, nx[i], k, ny[i]);
	f[i][i1d] = f_ana(x[i][j],y[i][k],kx,ky);
      }
    }
  }
  for (int i = 0; i < NUM_RESOLUTIONS; i++) {
    f_interp[i] = (double*)malloc(sizeof(double)*nx[NUM_GRIDS-1]*ny[NUM_GRIDS-1]);
    f_error[i] = (double*)malloc(sizeof(double)*nx[NUM_GRIDS-1]*ny[NUM_GRIDS-1]);
  }
  printf("Interpolating and calculating error\n");
  for ( o = 0; o < NUM_ORDERS; ++o) {
    printf("\tCalculating for order %d\n",order[o]);
    for (i = 0; i < NUM_RESOLUTIONS; ++i) {
      printf("\t\tfor resolution [nx,ny] = [%d,%d]\n",nx[i],ny[i]);
      for (j = 0; j < nx[NUM_GRIDS-1]; ++j) {
	for (k = 0; k < ny[NUM_GRIDS-1]; ++k) {
	  i1d = index_2D_to_1D(j, nx[NUM_GRIDS-1],
			       k, ny[NUM_GRIDS-1]);
	  f_interp[i][i1d] = lagrange_interp_2Dfo(x[NUM_GRIDS-1][j],
						  y[NUM_GRIDS-1][k],
						  order[o],order[o],
						  x[i],nx[i],
						  y[i],ny[i],
						  f[i]);
	  f_error[i][i1d] = f_interp[i][i1d] - f[NUM_GRIDS-1][i1d];
	}
      }
      if (norm == L2_NORM) {
	dx = x[i][1] - x[i][0];
	dy = y[i][1] - y[i][0];
	norm_error[o][i] = l2_norm(f_error[i],dx,dy,
				   nx[NUM_GRIDS-1]*ny[NUM_GRIDS-1]);
      } else if (norm == INFINITY_NORM) {
	norm_error[o][i] = infinity_norm(f_error[i],
					 nx[NUM_GRIDS-1]*ny[NUM_GRIDS-1]);
      } else {
	printf("ERROR norm not supported.\n");
	assert (norm == L2_NORM || norm == INFINITY_NORM);
      }
    }
  }
  printf("Saving to file.\n");
  fp = fopen(FILENAME_NORM,"w");
  fprintf(fp,"#nx ny error,order =");
  for (o = 0; o < NUM_ORDERS; ++o) {
    fprintf(fp," %d",order[o]);
  }
  fprintf(fp,"\n");
  for (i = 0; i < NUM_RESOLUTIONS; ++i) {
    fprintf(fp,"%d %d",nx[i],ny[i]);
    for (o = 0; o < NUM_ORDERS; ++o) {
      fprintf(fp," %lf",norm_error[o][i]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  printf("Cleaning up.\n");
  for (i = 0; i < NUM_GRIDS; ++i) {
    free(x[i]); free(y[i]); free(f[i]);
  }
  for (i = 0; i < NUM_RESOLUTIONS; ++i) {
    free(f_interp[i]); free(f_error[i]);
  }
}
			 

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("%s ORDER_1D\n",argv[0]);
    return 1;
  }
  int order_1d = atoi(argv[1]);

  printf("Beginning convergence test\n");
  printf("Interpolating:\n\tsin(kx*x)*cos(ky*y)\non domain:\n\t[-pi, pi]^2\n");
  // 1D test
  //int nx[NUM_GRIDS] = {5, 9, 17, 33};
  int nx[NUM_GRIDS] = {11, 21, 41, 81};
  int ny[NUM_GRIDS] = {9, 17, 33, 81};
  int order[NUM_ORDERS] = {2, 4, 6};
  int kx = 5;
  int ky = 3;
  test_convergence_1d(nx,kx,order_1d);
  test_convergence_2d(nx,ny,order,kx,ky);

  return 0;
}
