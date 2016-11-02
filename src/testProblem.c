#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "ellipsoid/ellipsoid.h"
#include "ellipsoid/ellSolv.h"
#include "sphere/sphere.h"
#include "sphere/sphSolv.h"
#include "constants.h"


void runTest1() {
  
  const int nCharges = 1;

  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 1.5, 1.25, 1.0);

  //grid dimensions
  int nx = 10;
  int ny = 10;
  int nz = 10;
  int nPoints = nx*ny*nz;
  double xa = -10;
  double xb = 10;
  double ya = -10;
  double yb = 10;
  double za = -10;
  double zb = 10;
  double xh = (xb - xa)/(nx-1);
  double yh = (yb - ya)/(ny-1);
  double zh = (zb - za)/(nz-1);


  double charges[nCharges];
  for(int i=0; i<nCharges; ++i)
    charges[i] = 1.0;

  //source charges
  SPoint sPoints[nCharges];
  Point  ePoints[nCharges];
  sPoints[0].x1 = 2.02;
  sPoints[0].x2 = -0.4;
  sPoints[0].x3 = -.463;
  ePoints[0].x1 = sPoints[0].x1;
  ePoints[0].x2 = sPoints[0].x2;
  ePoints[0].x3 = sPoints[0].x3;

  //permitivity constants
  double e1 = 1.0;
  double e2 = 1.0;

  
  printf("entering\n");
  //solution vector
  SPoint *rs = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  Point  *re = (Point*)  malloc(sizeof(Point) *nPoints);
  SPoint *ro = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  double *exact = (double*) malloc(sizeof(double)*nPoints);
  for(int i=0; i<nx; ++i) {
    for(int j=0; j<ny; ++j) {
      for(int k=0; k<nz; ++k) {
	int index = k*ny*nx + j*nx + i;
	rs[index].x1 = xa + i*xh;
	rs[index].x2 = ya + j*yh;
	rs[index].x3 = za + k*zh;
	re[index].x1 = rs[index].x1;
	re[index].x2 = rs[index].x2;
	re[index].x3 = rs[index].x3;
	ro[index].x1 = re[index].x1;
	ro[index].x2 = re[index].x2;
	ro[index].x3 = re[index].x3;

	exact[index] = 0.0;
	for(int num=0; num<nCharges; ++num)
	  exact[index] += charges[num]/(e1*sqrt((sPoints[num].x1 - ro[index].x1)*(sPoints[num].x1 - ro[index].x1) +
					    (sPoints[num].x2 - ro[index].x2)*(sPoints[num].x2 - ro[index].x2) +
					     (sPoints[num].x3 - ro[index].x3)*(sPoints[num].x3 - ro[index].x3)));
      }
    }
  }
  printf("exiting\n");

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nCharges; ++i) {
    cartesianToSpherical(sPoints+i);
    cartesianToEllipsoidal(&e, ePoints+i);
  }

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nPoints; ++i) {
    cartesianToSpherical(rs + i);
    cartesianToEllipsoidal(&e, re + i);
  }




  //set up problem
  SProblem sProblem = {.positions = sPoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .b = 3.0};
  Problem  eProblem = {.positions = ePoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .e = &e};
  
 
  double *ssolution = (double*) malloc(sizeof(double)*nPoints);
  double *esolution = (double*) malloc(sizeof(double)*nPoints);
  calcCoulombSphericalGrid(&sProblem, 14, rs, nPoints, ssolution);
  calcCoulombEllipsoidalGrid(&eProblem, 14, re, nPoints, esolution);

  //add coulomb inside for plotting. Add later to calculation
  for(int i=0; i<nPoints; ++i) {
    for(int k=0; k<nCharges; ++k) {
  //if inside sphere
  //
      //if(
    }
  }

  printf("exact[0] = %15.15f\n", exact[0]);
  //compare errors
  double *sError = (double*) malloc(sizeof(double)*nPoints);
  double *eError = (double*) malloc(sizeof(double)*nPoints);
  for(int i=0; i<nx; ++i) {
    for(int j=0; j<ny; ++j) {
      for(int k=0; k<nz; ++k) {
	int index = k*ny*nx + j*nx + i;
	sError[index] = fabs(ssolution[index] - exact[index]);
	eError[index] = fabs(esolution[index] - exact[index]);
	if(i==1)
	  printf("%d - s: %4.4e e: %4.4e\n", i, sError[index], eError[index]);
      }
    }
  }

}

void runTest2() {
  
  const int nCharges = 1;
  
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);

  //grid dimensions
  int nx = 50;
  int ny = 50;
  int nz = 1;
  int nPoints = nx*ny*nz;
  double xa = -10.1;
  double xb = 10.1;
  double ya = -10.1;
  double yb = 10.1;
  double za = 0.1;
  double zb = 0.1;
  double xh = (xb - xa)/(nx-1);
  double yh = (yb - ya)/(ny-1);
  double zh = (zb - za)/(nz-1); //zh = 0


  double charges[nCharges];
  for(int i=0; i<nCharges; ++i)
    charges[i] = 1.0;

  //source charges
  SPoint sPoints[nCharges];
  Point  ePoints[nCharges];
  sPoints[0].x1 = 1.02;
  sPoints[0].x2 = -0.4;
  sPoints[0].x3 = .105;
  ePoints[0].x1 = sPoints[0].x1;
  ePoints[0].x2 = sPoints[0].x2;
  ePoints[0].x3 = sPoints[0].x3;

  //permitivity constants
  double e1 = 1.0;
  double e2 = 1.0;

  
  printf("entering\n");
  //solution vector
  SPoint *rs = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  Point  *re = (Point*)  malloc(sizeof(Point) *nPoints);
  SPoint *ro = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  double *exact = (double*) malloc(sizeof(double)*nPoints);
  for(int i=0; i<nx; ++i) {
    for(int j=0; j<ny; ++j) {
      for(int k=0; k<nz; ++k) {
	int index = k*ny*nx + j*nx + i;
	rs[index].x1 = xa + i*xh;
	rs[index].x2 = ya + j*yh;
	rs[index].x3 = za + k*zh;
	re[index].x1 = rs[index].x1;
	re[index].x2 = rs[index].x2;
	re[index].x3 = rs[index].x3;
	ro[index].x1 = re[index].x1;
	ro[index].x2 = re[index].x2;
	ro[index].x3 = re[index].x3;

	exact[index] = 0.0;
	for(int num=0; num<nCharges; ++num)
	  exact[index] += charges[num]/(e1*sqrt((sPoints[num].x1 - ro[index].x1)*(sPoints[num].x1 - ro[index].x1) +
					    (sPoints[num].x2 - ro[index].x2)*(sPoints[num].x2 - ro[index].x2) +
					     (sPoints[num].x3 - ro[index].x3)*(sPoints[num].x3 - ro[index].x3)));
      }
    }
  }
  printf("exiting\n");

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nCharges; ++i) {
    cartesianToSpherical(sPoints+i);
    cartesianToEllipsoidal(&e, ePoints+i);
  }

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nPoints; ++i) {
    cartesianToSpherical(rs + i);
    cartesianToEllipsoidal(&e, re + i);
  }




  //set up problem
  SProblem sProblem = {.positions = sPoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .b = 3.0};
  Problem  eProblem = {.positions = ePoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .e = &e};
  
 
  double *ssolution = (double*) malloc(sizeof(double)*nPoints);
  double *esolution = (double*) malloc(sizeof(double)*nPoints);
  calcCoulombSphericalGrid(&sProblem, 14, rs, nPoints, ssolution);
  calcCoulombEllipsoidalGrid(&eProblem, 14, re, nPoints, esolution);

  //add coulomb inside for plotting. Add later to calculation
  for(int i=0; i<nPoints; ++i) {
    for(int k=0; k<nCharges; ++k) {
  //if inside sphere
  //
      //if(
    }
  }

  printf("exact[0] = %15.15f\n", exact[0]);
  //compare errors
  double *sError = (double*) malloc(sizeof(double)*nPoints);
  double *eError = (double*) malloc(sizeof(double)*nPoints);
  //for(int i=0; i<nx; ++i) {
  //for(int j=0; j<ny; ++j) {
  //for(int k=0; k<nz; ++k) {
  //int index = k*ny*nx + j*nx + i;
  for(int i=0; i<nPoints; ++i) {
    if( (ro[i].x1*ro[i].x1)/(a*a) + (ro[i].x2*ro[i].x2)/(b*b) + (ro[i].x3*ro[i].x3)/(c*c) < 1) {
      esolution[i] = 0.0;
      ssolution[i] = 0.0;
    }
    sError[i] = fabs(ssolution[i] - exact[i]);
    eError[i] = fabs(esolution[i] - exact[i]);
    printf("x: %2.2f y: %2.2f - s: %4.4e e: %4.4e exact: %5.5f\n", ro[i].x1, ro[i].x2, sError[i], eError[i], exact[i]);
  }
  //if(j==0)
  //printf("%d - s: %4.4e e: %4.4e\n", index, sError[index], eError[index]);
  //      }
  //}
  //}


  FILE *fp  = fopen("plotelltest2.txt", "w");
  FILE *fp2 = fopen("plotsphtest2.txt", "w");
  for(int j=0; j<ny; ++j) {
    for(int i=0; i<nx; ++i) {
      int index = j*nx + i;
      fprintf(fp,  "%15.15f ", esolution[index]);
      fprintf(fp2, "%15.15f ", ssolution[index]);
    }
    fprintf(fp,  "\n");
    fprintf(fp2, "\n");
  }

  fclose(fp);
  fclose(fp2);
}

void runTest3() {
  
  const int nCharges = 1;
  
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);

  //grid dimensions
  int nx = 50;
  int ny = 50;
  int nz = 1;
  int nPoints = nx*ny*nz;
  double xa = -10.1;
  double xb = 10.1;
  double ya = -10.1;
  double yb = 10.1;
  double za = 0.1;
  double zb = 0.1;
  double xh = (xb - xa)/(nx-1);
  double yh = (yb - ya)/(ny-1);
  double zh = (zb - za)/(nz-1); //zh = 0


  double charges[nCharges];
  for(int i=0; i<nCharges; ++i)
    charges[i] = 1.0;

  //source charges
  SPoint sPoints[nCharges];
  Point  ePoints[nCharges];
  sPoints[0].x1 = 2.02;
  sPoints[0].x2 = -0.8;
  sPoints[0].x3 = .105;
  ePoints[0].x1 = sPoints[0].x1;
  ePoints[0].x2 = sPoints[0].x2;
  ePoints[0].x3 = sPoints[0].x3;

  //permitivity constants
  double e1 = 1.0;
  double e2 = 2.0;

  
  printf("entering\n");
  //solution vector
  SPoint *rs = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  Point  *re = (Point*)  malloc(sizeof(Point) *nPoints);
  SPoint *ro = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  double *exact = (double*) malloc(sizeof(double)*nPoints);
  for(int i=0; i<nx; ++i) {
    for(int j=0; j<ny; ++j) {
      for(int k=0; k<nz; ++k) {
	int index = k*ny*nx + j*nx + i;
	rs[index].x1 = xa + i*xh;
	rs[index].x2 = ya + j*yh;
	rs[index].x3 = za + k*zh;
	re[index].x1 = rs[index].x1;
	re[index].x2 = rs[index].x2;
	re[index].x3 = rs[index].x3;
	ro[index].x1 = re[index].x1;
	ro[index].x2 = re[index].x2;
	ro[index].x3 = re[index].x3;

	exact[index] = 0.0;
	for(int num=0; num<nCharges; ++num)
	  exact[index] += charges[num]/(e1*sqrt((sPoints[num].x1 - ro[index].x1)*(sPoints[num].x1 - ro[index].x1) +
					    (sPoints[num].x2 - ro[index].x2)*(sPoints[num].x2 - ro[index].x2) +
					     (sPoints[num].x3 - ro[index].x3)*(sPoints[num].x3 - ro[index].x3)));
      }
    }
  }
  printf("exiting\n");

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nCharges; ++i) {
    cartesianToSpherical(sPoints+i);
    cartesianToEllipsoidal(&e, ePoints+i);
  }

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nPoints; ++i) {
    cartesianToSpherical(rs + i);
    cartesianToEllipsoidal(&e, re + i);
  }




  //set up problem
  SProblem sProblem = {.positions = sPoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .b = e.a};
  Problem  eProblem = {.positions = ePoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .e = &e};
  
 
  double *ssolution = (double*) malloc(sizeof(double)*nPoints);
  double *esolution = (double*) malloc(sizeof(double)*nPoints);
  calcCoulombSphericalGrid(&sProblem, 4, rs, nPoints, ssolution);
  calcCoulombEllipsoidalGrid(&eProblem, 4, re, nPoints, esolution);

  //add coulomb inside for plotting. Add later to calculation
  for(int i=0; i<nPoints; ++i) {
    for(int k=0; k<nCharges; ++k) {
  //if inside sphere
  //
      //if(
    }
  }

  printf("exact[0] = %15.15f\n", exact[0]);
  //compare errors
  double *sError = (double*) malloc(sizeof(double)*nPoints);
  double *eError = (double*) malloc(sizeof(double)*nPoints);
  //for(int i=0; i<nx; ++i) {
  //for(int j=0; j<ny; ++j) {
  //for(int k=0; k<nz; ++k) {
  //int index = k*ny*nx + j*nx + i;
  for(int i=0; i<nPoints; ++i) {
    if( (ro[i].x1*ro[i].x1)/(a*a) + (ro[i].x2*ro[i].x2)/(b*b) + (ro[i].x3*ro[i].x3)/(c*c) < 1) {
      esolution[i] += exact[i];
    }
    if( ro[i].x1*ro[i].x1 + ro[i].x2*ro[i].x2 + ro[i].x3*ro[i].x3 < a*a) {
      ssolution[i] += exact[i];
    }
    sError[i] = fabs(ssolution[i] - exact[i]);
    eError[i] = fabs(esolution[i] - exact[i]);
    printf("x: %2.2f y: %2.2f - s: %4.4e e: %4.4e exact: %5.5f\n", ro[i].x1, ro[i].x2, sError[i], eError[i], exact[i]);
  }
  //if(j==0)
  //printf("%d - s: %4.4e e: %4.4e\n", index, sError[index], eError[index]);
  //      }
  //}
  //}


  FILE *fp  = fopen("plotelltest3.txt", "w");
  FILE *fp2 = fopen("plotsphtest3.txt", "w");

  for(int j=0; j<ny; ++j) {
    for(int i=0; i<nx; ++i) {
      int index = j*nx + i;
      fprintf(fp, "%15.15f ", esolution[index]);
      fprintf(fp, "%15.15f ", ssolution[index]);
    }
    fprintf(fp,  "\n");
    fprintf(fp2, "\n");
  }

  fclose(fp);
  fclose(fp2);
}

void runTest4() {
  
  const int nCharges = 2;
  
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);

  //grid dimensions
  int nx = 50;
  int ny = 50;
  int nz = 1;
  int nPoints = nx*ny*nz;
  double xa = -10.1;
  double xb = 10.1;
  double ya = -10.1;
  double yb = 10.1;
  double za = 0.1;
  double zb = 0.1;
  double xh = (xb - xa)/(nx-1);
  double yh = (yb - ya)/(ny-1);
  double zh = (zb - za)/(nz-1); //zh=0

  double charges[nCharges];
  for(int i=0; i<nCharges; ++i)
    charges[i] = 1.0;

  //source charges
  SPoint sPoints[nCharges];
  Point  ePoints[nCharges];
  sPoints[0].x1 = 2.02;
  sPoints[0].x2 = -0.8;
  sPoints[0].x3 = .105;
  sPoints[1].x1 = -2.2;
  sPoints[1].x2 = .3;
  sPoints[1].x3 = -.08;
  ePoints[0].x1 = sPoints[0].x1;
  ePoints[0].x2 = sPoints[0].x2;
  ePoints[0].x3 = sPoints[0].x3;
  ePoints[1].x1 = sPoints[1].x1;
  ePoints[1].x2 = sPoints[1].x2;
  ePoints[1].x3 = sPoints[1].x3;

  //permitivity constants
  double e1 = 1.0;
  double e2 = 2.0;

  
  printf("entering\n");
  //solution vector
  SPoint *rs = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  Point  *re = (Point*)  malloc(sizeof(Point) *nPoints);
  SPoint *ro = (SPoint*) malloc(sizeof(SPoint)*nPoints);
  double *exact = (double*) malloc(sizeof(double)*nPoints);
  for(int i=0; i<nx; ++i) {
    for(int j=0; j<ny; ++j) {
      for(int k=0; k<nz; ++k) {
	int index = k*ny*nx + j*nx + i;
	rs[index].x1 = xa + i*xh;
	rs[index].x2 = ya + j*yh;
	rs[index].x3 = za + k*zh;
	re[index].x1 = rs[index].x1;
	re[index].x2 = rs[index].x2;
	re[index].x3 = rs[index].x3;
	ro[index].x1 = re[index].x1;
	ro[index].x2 = re[index].x2;
	ro[index].x3 = re[index].x3;

	exact[index] = 0.0;
	for(int num=0; num<nCharges; ++num)
	  exact[index] += charges[num]/(e1*sqrt((sPoints[num].x1 - ro[index].x1)*(sPoints[num].x1 - ro[index].x1) +
					    (sPoints[num].x2 - ro[index].x2)*(sPoints[num].x2 - ro[index].x2) +
					     (sPoints[num].x3 - ro[index].x3)*(sPoints[num].x3 - ro[index].x3)));
      }
    }
  }
  printf("exiting\n");

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nCharges; ++i) {
    cartesianToSpherical(sPoints+i);
    cartesianToEllipsoidal(&e, ePoints+i);
  }

  //convert to spherical/ellipsoidal coordinates
  for(int i=0; i<nPoints; ++i) {
    cartesianToSpherical(rs + i);
    cartesianToEllipsoidal(&e, re + i);
  }




  //set up problem
  SProblem sProblem = {.positions = sPoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .b = e.a};
  Problem  eProblem = {.positions = ePoints, .charges = charges, .nCharges = nCharges, .e1 = e1, .e2 = e2, .e = &e};
  
 
  double *ssolution = (double*) malloc(sizeof(double)*nPoints);
  double *esolution = (double*) malloc(sizeof(double)*nPoints);
  calcCoulombSphericalGrid(&sProblem, 4, rs, nPoints, ssolution);
  calcCoulombEllipsoidalGrid(&eProblem, 4, re, nPoints, esolution);

  //add coulomb inside for plotting. Add later to calculation
  for(int i=0; i<nPoints; ++i) {
    for(int k=0; k<nCharges; ++k) {
  //if inside sphere
  //
      //if(
    }
  }

  printf("exact[0] = %15.15f\n", exact[0]);
  //compare errors
  double *sError = (double*) malloc(sizeof(double)*nPoints);
  double *eError = (double*) malloc(sizeof(double)*nPoints);
  //for(int i=0; i<nx; ++i) {
  //for(int j=0; j<ny; ++j) {
  //for(int k=0; k<nz; ++k) {
  //int index = k*ny*nx + j*nx + i;
  for(int i=0; i<nPoints; ++i) {
    if( (ro[i].x1*ro[i].x1)/(a*a) + (ro[i].x2*ro[i].x2)/(b*b) + (ro[i].x3*ro[i].x3)/(c*c) < 1) {
      //esolution[i] += exact[i];
    }
    if( ro[i].x1*ro[i].x1 + ro[i].x2*ro[i].x2 + ro[i].x3*ro[i].x3 < a*a) {
      //ssolution[i] += exact[i];
    }
    sError[i] = fabs(ssolution[i] - exact[i]);
    eError[i] = fabs(esolution[i] - exact[i]);
    printf("x: %2.2f y: %2.2f - s: %4.4e e: %4.4e exact: %5.5f\n", ro[i].x1, ro[i].x2, sError[i], eError[i], exact[i]);
  }
  //if(j==0)
  //printf("%d - s: %4.4e e: %4.4e\n", index, sError[index], eError[index]);
  //      }
  //}
  //}


  FILE *fp  = fopen("plotelltest4.txt", "w");
  FILE *fp2 = fopen("plotsphtest4.txt", "w");

  for(int j=0; j<ny; ++j) {
    for(int i=0; i<nx; ++i) {
      int index = j*nx + i;
      fprintf(fp,  "%15.15f ", esolution[index]);
      fprintf(fp2, "%15.15f ", ssolution[index]);
    }
    fprintf(fp,  "\n");
    fprintf(fp2, "\n");
  }

  fclose(fp);
  fclose(fp2);
}

void test1sphere() {

  const int npoints = 1;

  //calculation point
  SPoint r = {.x1 = -3.0, .x2 = -1.6, .x3 = 0.2, .type = 'c'};

  //source charges
  SPoint points[npoints];
  points[0].x1 = 2.0;
  points[0].x2 = 0.3;
  points[0].x3 = -0.1;
  
  double charges[npoints];
  charges[0] = 1.0;

  double exact = 0;
  for(int i=0; i<npoints; i++)
    exact += charges[i]/sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

  printf("the exact result is %15.15f\n", exact);


  cartesianToSpherical(&r);
  for(int i=0; i<npoints; ++i)
    cartesianToSpherical(points+i);

  

  printf("r.x1 = %15.15f\n", r.x1);
  printf("r.x2 = %15.15f\n", r.x2);
  printf("r.x3 = %15.15f\n", r.x3);


  for(int i=0; i<npoints; ++i) {
    printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }

  



  
  //set up problem
  SProblem problem = {.positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0, .b = 1};


  calcCoulombSpherical(&problem, 50, &r);
  
  printf("the exact result is %15.15f\n", exact);
  //speedTester(&problem, 10, &r);
  //printf("norm for %d,%d: %8.8e\n", 25, 39, calcNormalization(&e, 25, 39));
}

void test2() {
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  

  const int npoints = 1;
  
  //calculation point
  Point r = {.x1 = -1.1, .x2 = 1.5, .x3 = 4.1, .type = 'c'};
  
  //source charges
  Point points[npoints];
  points[0].x1 = 0.2;
  points[0].x2 = 0.3;
  points[0].x3 = -0.5;

  
  double charges[npoints];
  charges[0] = 1.0;
  //charges[1] = 1.0;


  double exact = 0;
  for(int i=0; i<npoints; i++)
    exact += charges[i]/sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

  printf("the exact result is %15.15f\n", exact);


  cartesianToEllipsoidal(&e, &r);
  //cartesianToEllipsoidal(&e, &r);
  for(int i=0; i<npoints; ++i) {
    cartesianToEllipsoidal(&e, points+i);
    //cartesianToEllipsoidal(&e, points+i);
  }

  printf("r.x1 = %15.15f\n", r.x1);
  printf("r.x2 = %15.15f\n", r.x2);
  printf("r.x3 = %15.15f\n", r.x3);

  for(int i=0; i<npoints; ++i) {
    printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }


  //set up problem
  Problem problem = {.positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0, .e = &e};


  
  double approx = calcCoulomb(&problem, 15, &r);

  printf("the error is %15.15f\n", fabs(exact-approx));
  //speedTester(&problem, 10, &r);
  //printf("norm for %d,%d: %8.8e\n", 25, 39, calcNormalization(&e, 25, 39));
}

void test2sphere() {

  const int npoints = 2;

  //calculation point
  SPoint r = {.x1 = 1.2, .x2 = -0.5, .x3 = 3.0, .type = 'c'};

  //source charges
  SPoint points[npoints];
  points[0].x1 = -0.2;
  points[0].x2 = 1.0;
  points[0].x3 = -.8;
  points[1].x1 = 0.3;
  points[1].x2 = -1.2;
  points[1].x3 = 0.5;
  
  double charges[npoints];
  charges[0] = 1.0;
  charges[1] = -1.0;

  double exact = 0;
  for(int i=0; i<npoints; i++)
    exact += charges[i]/sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

  printf("the exact result is %15.15f\n", exact);
  
  cartesianToSpherical(&r);
  for(int i=0; i<npoints; ++i)
    cartesianToSpherical(points+i);

  printf("r.x1 = %15.15f\n", r.x1);
  printf("r.x2 = %15.15f\n", r.x2);
  printf("r.x3 = %15.15f\n", r.x3);

  for(int i=0; i<npoints; ++i) {
    printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }

  
  //set up problem
  SProblem problem = {.positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0};

  calcCoulombSpherical(&problem, 300, &r);
  
  //calcCoulombSpherical(&problem, 30, &r);

}

void test3() {
  EllipsoidalSystem e;
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  initEllipsoidalSystem(&e, a, b, c);
  

  const int npoints = 4;

  //calculation point
  Point r = {.x1 = 1.0, .x2 = -1.7, .x3 = 4.0, .type = 'c'};

  //source charges
  Point points[npoints*npoints];
  //generate points on the surface of a 2.9, 1.9, .9 ellipse
  double distIn = .2;
  double pointsA = a-distIn;
  double pointsB = b-distIn;
  double pointsC = c-distIn;
  double delU = (2*M_PI-.01)/npoints;
  double delV = (M_PI-.01)/npoints;
  double u,v;
  for(int i=0; i<npoints; ++i) {
    u = (-M_PI/2+.01) + i*delU;
    for(int j=0; j<npoints; ++j) {
      v = (-M_PI+.01) + j*delV;
      points[i*npoints+j].x1 = pointsA*cos(u)*sin(v);
      points[i*npoints+j].x2 = pointsB*sin(u)*sin(v);
      points[i*npoints+j].x3 = pointsC*cos(v);
      points[i*npoints+j].type = 'c';
    }
  }

  //charges array, all set to 1
  double charges[npoints*npoints];
  for(int i=0; i<npoints*npoints; ++i)
    charges[i] = 1.0;

  //calculate exact solution ahead of time
  double exact = 0;
  for(int i=0; i<npoints*npoints; i++)
    exact += charges[i]/sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

  printf("the exact result is %15.15f\n", exact);


  //convert all points to ellipsoidal coordinates
  cartesianToEllipsoidal(&e, &r);
  for(int i=0; i<npoints*npoints; ++i)
    cartesianToEllipsoidal(&e, points+i);

  for(int i=0; i<npoints*npoints; ++i) {
    //printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    //printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    //printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }


  //set up problem
  Problem problem = {.positions = points, .charges = charges, .nCharges = npoints*npoints, .e1 = 1.0, .e2 = 1.0, .e = &e};

  clock_t start, diff;
  int msec;
  //calculate the coulomb potential
  start = clock();  
  double ellSol = calcCoulomb(&problem, 32, &r); ellSol *= 1.0;
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("The time for the solution to n=10 is %d seconds and %d milliseconds\n", msec/1000, msec%1000);
}

void test3sphere() {

  double a = 2.0;
  double b = 2.5;
  double c = 3.0;

  const int npoints = 11;

  //calculation points
  SPoint r = {.x1 = -2.0, .x2 = 3.5, .x3 = 4.7, .type = 'c'};

  //source charges
  SPoint points[npoints*npoints];
  //generate points on the surface of a 2.9, 1.9, .9 ellipse
  double distIn = .1;
  double pointsA = a-distIn;
  double pointsB = b-distIn;
  double pointsC = c-distIn;
  double delU = (M_PI-.01)/npoints;
  double delV = (2*M_PI-.01)/npoints;
  double u,v;
  for(int i=0; i<npoints; ++i) {
    u = (-M_PI/2+.01) + i*delU;
    for(int j=0; j<npoints; ++j) {
      v = (-M_PI+.01) + j*delV;
      points[i*npoints+j].x1 = pointsA*cos(u)*cos(v);
      points[i*npoints+j].x2 = pointsB*cos(u)*sin(v);
      points[i*npoints+j].x3 = pointsC*sin(u);
      points[i*npoints+j].type = 'c';
    }
  }


  //charges array, all set to 1
  double charges[npoints*npoints];
  for(int i=0; i<npoints*npoints; ++i)
    charges[i] = 1.0;

  //calculate exact solution ahead of time
  double exact = 0;
  for(int i=0; i<npoints*npoints; i++)
    exact += charges[i]/sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));
  
  printf("the exact result is %15.15f\n", exact);
  
  //convert everything to spherical coordinates
  cartesianToSpherical(&r);
  for(int i=0; i<npoints*npoints; ++i)
    cartesianToSpherical(points+i);


  printf("r.x1 = %15.15f\n", r.x1);
  printf("r.x2 = %15.15f\n", r.x2);
  printf("r.x3 = %15.15f\n", r.x3);
  /*
  for(int i=0; i<npoints*npoints; ++i) {
    cartesianToSpherical(points+i);
    printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }
  */


  //set up problem
  SProblem problem = {.positions = points, .charges = charges, .nCharges = npoints*npoints, .e1 = 1.0, .e2 = 1.0};

  double *solution = (double*) malloc(sizeof(double)*1);

  //calculate the coulomb potential
  calcCoulombSphericalGrid(&problem, 50, &r, 1, solution);
  printf("the exact result is %15.15f\n", exact);
  
}

void test4sphere()
{
  
  double a = 3.0;
  double b = 2.0;
  //double c = 1.0;

  double gridXMin = -2.0*a;
  double gridXMax = fabs(gridXMin);
  double gridYMin = -2.0*a;
  double gridYMax = fabs(gridYMin);

  const double srad = a;
  const int nCharges = 7;
  const double e1 = 2.5;
  const double e2 = 80.1;
  const double distIn = .1;
  const int xdim = 200;
  const int ydim = 200;
  const int dim  = xdim*ydim;

  //make grid of calculation points
  SPoint r = {.x1 = -2.0, .x2 = 3.5, .x3 = 0.0, .type = 'c'};
  SPoint *grid = (SPoint*) malloc(sizeof(SPoint)*dim);
  double delX = (gridXMax - gridXMin)/(xdim - 1.0);
  double delY = (gridYMax - gridYMin)/(ydim - 1.0);
  double xVal, yVal;
  for(int i=0; i<xdim; ++i) {
    xVal = gridXMin + i*delX;
    for(int j=0; j<ydim; ++j) {
      yVal = gridYMin + j*delY;
      grid[i*ydim + j].x1 = xVal;
      grid[i*ydim + j].x2 = yVal;
      grid[i*ydim + j].x3 = 0.0;
    }
  }

  //source charges
  SPoint chargePoints[nCharges];

  //generate points on the surface of a 2.9, 1.9, .9 ellipse
  double pointsA = a - distIn;
  double pointsB = b - distIn;
  double pointsC = 0.0;
  double delV = (2*M_PI-.01)/nCharges;
  double v;
  double pert;
  for(int k=0; k<nCharges; ++k) {
      v = .01 + k*delV;
      pert = ((double)rand()/(double)RAND_MAX);
      pert = 1.0;//pert = 0.5;//0.5*pert + 0.5;
      chargePoints[k].x1 = pert*pointsA*cos(v);
      chargePoints[k].x2 = pert*pointsB*sin(v);
      chargePoints[k].x3 = pert*pointsC;
      chargePoints[k].type = 'c';
  }

  //charges array, all set to 1
  double chargeValues[nCharges];
  for(int k=0; k<nCharges; ++k)
    chargeValues[k] = 1.0;
  
  double *solution = (double*) malloc(sizeof(double)*dim);
  double *coulombInside = (double*) malloc(sizeof(double)*dim);

  for(int k=0; k<dim; ++k) {
    coulombInside[k] = 0.0;
    if(sqrt(grid[k].x1*grid[k].x1 + grid[k].x2*grid[k].x2 + grid[k].x3*grid[k].x3) < srad) {
      for(int i=0; i<nCharges; ++i) {
	coulombInside[k] += chargeValues[i]/sqrt((chargePoints[i].x1-grid[k].x1)*(chargePoints[i].x1-grid[k].x1) + (chargePoints[i].x2-grid[k].x2)*(chargePoints[i].x2-grid[k].x2) + (chargePoints[i].x3-grid[k].x3)*(chargePoints[i].x3-grid[k].x3));
      }
      coulombInside[k] *= 1./e1;
    }
  }
  
  //calculate exact solution ahead of time
  double exact = 0;
  for(int k=0; k<nCharges; k++)
    exact += chargeValues[k]/1.0;
  
  printf("the exact result is %15.15f\n", exact);

  printf("grid[0].x1 = %15.15f\n", grid[0].x1);
  printf("grid[0].x2 = %15.15f\n", grid[0].x2);
  printf("grid[0].x3 = %15.15f\n", grid[0].x3);
  //convert everything to spherical coordinates
  cartesianToSpherical(&r);
  for(int i=0; i<nCharges; ++i)
    cartesianToSpherical(chargePoints+i);
  for(int i=0; i<dim; ++i)
    cartesianToSpherical(grid+i);


  printf("grid[0].x1 = %15.15f\n", grid[0].x1);
  printf("grid[0].x2 = %15.15f\n", grid[0].x2);
  printf("grid[0].x3 = %15.15f\n", grid[0].x3);
  /*
  for(int i=0; i<nCharges; ++i) {
    cartesianToSpherical(chargePoints+i);
    printf("chargePoints[%d].x1 = %15.15f\n", i, chargePoints[i].x1);
    printf("chargePoints[%d].x2 = %15.15f\n", i, chargePoints[i].x2);
    printf("chargePoints[%d].x3 = %15.15f\n", i, chargePoints[i].x3);
  }
  */


  //set up problem
  SProblem problem = {.positions = chargePoints, .charges = chargeValues, .nCharges = nCharges, .e1 = e1, .e2 = e2, .b = srad};


  //calculate the coulomb potential
  calcCoulombSphericalGrid(&problem, 50, grid, dim, solution);
  printf("pisser\n");
  printf("the exact result is %15.15f\n", solution[0]);
  
  //OUTPUT TO FILE contourplot.txt
  FILE *fp = fopen("contourplot.txt", "w");
  int k;
  double outSol;
  double cutOff = 2.0;
  for(int i=0; i<xdim; ++i) {
    for(int j=0; j<ydim; ++j) {
      k = i*ydim + j;
      outSol = solution[k];
      if(fabs(outSol) > cutOff)
	outSol = cutOff*(outSol/fabs(outSol));
      //if(solution[k] < 0)
      //solution[k] = 1e1;
      //fprintf(fp, "%8.8f ", solution[k]);
      fprintf(fp, "%16.16e", outSol);
      if(j != ydim-1)
	fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

#undef __FUNCT__
#define __FUNCT__ "RunArg"
PetscErrorCode RunArg()
{
  PetscErrorCode ierr;
  char buff[64];
  int npts;
  double a, b, c;
  double *xyz;
  double *chargeValues;
  double *solution;
  double freeEnergy;

  PetscFunctionBegin;
  
  //constants
  const double q     = ELECTRON_CHARGE;
  const double Na    = AVOGADRO_NUMBER;
  const double JperC = 4.184; /* Jouled/Calorie */
  const double cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/M_PI; /* kcal ang/mol */


  FILE *fp = fopen("../../pointbem/src/ellData.txt", "r");
  // READ NUMBER OF POINTS
  fscanf(fp, "%s", buff); //skip "size: "
  fscanf(fp, "%d", &npts);

  //initialize local charges/coordinates arrays
  xyz     = (double*) malloc(sizeof(double)*3*npts);
  chargeValues = (double*) malloc(sizeof(double)*npts);
  // READ A,B,C ELLIPSOID VALUES
  fscanf(fp, "%s", buff); //skip "a: "
  fscanf(fp, "%lf", &a); //read a value
  fscanf(fp, "%s", buff); //skip "b: "
  fscanf(fp, "%lf", &b); //read b value
  fscanf(fp, "%s", buff); //skip "c: "
  fscanf(fp, "%lf", &c);
  // READ CHARGE X,Y,Z AND CHARGE DATA
  for(int i=0; i<npts; ++i) {
    fscanf(fp, "%lf", xyz+3*i+0); //read x coordinate
    fscanf(fp, "%lf", xyz+3*i+1); //read y coordinate
    fscanf(fp, "%lf", xyz+3*i+2); //read z coordinate
    fscanf(fp, "%lf", chargeValues+i); //read charge
    printf("%lf %lf %lf %lf\n", xyz[3*i+0], xyz[3*i+1], xyz[3*i+2], chargeValues[i]);
  }
  fclose(fp);


  //initialize ellipsoidal system
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);
  
  //set charge positions and calculation points to be the same
  Point *solPoints = (Point*) malloc(sizeof(Point)*npts);
  Point *chargePoints = (Point*) malloc(sizeof(Point)*npts);
  for(int i=0; i<npts; ++i) {
    solPoints[i].x1 = xyz[3*i+0];
    solPoints[i].x2 = xyz[3*i+1];
    solPoints[i].x3 = xyz[3*i+2];
    solPoints[i].type = 'c';
    chargePoints[i] = solPoints[i];
  }

  //convert to ellipsoidal coordinates
  for(int i=0; i<npts; ++i) {
    cartesianToEllipsoidal(&e, solPoints+i);
    cartesianToEllipsoidal(&e, chargePoints+i);
  }
  
  //initialize ellipsoidal problem context
  Problem prob;
  prob.e = &e;
  prob.positions = chargePoints;
  prob.charges = chargeValues;
  prob.nCharges = npts;
  prob.e1 = 4.0;
  prob.e2 = 80.0;
  
  //solution vector
  solution = (double*) malloc(sizeof(double)*npts);
  
  //calculate potential on solPoints
  calcCoulombEllipsoidalGrid(&prob, 10, solPoints, npts, solution);
  printf("wow!\n");
  for(int i=0; i<npts; ++i)
    printf("solution[%d] = %15.15f\n", i, solution[i]);
  for(int i=0; i<npts; ++i)
    printf("chargeValues[%d] = %15.15f\n", i, chargeValues[i]);
  
  //calculate free energy
  freeEnergy = 0;
  for(int i=0; i<npts; ++i) {
    freeEnergy += solution[i]*chargeValues[i];
    printf("freeEnergy[%d] = %15.15f\n", i, freeEnergy*.5*cf);
  }

  
  free(xyz);
  free(chargeValues);
  free(solPoints);
  free(solution);
  free(chargePoints);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RunArgTester"
PetscErrorCode RunArgTester()
{
  PetscErrorCode ierr;
  char buff[64];
  int npts;
  double a, b, c;
  double *xyz;
  double *chargeValues;
  double *solution;
  double freeEnergy;

  PetscFunctionBegin;
  
  //constants
  const double q     = ELECTRON_CHARGE;
  const double Na    = AVOGADRO_NUMBER;
  const double JperC = 4.184; /* Jouled/Calorie */
  const double cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/M_PI; /* kcal ang/mol */


  FILE *fp = fopen("../../pointbem/src/ellData.txt", "r");
  // READ NUMBER OF POINTS
  fscanf(fp, "%s", buff); //skip "size: "
  fscanf(fp, "%d", &npts);

  //initialize local charges/coordinates arrays
  xyz     = (double*) malloc(sizeof(double)*3*npts);
  chargeValues = (double*) malloc(sizeof(double)*npts);
  // READ A,B,C ELLIPSOID VALUES
  fscanf(fp, "%s", buff); //skip "a: "
  fscanf(fp, "%lf", &a); //read a value
  fscanf(fp, "%s", buff); //skip "b: "
  fscanf(fp, "%lf", &b); //read b value
  fscanf(fp, "%s", buff); //skip "c: "
  fscanf(fp, "%lf", &c);
  // READ CHARGE X,Y,Z AND CHARGE DATA
  for(int i=0; i<npts; ++i) {
    fscanf(fp, "%lf", xyz+3*i+0); //read x coordinate
    fscanf(fp, "%lf", xyz+3*i+1); //read y coordinate
    fscanf(fp, "%lf", xyz+3*i+2); //read z coordinate
    fscanf(fp, "%lf", chargeValues+i); //read charge
    printf("%lf %lf %lf %lf\n", xyz[3*i+0], xyz[3*i+1], xyz[3*i+2], chargeValues[i]);
  }
  fclose(fp);


  //initialize ellipsoidal system
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);
  
  //set charge positions and calculation points to be the same
  Point *solPoints = (Point*) malloc(sizeof(Point)*npts);
  Point *chargePoints = (Point*) malloc(sizeof(Point)*npts);
  for(int i=0; i<npts; ++i) {
    solPoints[i].x1 = xyz[3*i+0];
    solPoints[i].x2 = xyz[3*i+1];
    solPoints[i].x3 = xyz[3*i+2];
    solPoints[i].type = 'c';
    chargePoints[i] = solPoints[i];
  }

  Vec valsbef; //testing
  PetscScalar* valsbefArray; //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, npts, &valsbef); CHKERRQ(ierr); //testing
  ierr = VecGetArray(valsbef, &valsbefArray); CHKERRQ(ierr); //testing
  //convert to ellipsoidal coordinates
  for(int i=0; i<npts; ++i) {
    cartesianToEllipsoidal(&e, solPoints+i);
    ierr = CalcSolidInteriorHarmonic(&e, solPoints[i].x1, solPoints[i].x2, solPoints[i].x3, 2, 3, valsbefArray+i); //testing
    cartesianToEllipsoidal(&e, chargePoints+i);
  }
  ierr = VecRestoreArray(valsbef, &valsbefArray); CHKERRQ(ierr); //testing
  Vec ellPoints, vals; //testing
  PetscScalar* ellPointsArray; //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*npts, &ellPoints); CHKERRQ(ierr); //testing
  ierr = VecGetArray(ellPoints, &ellPointsArray); CHKERRQ(ierr); //testing
  for(PetscInt i=0; i<npts; ++i) { //testing
    ellPointsArray[3*i+0] = solPoints[i].x1; //testing
    ellPointsArray[3*i+1] = solPoints[i].x2; //testing
    ellPointsArray[3*i+2] = solPoints[i].x3; //testing 
  }
  ierr = VecRestoreArray(ellPoints, &ellPointsArray); CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, npts, &vals); CHKERRQ(ierr); //testing
  ierr = CalcSolidInteriorHarmonicVec(&e, npts, ellPoints, 2, 3, vals); //testing
  PetscReal diff, val1, val2; //testing
  diff = 0;
  for(PetscInt i=0; i<npts; ++i) { //testing
    ierr = VecGetValues(vals, 1, &i, &val1); CHKERRQ(ierr); //testing
    ierr = VecGetValues(valsbef, 1, &i, &val2); CHKERRQ(ierr); //testing
    diff += PetscSqrtReal(PetscAbsReal(val1 - val2)); //testing
  }
  printf("Solid Harmonic Error: %3.3e\n", diff); //testing

  //initialize ellipsoidal problem context
  Problem prob;
  prob.e = &e;
  prob.positions = chargePoints;
  prob.charges = chargeValues;
  prob.nCharges = npts;
  prob.e1 = 4.0;
  prob.e2 = 80.0;

  
  /* ---------------------------- */
  /* NOW WE TEST THE Gnp FUNCTION */
  /* ---------------------------- */
  PetscInt NMAX = 5;
  PetscInt count = 0;

  for(PetscInt n=0; n<=NMAX; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p)
      count++;
  }
  Vec srcCharges, GnpVals, GnpValsOld, BnpVals, BnpValsOld, CnpVals, CnpValsOld;
  Vec errorV;
  ierr = VecCreateSeq(PETSC_COMM_SELF, npts, &srcCharges);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &GnpValsOld);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &BnpVals);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &BnpValsOld);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &CnpVals);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &CnpValsOld);CHKERRQ(ierr); //testing
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, &errorV);CHKERRQ(ierr); //testing
  for(PetscInt i=0; i<npts; ++i)
    ierr = VecSetValues(srcCharges, 1, &i, chargeValues+i, INSERT_VALUES);CHKERRQ(ierr); //testing
  ierr = VecAssemblyBegin(srcCharges);CHKERRQ(ierr); ierr = VecAssemblyEnd(srcCharges);CHKERRQ(ierr); //testing

  //WOULD ACTUALLY USE CHARGE POINTS, NOT TARGET POINTS
  //Calculate values with new function
  ierr = CalcCoulombEllCoefs(&e, npts, ellPoints, srcCharges, NMAX, &GnpVals);CHKERRQ(ierr); //testing
  ierr = CalcReactAndExtCoefsFromCoulomb(&e, prob.e1, prob.e2, NMAX, GnpVals, BnpVals, CnpVals);CHKERRQ(ierr);

  //Calculate values with old function
  PetscReal GnpOld, BnpOld, CnpOld;
  count = 0;
  for(PetscInt n=0; n<=NMAX; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p) {
      GnpOld = calcGnp(&e, chargePoints, chargeValues, npts, n, p);
      calcBnpAndCnpFromGnp(&prob, n, p, GnpOld, &BnpOld, &CnpOld);
      ierr = VecSetValues(GnpValsOld, 1, &count, &GnpOld, INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(BnpValsOld, 1, &count, &BnpOld, INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(CnpValsOld, 1, &count, &CnpOld, INSERT_VALUES);CHKERRQ(ierr);

      count++;
    }
  }
  ierr = VecAssemblyBegin(GnpValsOld);CHKERRQ(ierr); ierr = VecAssemblyEnd(GnpValsOld);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(BnpValsOld);CHKERRQ(ierr); ierr = VecAssemblyEnd(BnpValsOld);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(CnpValsOld);CHKERRQ(ierr); ierr = VecAssemblyEnd(CnpValsOld);CHKERRQ(ierr);

  PetscReal GnpErr, BnpErr, CnpErr;
  ierr = VecWAXPY(errorV, -1.0, GnpVals, GnpValsOld);CHKERRQ(ierr);
  ierr = VecNorm(errorV, NORM_2, &GnpErr);
  ierr = VecWAXPY(errorV, -1.0, BnpVals, BnpValsOld);CHKERRQ(ierr);
  ierr = VecNorm(errorV, NORM_2, &BnpErr);
  ierr = VecWAXPY(errorV, -1.0, CnpVals, CnpValsOld);CHKERRQ(ierr);
  ierr = VecNorm(errorV, NORM_2, &CnpErr);
  printf("GnpErr: %3.3e\n", GnpErr);
  printf("BnpErr: %3.3e\n", BnpErr);
  printf("CnpErr: %3.3e\n", CnpErr);
  
  
  //solution vector
  solution = (double*) malloc(sizeof(double)*npts);
  
  //calculate potential on solPoints
  calcCoulombEllipsoidalGrid(&prob, 10, solPoints, npts, solution);
  printf("wow!\n");
  for(int i=0; i<npts; ++i)
    printf("solution[%d] = %15.15f\n", i, solution[i]);
  for(int i=0; i<npts; ++i)
    printf("chargeValues[%d] = %15.15f\n", i, chargeValues[i]);
  
  //calculate free energy
  freeEnergy = 0;
  for(int i=0; i<npts; ++i) {
    freeEnergy += solution[i]*chargeValues[i];
    printf("freeEnergy[%d] = %15.15f\n", i, freeEnergy*.5*cf);
  }

  
  free(xyz);
  free(chargeValues);
  free(solPoints);
  free(solution);
  free(chargePoints);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__
int main( int argc, char **argv ) {

  PetscErrorCode ierr;

  //initialize Petsc
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);  

  
  //testSphericalCoordinates();
  //runTest1();
  //testLegendre();
  //runTest4();
  //test4sphere();
  //RunArg();
  ierr = RunArgTester(); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
