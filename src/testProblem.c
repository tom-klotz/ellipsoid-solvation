#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "ellipsoid.h"

typedef struct Problem {
  EllipsoidalSystem *e;
  Point *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
} Problem;

/*
//Uses direct computation of polynomial coefficients rather than recursive
double legendre(int l, int m, double x) {
  
  int maxcf;       //largest polynomial coefficient
  double cfnm = 1; //proportionality constant for m<0 polynomials compared to m>0
  double cl, x2, p, f1, f2, f3;
  int px, j;
  
  if( l<0 ) l = -(l+1);
  if( m<0 ) {
    double num, den;

    m = -m;
    num = factorial(l-m);
    den = factorial(l+m);
    cfnm = powInt(-1,m)*num/den;
  }
  //calculate coef of maximum degree in x from the explicit analytical formula
  f1 = factorial(2*l);
  f2 = factorial(l);
  f3 = factorial(l-m);
  cl = powInt(-1,m) * cfnm * f1/((1 << l)*f2*f3);
  
  return 0;
}
*/

//
void legendre2(int l, int nz, double z, double *leg)
{
  double sqz2   = sqrt(1.0 - z*z);
  double hsqz2  = 0.5*sqz2;
  double ihsqz2 = fabs(z)/hsqz2;
  double fac    = 1.0;
  int pre       = l % 2 ? -1 : 1;
  int m = 0;

  for(m = 2; m <= l; ++m) fac *= m;
  if(!l) {
    leg[0] = 1.0;
  } else if(l == 1) {
    leg[0] = -hsqz2;
    leg[1] = z;
    leg[2] = sqz2;
  } else {
    leg[0] = (1.0 - 2.0*fabs(l - 2.0*floor(l/2.0)))*pow(hsqz2, l)/fac;
    leg[1] = -leg[0]*l*ihsqz2;
    for(m = 1; m < 2*l; ++m) leg[m+1] = (m - l)*ihsqz2*leg[m] - (2*l - m + 1)*m*leg[m-1];
  }
  for(m = 0; m <= 2*l; ++m, pre = -pre) leg[m] *= pre;
  
}

double calcEnp(EllipsoidalSystem *e, Point *point, int n, int p) {
  if(point->type == 'c')
    cartesianToEllipsoidal(e, point);
  
  double eL = calcLame(e, n, p, point->x1);
  double eM = calcLame(e, n, p, point->x2);
  double eN = calcLame(e, n, p, point->x3);

  return eL*eM*eN;
}


double calcGnp(EllipsoidalSystem *e, Point *positions, double *charges, int nCharges, int n, int p)
{
  double normConstant = calcNormalization(e, n, p);
  //double normConstant = 1.0;
  double sum = 0;
  for(int k=0; k<nCharges; ++k) {
    sum += charges[k] * calcEnp(e, positions+k, n, p);
  }

  sum *= (4*M_PI)/((2*n+1)*normConstant);

  return sum;
}

double calcFnp(EllipsoidalSystem *e, Point *point, int n, int p)
{ 
  return (2*n + 1) * calcEnp(e, point, n, p) * calcI(e, n, p, point->x1);
}

double calcBnpFromGnp(EllipsoidalSystem *e, int n, int p, double Gnp)
{
  double Fa = 1;
  double Ea = calcLame(e, n, p, e->a);
  double EaDer = calcLameDerivative(e, n, p, e->a);
  double FaDer = (2*n+1) * (calcLame(e, n, p, e->a) * calcIDerivative(e, n, p, e->a) +
			    EaDer                   * calcI(e, n, p, e->a)); //chain rule?
  return 0;
}

double calcCoulomb(Problem *problem, int maxDeg, Point *r)
{
  //FILE *fp = fopen("fnpnorm.txt", "w");
  double Gnp, Fnp, Enp;
  double sum = 0;
  int count = 0;
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p < 2*n+1; ++p) {

      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);

      Fnp = calcFnp(problem->e, r, n, p);
      
      //fprintf(fp, "%8.8e %8.8e\n", (double)count, fabs(Fnp));

      sum += (Gnp*Fnp)/problem->e1;

      count++;
    }
      printf("n: %d\n", n);
      printf("Gnp: %16.16f\n", Gnp);
      printf("Fnp: %16.16f\n", Fnp);
      printf("enp: %16.16f\n", Enp);
      printf("sum: %16.16f\n", sum);

  }
  //fclose(fp);
  return sum;
}

double calcCoulombSpherical(Problem *problem, int maxDeg, Point *r) {
  //  for(int i=0; i<10; ++i)
  //printf("P%d^1(.7) = %15.15f\n", i, legendre((unsigned)i, .7));
  return 0.0;
}

void speedTester(Problem *problem, int maxDeg, Point *r)
{
  double Gnp, Fnp, Enp;
  printf("gnp starting\n");
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);
      //Gnp = calcI(problem->e, n, p, .67);
    }
  }
  printf("gnp done\n");

  printf("fnp starting\n");  
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      //calcI(problem->e, n, p, .67);
      //Fnp = (2*n + 1) * calcEnp(problem->e, problem->positions, n, p) * calcI(problem->e, n, p, problem->positions->x1);
      Fnp = calcFnp(problem->e, r, n, p);
    }
  }
  printf("fnp done\n");
}

void test1() {
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  

  const int npoints = 2;

  //calculation point
  Point r = {.x1 = 0.0, .x2 = 0.0, .x3 = 10.0, .type = 'c'};


  //source charges
  Point points[npoints];
  points[0].x1 = 0.3;
  points[0].x2 = 0.2;
  points[0].x3 = 0.5;
  points[1].x1 = 0.01;
  points[1].x2 = 0.1;
  points[1].x3 = 0.1;
  
  double charges[npoints];
  charges[0] = 1.0;
  charges[1] = 1.0;


  double exact = 0;
  for(int i=0; i<npoints; i++)
    exact += 1./sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

  printf("the exact result is %15.15f\n", exact);


  cartesianToEllipsoidal(&e, &r);
  cartesianToEllipsoidal(&e, points);
  cartesianToEllipsoidal(&e, points+1);

  for(int i=0; i<npoints; ++i) {
    printf("points[%d].x1 = %15.15f\n", i, points[i].x1);
    printf("points[%d].x2 = %15.15f\n", i, points[i].x2);
    printf("points[%d].x3 = %15.15f\n", i, points[i].x3);
  }

  double charge = 1.0;
  //set up problem
  Problem problem = {.e = &e, .positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0};


  int iter = 99;
  
  calcCoulomb(&problem, 12, &r);


  //speedTester(&problem, 10, &r);
  //printf("norm for %d,%d: %8.8e\n", 25, 39, calcNormalization(&e, 25, 39));
}

void test2() {
  EllipsoidalSystem e;
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  initEllipsoidalSystem(&e, a, b, c);
  

  const int npoints = 5;

  //calculation point
  Point r = {.x1 = -0.1, .x2 = 0.1, .x3 = 3.5, .type = 'c'};

  //source charges
  Point points[npoints*npoints];
  //generate points on the surface of a 2.9, 1.9, .9 ellipse
  double distIn = .2;
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
      points[i*5+j].x1 = pointsA*cos(u)*cos(v);
      points[i*5+j].x2 = pointsB*cos(u)*sin(v);
      points[i*5+j].x3 = pointsC*sin(u);
      points[i*5+j].type = 'c';
    }
  }

  //charges array, all set to 1
  double charges[npoints*npoints];
  for(int i=0; i<npoints*npoints; ++i)
    charges[i] = 1.0;

  //calculate exact solution ahead of time
  double exact = 0;
  for(int i=0; i<npoints*npoints; i++)
    exact += 1./sqrt((points[i].x1-r.x1)*(points[i].x1-r.x1) + (points[i].x2-r.x2)*(points[i].x2-r.x2) + (points[i].x3-r.x3)*(points[i].x3-r.x3));

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
  Problem problem = {.e = &e, .positions = points, .charges = charges, .nCharges = npoints*npoints, .e1 = 1.0, .e2 = 1.0};

  //calculate the coulomb potential
  double ellSol = calcCoulomb(&problem, 10, &r);
}


int main(){
  test2();
  
  return 0;
}
