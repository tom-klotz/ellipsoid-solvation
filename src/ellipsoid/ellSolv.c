#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <petsc.h>
#include "ellipsoid.h"
#include "../sphere/sphere.h"
#include "ellSolv.h"


double calcEnp(EllipsoidalSystem *e, Point *point, int n, int p) {
  //if(point->type == 'c')
  //cartesianToEllipsoidal(e, point);
  int signm = 1;
  int signn = 1;
  if(point->x2 < 0)
    signm = -1;
  if(point->x3 < 0)
    signn = -1;
  double eL, eM, eN;
  calcLame(e, n, p, point->x1, signm, signn, &eL);
  calcLame(e, n, p, point->x2, signm, signn, &eM);
  calcLame(e, n, p, point->x3, signm, signn, &eN);

  return eL*eM*eN;
}

double calcGnp(EllipsoidalSystem *e, Point *positions, double *charges, int nCharges, int n, int p)
{
  clock_t start, diff;
  int msec;
  start = clock();
  double normConstant = calcNormalization(e, n, p);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("normConstant took %d seconds and %d milliseconds\n", msec/1000, msec%1000);
  //double normConstant = 1.0;
  double sum = 0;
  double Enp;

  for(int k=0; k<nCharges; ++k) {

    Enp = calcEnp(e, positions+k, n, p);



    sum += charges[k] * Enp;
  }

  sum *= (4*M_PI)/((2.0*n+1.0)*normConstant);

  return sum;
}


double calcFnp(EllipsoidalSystem *e, Point *point, int n, int p)
{ 
  int signm = 1;
  int signn = 1;
  if(point->x2 < 0)
    signm = -1;
  if(point->x3 < 0)
    signn = -1;

  double enpval = calcEnp(e, point, n, p);
  double ival;
  calcI(e, n, p, point->x1, signm, signn, &ival);
  //printf("ENP: %15.15f\nI: %15.15f\n", enpval, ival);
  return (2*n + 1) * enpval * ival;
}


void calcBnpAndCnpFromGnp(Problem *problem, int n, int p, double Gnp, double *Bnp, double *Cnp)
{
  EllipsoidalSystem *e = problem->e;
  double Ea;
  calcLame(e, n, p, e->a, 1, 1, &Ea);
  double Ia;
  calcI(problem->e, n, p, e->a, 1, 1, &Ia);
  double Fa = (2*n + 1) * Ea * Ia;
  double EaDer = calcLameDerivative(e, n, p, e->a, 1, 1);
  double IaDer = calcIDerivative(e, n, p, e->a, 1, 1);
  double FaDer = (2*n+1) * (Ea*IaDer + EaDer*Ia); //chain rule?
  double e1 = problem->e1;
  double e2 = problem->e2;
  
  //calc Bnp
  double first = (Fa/Ea)*(e1 - e2)/(e1*e2);
  first /= (1-(e1/e2)*((EaDer*Fa)/(FaDer*Ea)));
  *Bnp = first*Gnp;

  //calc Cnp
  first = (e1/e2)*(EaDer/FaDer)*(*Bnp);
  first += (Gnp/e2);
  *Cnp = first;
}



double calcCoulomb(Problem *problem, int maxDeg, Point *r)
{
  //FILE *fp = fopen("fnpnorm.txt", "w");
  double Gnp, Fnp, Enp;
  double sum = 0;
  int count = 0;
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p < 2*n+1; ++p) {
      clock_t start = clock(), diff;
      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);
      diff = clock() - start;
      int msec = diff * 1000 / CLOCKS_PER_SEC;
      //printf("Gnp took %d seconds %d milliseconds\n", msec/1000, msec%1000);
      start = clock();
      Fnp = calcFnp(problem->e, r, n, p);
      diff = clock() - start;
      msec = diff * 1000 / CLOCKS_PER_SEC;
      //printf("Fnp took %d seconds %d milliseconds\n\n", msec/1000, msec%1000);
      //fprintf(fp, "%8.8e %8.8e\n", (double)count, fabs(Fnp));

      sum += (Gnp*Fnp)/problem->e1;

      count++;
    }
    printf("n: %d\n", n);
    printf("Gnp: %16.16f\n", Gnp);
    printf("Fnp: %16.16f\n", Fnp);
    //printf("enp: %16.16f\n", Enp);
    printf("sum: %16.16f\n", sum);
    
  }
  //fclose(fp);
  return sum;
}







double calcCoulombEllipsoidalGrid(Problem *problem, int maxDeg, Point *r, int nPoints, double *solution) {

  EllipsoidalSystem *e = problem->e;
  double Gnp, Bnp, Fnp, Enp, Cnp;
  for(int i=0; i<nPoints; ++i)
    solution[i] = 0;
  int count = 0;
  double *fnpvals = (double*) malloc(sizeof(double)*nPoints);
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p < 2*n+1; ++p) {
      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);
      calcBnpAndCnpFromGnp(problem, n, p, Gnp, &Bnp, &Cnp);
     
      for(int k=0; k<nPoints; ++k) {
	//if on exterior
	if(fabs(r[k].x1) >= e->a) {
	  Fnp = calcFnp(problem->e, r+k, n, p);
	  solution[k] += Cnp*Fnp;
	}
	else {
	  Enp = calcEnp(problem->e, r+k, n, p);
	  solution[k] += Bnp*Enp;
	}
      }



      count++;
    }
    printf("n: %d\n", n);
    printf("Gnp: %16.16f\n", Gnp);
    printf("Fnp: %16.16f\n", Fnp);
    printf("enp: %16.16f\n", Enp);
    printf("sum: %16.16f\n", solution[0]);
    
  }
  //fclose(fp);
  return solution[0];
}




