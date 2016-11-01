#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include "ellipsoid.h"
#include "constants.h"
#include <boost/math/special_functions/legendre.hpp>
//#include "legendre_polynomial.h"

typedef struct SPoint {

  double x1;
  double x2;
  double x3;
  char type;

} SPoint;

typedef struct Problem {
  Point *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
  EllipsoidalSystem *e;
} Problem;

typedef struct SProblem {
  SPoint *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
  double b;  //radius of molecule
} SProblem;

//output leg is vector of length nz*(2*n+1)
void legendre2(int n, int nz, double z, double *leg)
{
  double fac = 1.0;
  double sqz2   = sqrt(1.0 - z*z);
  double hsqz2  = 0.5*sqz2;
  double ihsqz2 = z/hsqz2;
  for(int i=2; i<=n; ++i)
    fac *= i;
  double hsqz2ton = 1.0;
  for(int i=0; i<n; ++i)
    hsqz2ton *= hsqz2;
  int pre = n % 2 ? -1 : 1;

  
  if(sqz2 == 0) {
    //printf("WOW!\n");
    for(int mr = 0; mr < 2*n+1; ++mr)
      leg[mr] = boost::math::legendre_p(n, mr-n, z);
  }
  else if(n==0) {
    leg[0] = 1.0;
  }
  else if(n==1) {
    leg[0] = -hsqz2;
    leg[1] = z;
    leg[2] = sqz2;
  }
  else {
    leg[0] = (1.0 - 2.0*fabs(n - 2.0*floor(n/2.))) * hsqz2ton/fac;
    leg[1] = -leg[0]*n*ihsqz2;
    for(int mr = 1; mr < 2*n; ++mr)
      leg[mr+1] = (mr - n)*ihsqz2*leg[mr] - (2.0*n - mr + 1)*mr*leg[mr-1];
  }
  for(int m=0; m<2*n; ++m, pre = -pre) leg[m] *= pre;
}


void cartesianToSpherical(SPoint *point) {
  
  double r, theta, phi;

  r = sqrt((point->x1)*(point->x1) + (point->x2)*(point->x2) + (point->x3)*(point->x3));
  
  //if(point->x3/r == 1)
  //theta = 0;
  //else if(point->x3/r == -1)
  //theta = M_PI;
  //else
  theta = acos(point->x3/r);


  phi = atan2((point->x2),(point->x1));
  //if(phi < 0)
  //phi = phi+2*M_PI;

  point->x1 = r;
  point->x2 = 0;
  point->x3 = phi;

  if(r > 0) point->x2 = theta;

  point->type = 's';

}

void sphericalToCartesian(SPoint *point) {

  double r = point->x1;
  double theta = point->x2;
  double phi = point->x3;

  double x = r*sin(theta)*cos(phi);
  double y = r*sin(theta)*sin(phi);
  double z = r*cos(theta);

  point->x1 = x;
  point->x2 = y;
  point->x3 = z;

}

double calcEnp(EllipsoidalSystem *e, Point *point, int n, int p) {
  //if(point->type == 'c')
  //cartesianToEllipsoidal(e, point);
  int signm = 1;
  int signn = 1;
  if(point->x2 < 0)
    signm = -1;
  if(point->x3 < 0)
    signn = -1;
  double eL = calcLame(e, n, p, point->x1, signm, signn);
  double eM = calcLame(e, n, p, point->x2, signm, signn);
  double eN = calcLame(e, n, p, point->x3, signm, signn);

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
  double ival   = calcI(e, n, p, point->x1, signm, signn);
  //printf("ENP: %15.15f\nI: %15.15f\n", enpval, ival);
  return (2*n + 1) * enpval * ival;
}

void calcBnpAndCnpFromGnp(Problem *problem, int n, int p, double Gnp, double *Bnp, double *Cnp)
{
  EllipsoidalSystem *e = problem->e;
  double Ea = calcLame(e, n, p, e->a, 1, 1);
  double Ia = calcI(problem->e, n, p, e->a, 1, 1);
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


void calcEnm(SPoint *positions, double *charges, int nCharges, int n, int m, double *rval, double *ival) {
  
  int print = 0;
  if(n < 2) {
    printf("\n\nn = %d\nm = %d\n", n, m);
    print = 0;
  }

  int absm = abs(m);
  double coef = 1;
  //compute factorial coefficient
  for(int i = n - absm + 1; i <= n + absm; ++i)
    coef *= i;
  coef = 1./coef;

  *rval = 0;
  *ival = 0;



  double qk, rk, rkton, thetak, phik, lterm, rterm, iterm;
  //double *legvals = (double*) malloc(sizeof(double)*(2*n+1));
  //sum over charges
  for(int k=0; k<nCharges; ++k) {
    
    qk = charges[k];
    rk = positions[k].x1;
    thetak = positions[k].x2;
    phik = positions[k].x3;
    

    rkton = 1;
    for(int i=0; i<n; ++i)
      rkton = rkton*rk;
    
    //legendre2(n, 1, cos(thetak), legvals);

    lterm = boost::math::legendre_p(n, absm, cos(thetak));//legvals[absm+n];//boost::math::legendre_p(n, m, cos(thetak));
    //if(m%2 == 1)
    //lterm *= -1;
    //eterm = cos(m*phi)*cos(m*phik) + sin(m*phi)*sin(m*phik);
    rterm = cos(m*phik);
    iterm = -sin(m*phik);
    //if(sin(m*phi)*sin(m*phik) != 0)
    //printf("WOW: %15.15f\n", sin(m*phi)*sin(m*phik));

    if(print) {
      printf("QK: %15.15f\n", qk);
      printf("rkton: %15.15f\n", rkton);
      printf("ltern: %15.15f\n", lterm);
      //printf("eterm; %15.15f\n", eterm);
      
    }

    //update val
    *rval += qk * coef * rkton * lterm * rterm;
    *ival += qk * coef * rkton * lterm * iterm;
    //if(print)
    //printf("VAL after %dth point: %15.15f\n", k, val);
  }

  
}

void calcBnmFromEnm(int n, int m, double Enmrval, double Enmival, SProblem *problem, double *Bnmrval, double *Bnmival) {
  double e1 = problem->e1;
  double e2 = problem->e2;
  double b  = problem->b;
  
  double val = 1.0;
  for(int i=0; i<2*n+1; ++i) {
    val /= b;
  }

  val /= e1;
  
  val *= (e1-e2)*(n+1);
  val /= (e1*n + e2*(n+1));
  //printf("VAL: %15.15f\n", val);
  *Bnmrval = val*Enmrval;
  *Bnmival = val*Enmival;

}

void calcCnmFromEnm(int n, int m, double Enmrval, double Enmival, SProblem *problem, double *Cnmrval, double *Cnmival) {
  double e1 = problem->e1;
  double e2 = problem->e2;
  
  double val;
  
  val = n*(e1-e2)/(n*e1+e2*(n+1));
  val = (1./e2)*(1-val);
  
  *Cnmrval = val * Enmrval;
  *Cnmival = val * Enmival;

}


double calcCoulombSpherical(SProblem *problem, int maxDeg, SPoint *r) {
  double val = 0;
  double Enm;
  double e1 = problem->e1;
  double rad = r->x1;
  double theta = r->x2;
  double phi = r->x3;
  double radpow = 1.0;
  double lterm;
  double Enmrterm, Enmiterm;
  double rterm, iterm;
  double *legvals = (double*) malloc(sizeof(double)*(maxDeg*2 + 1));

  for(int n=0; n<=maxDeg; ++n) {    
    //printf("cos(theta) = %15.15f\n", cos(theta));
    //legendre2(n, 1, cos(theta), legvals);
    radpow *= rad;
    for(int m=-n; m<=n; ++m) {
      calcEnm(problem->positions, problem->charges, problem->nCharges, n, m, &Enmrterm, &Enmiterm);
      //printf("Enm %d %d %15.15f\n", n, m, Enm);

      //printf("cos(theta): %15.15f\n", theta);
      lterm = boost::math::legendre_p(n, abs(m), cos(theta));//legvals[n+abs(m)];//boost::math::legendre_p(n, m, cos(theta));
      
      rterm = cos(m*phi);
      iterm = sin(m*phi);
      
      val += ((Enmrterm*rterm - Enmiterm*iterm)*lterm)/(e1*radpow);
      //printf("lterm: %15.15f\n", lterm);
      //printf("cos(mphi) = %15.15f\n\n", cos(m*phi));
    }
    printf("n = %d\nval = %15.15f\n", n, val);
  }

  return 0.0;
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
      //clock_t start = clock(), diff;
      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);
      calcBnpAndCnpFromGnp(problem, n, p, Gnp, &Bnp, &Cnp);
      //diff = clock() - start;
      //int msec = diff * 1000 / CLOCKS_PER_SEC;
      //printf("Gnp took %d seconds %d milliseconds\n", msec/1000, msec%1000);
      //start = clock();
     
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
      //diff = clock() - start;
      //msec = diff * 1000 / CLOCKS_PER_SEC;
      //printf("Fnp took %d seconds %d milliseconds\n\n", msec/1000, msec%1000);
      //fprintf(fp, "%8.8e %8.8e\n", (double)count, fabs(Fnp));



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


double calcCoulombSphericalGrid(SProblem *problem, int maxDeg, SPoint *r, int nPoints, double *solution) {
  double val = 0;
  double Enm;
  double e1 = problem->e1;
  
  double rad, theta, phi;
  double lterm;
  double Enmrval, Enmival;
  double Bnmrval, Bnmival;
  double Cnmrval, Cnmival;
  double rval, ival;

  for(int i=0; i<nPoints; ++i)
    solution[i] = 0.0;
  for(int n=0; n<=maxDeg; ++n) {    
    for(int m=-n; m<=n; ++m) {

      calcEnm(problem->positions, problem->charges, problem->nCharges, n, m, &Enmrval, &Enmival);
      calcBnmFromEnm(n, m, Enmrval, Enmival, problem, &Bnmrval, &Bnmival);
      calcCnmFromEnm(n, m, Enmrval, Enmival, problem, &Cnmrval, &Cnmival);

      for(int k=0; k<nPoints; ++k) {
	
	rad   = r[k].x1;
	theta = r[k].x2;
	phi   = r[k].x3;

	lterm = boost::math::legendre_p(n, abs(m), cos(theta));

	rval = cos(m*phi);
	ival = sin(m*phi);

	if(rad < problem->b) {

	  
	  solution[k] += (Bnmrval*rval - Bnmival*ival)*lterm*pow(rad,n);

	}
	else if(rad >  problem->b) {
	  


	  solution[k] += ((Cnmrval*rval - Cnmival*ival)*lterm)/pow(rad,n+1);

	}
	else {
	  
	  printf("CALCULATING ON BOUNDARY\n");
	  solution[k] += ((Enmrval*rval - Enmival*ival)*lterm)/(e1*pow(rad,n+1));
	}
      }
    }
    printf("n = %d\nval = %15.15f\n", n, solution[0]);
  }

  return 0.0;
}

double calcRegion2(SProblem *problem, int maxDeg, SPoint *r) {
  double val = 0;
  double Enm, Cnm;
  double e1 = problem->e1;
  double e2 = problem->e2;
  double rad = r->x1;
  double theta = r->x2;
  double phi = r->x3;
  double radpow = 1;
  double lterm;
  double rterm, iterm;
  for(int n=0; n<=maxDeg; ++n) {
    radpow *= rad;
    for(int m=-n; m<=n; ++m) {
      calcEnm(problem->positions, problem->charges, problem->nCharges, n, m, &rterm, &iterm);
      //Cnm = calcCnmFromEnm(m, n, rterm, problem);
      Cnm = 0;
      if(fabs(Enm-Cnm) > 1e-14)
	printf("\n\n\nWHATWHATWHAT\n\n\n");
      lterm = boost::math::legendre_p(n, m, cos(theta));
      //if(m%2 == 1)
      //lterm *= -1;
      val += Cnm*lterm/radpow;
    }
    printf("n: %d\nval: %15.15f\n\n", n, val);

  }
  return val;
}

double calcReactionPotential(SProblem *problem, int maxDeg, SPoint *r) {
  double val = 0;
  double Enm, Bnm;
  double e1 = problem->e1;
  double e2 = problem->e2;
  double rad = r->x1;
  double theta = r->x2;
  double phi = r->x3;
  
  
  return 0;
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

void testSphericalCoordinates() {

  const int numRandPoints = 20;

  SPoint testPoints[numRandPoints];
  SPoint testPointsCpy[numRandPoints];
  
  srand(time(NULL));
  for(int i=0; i < numRandPoints; ++i) {
    testPoints[i].x1 = (double)rand()/(double)RAND_MAX;
    testPoints[i].x2 = (double)rand()/(double)RAND_MAX;
    testPoints[i].x3 = (double)rand()/(double)RAND_MAX;
    testPointsCpy[i].x1 = testPoints[i].x1;
    testPointsCpy[i].x2 = testPoints[i].x2;
    testPointsCpy[i].x3 = testPoints[i].x3;
  }
  
  //convert to spherical and back to cartesian and print error
  for(int i=0; i < numRandPoints; ++i) {
    cartesianToSpherical(testPoints+i);
    sphericalToCartesian(testPoints+i);

    printf("error[%d].x1 = %8.8e\n", i, fabs(testPoints[i].x1-testPointsCpy[i].x1));
    printf("error[%d].x2 = %8.8e\n", i, fabs(testPoints[i].x2-testPointsCpy[i].x2));
    printf("error[%d].x3 = %8.8e\n", i, fabs(testPoints[i].x3-testPointsCpy[i].x3));
  }
  

}


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
  SProblem sProblem = {.positions = sPoints, .charges = charges, .nCharges = nCharges, .e1 = 1.0, .e2 = 1.0, .b = 3.0};
  Problem  eProblem = {.positions = ePoints, .charges = charges, .nCharges = nCharges, .e1 = 1.0, .e2 = 1.0, .e = &e};
  
 
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
  //double zh = (zb - za)/(nz-1);
  double zh = 0;


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
  //double zh = (zb - za)/(nz-1);
  double zh = 0;


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
  //double zh = (zb - za)/(nz-1);
  double zh = 0;


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

  double charge = 1.0;
  //set up problem
  Problem problem = {.positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0, .e = &e};


  int iter = 99;
  
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

  


  double charge = 1.0;
  //set up problem
  SProblem problem = {.positions = points, .charges = charges, .nCharges = npoints, .e1 = 1.0, .e2 = 1.0};


  calcCoulombSpherical(&problem, 300, &r);
  
  //printf("\n\n\nPART 2 PART 2 PART 2 PART 2\n\n\n");
  
  //calcCoulombSpherical(&problem, 30, &r);

  //speedTester(&problem, 10, &r);
  //printf("norm for %d,%d: %8.8e\n", 25, 39, calcNormalization(&e, 25, 39));
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
  double ellSol = calcCoulomb(&problem, 32, &r);
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
  double c = 1.0;


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


void testLegendre()
{
  int Nmax = 33;
  double xval = .564;
  double lterm;
  for(int n=0; n<Nmax; ++n) {
    for(int m=-n; m<=n; ++m) {
      lterm = boost::math::legendre_p(n, m, xval);
      printf("%d,%d: %15.15f\n", n, m, lterm);
    }
  }
}


int runArg()
{
  char buff[64];
  int npts;
  double a, b, c;
  double *xyz;
  double *chargeValues;
  double *solution;
  double freeEnergy;

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
  return 0;
}


int main(){
  /*
  int n = 10;
  double val = -.6;
  double *legvals = (double*) malloc(sizeof(double)*(2*n+1));
  legendre2(n, 1, val, legvals);
  for(int i=0; i<2*n+1; ++i) {
    printf("legendre2(%d, %d, %3.3e) = %15.15f\n", n, i-n, val, legvals[i]);
    printf("boost:                     %15.15f\n", boost::math::legendre_p(n, i-n, val));
  }
  */
  //testSphericalCoordinates();
  //runTest1();
  //testLegendre();
  //runTest4();
  //test4sphere();
  runArg();
  return 0;
}
