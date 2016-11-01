#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include "sphere.h"
#include "sphSolv.h"

void calcEnm(SPoint *positions, double *charges, int nCharges, int n, int m, double *rval, double *ival)
{
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
    

    //seems like we don't have to worry about (-1)^m? if(m%2 ==1) lterm*=-1; not necessary
    lterm = gsl_sf_legendre_Plm(n, absm, cos(thetak));//boost::math::legendre_p(n, absm, cos(thetak));
    
    rterm = cos(m*phik);
    iterm = -sin(m*phik);
    
    if(print) {
      printf("QK: %15.15f\n", qk);
      printf("rkton: %15.15f\n", rkton);
      printf("ltern: %15.15f\n", lterm);      
    }
    
    //update val
    *rval += qk * coef * rkton * lterm * rterm;
    *ival += qk * coef * rkton * lterm * iterm;
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
  double e1 = problem->e1;
  double rad = r->x1;
  double theta = r->x2;
  double phi = r->x3;
  double radpow = 1.0;
  double lterm;
  double Enmrterm, Enmiterm;
  double rterm, iterm;

  for(int n=0; n<=maxDeg; ++n) {    
    radpow *= rad;
    for(int m=-n; m<=n; ++m) {
      calcEnm(problem->positions, problem->charges, problem->nCharges, n, m, &Enmrterm, &Enmiterm);
      //printf("Enm %d %d %15.15f\n", n, m, Enm);

      //printf("cos(theta): %15.15f\n", theta);
      lterm = gsl_sf_legendre_Plm(n, abs(m), cos(theta));//boost::math::legendre_p(n, abs(m), cos(theta));//legvals[n+abs(m)];//boost::math::legendre_p(n, m, cos(theta));
      
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


double calcCoulombSphericalGrid(SProblem *problem, int maxDeg, SPoint *r, int nPoints, double *solution) {
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

	lterm = gsl_sf_legendre_Plm(n, abs(m), cos(theta));//boost::math::legendre_p(n, abs(m), cos(theta));

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
