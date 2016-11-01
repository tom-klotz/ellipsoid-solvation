#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "ellipsoid.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

typedef struct Context {
  int l;
  void (*funcPtr)(mpfr_t*,mpfr_t*);
} Context;

/*
double chargeDensity(double theta) {
  return cos(2*theta);
}
*/
void chargeDensity(mpfr_t *x, mpfr_t *val)
{
  mpfr_mul_d(*val, *x, 2.0, MPFR_RNDN);
  mpfr_cos(*val, *val, MPFR_RNDN);
}

void rFunc(mpfr_t *x, mpfr_t *val, int *l)
{
  //calculates x^l
  mpfr_log(*val, *x, MPFR_RNDN);
  mpfr_mul_d(*val, *val, (double) *l, MPFR_RNDN);
  mpfr_exp(*val, *val, MPFR_RNDN);
}

void tFunc(mpfr_t *x, mpfr_t *val, int *l)
{
  //calculates P_l(cos(x)) * rho(x)
  mpfr_t temp;
  mpfr_init2(temp, mpfr_get_prec(*x));

  mpfr_cos(*val, *x, MPFR_RNDN);
  mpfr_set_d(*val, boost::math::legendre_p(*l, mpfr_get_d(*val, MPFR_RNDN)), MPFR_RNDN);
  
  //multiply by charge density
  chargeDensity(x, &temp);
  mpfr_mul(*val, *val, temp, MPFR_RNDN);
}

double computeExpansionCoefficient(int l)
{
  mpfr_t a, b;
  mpfr_inits(a, b, NULL);
  
  //a=0, b=2*pi
  mpfr_set_d(a, 0.0, MPFR_RNDN);
  mpfr_const_pi(b, MPFR_RNDN);
  mpfr_mul_d(b, b, 2.0, MPFR_RNDN);
  
  double tintegral;
  //perform theta integration
  integrateMPFR((void (*)(mpfr_t *,mpfr_t*,void*)) tFunc, a, b, 16, &tintegral, (void*) &l);
  
  //a=0, b=1
  mpfr_set_d(b, 1.0, MPFR_RNDN);
  
  double rintegral;
  //perform r integration
  integrateMPFR((void (*)(mpfr_t*,mpfr_t*,void*)) rFunc, a, b, 16, &rintegral, (void*) &l);

  return tintegral*rintegral;
}

int main(int argc, char **argv)
{
  //terms in expansion
  int terms = 1000;
  
  //radius of sphere
  double R = 1.0;

  //size of grid
  int rdim = 25;
  int tdim = 100;
  
  //grid boundary values
  double rmin = R;
  double rmax = 5*R;
  double tmin = 0;
  double tmax = 2*M_PI;

  //initialize r,theta values
  double *rvals = (double*) malloc(sizeof(double)*rdim);
  double *tvals = (double*) malloc(sizeof(double)*tdim);
  for(int i=0; i<rdim; ++i)
    rvals[i] = rmin + (rmax-rmin)*i/(rdim-1);
  for(int i=0; i<tdim; ++i)
    tvals[i] = tmin + (tmax-tmin)*i/(tdim-1);

  //set up solution matrix
  double *sol = (double*) calloc(sizeof(double),rdim*tdim);
  //error vector
  double *oldsol = (double*) calloc(sizeof(double), rdim*tdim);


  //compute expansion coefficients
  double *coeffs = (double*) malloc(sizeof(double)*terms);
  for(int l=0; l<terms; ++l) {
    coeffs[l] = computeExpansionCoefficient(l);
    printf("expansion coefficient %d calculated as %8.8f\n", l, coeffs[l]);
  }

  double error = 0;
  //series loop
  for(int l=0; l<terms; ++l) {
    if(l%2==0)
      error = 0;
    memcpy(oldsol, sol, sizeof(double)*rdim*tdim);
    //update for evaluation point
    for(int i=0; i<rdim; ++i) {
      for(int j=0; j<tdim; ++j) {
	sol[i*rdim+j] += (coeffs[l]/rvals[i]) * boost::math::legendre_p(l,cos(tvals[j]));
	if(l%2==0)
	  error += fabs(oldsol[i*rdim+j] - sol[i*rdim+j]) * fabs(oldsol[i*rdim+j] - sol[i*rdim+j]);
      }
    }
    error = sqrt(error);
    printf("error at level %d is %8.8e\n", l, error);
  }
  
  printf("hoorah!\n");
  free(rvals);
  free(tvals);
  free(coeffs);
  free(sol);
  return 0;
}
