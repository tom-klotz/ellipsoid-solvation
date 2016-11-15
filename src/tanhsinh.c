#include <math.h>
#include <stdio.h>
#include <mpfr.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <petsc.h>


#undef __FUNCT__
#define __FUNCT__ "DEwk"
PetscErrorCode DEwk(PetscInt k, mpfr_t h, mpfr_t *wk)
{
  PetscErrorCode ierr;
  mpfr_t kh, msinh, mcosh, piOver2, tmp;
  PetscFunctionBegin;
  
  mpfr_inits(kh, msinh, mcosh, piOver2, tmp, NULL);

  //init piOver2
  mpfr_const_pi(piOver2, MPFR_RNDN);
  mpfr_mul_d(piOver2, piOver2, 0.5, MPFR_RNDN);
  
  mpfr_set_d(kh, (double) k, MPFR_RNDN);
  mpfr_mul(kh, kh, h, MPFR_RNDN);
  mpfr_set(*wk, h, MPFR_RNDN);
  mpfr_sinh_cosh(msinh, mcosh, kh, MPFR_RNDN);
  mpfr_mul(msinh, msinh, piOver2, MPFR_RNDN);
  mpfr_mul(mcosh, mcosh, piOver2, MPFR_RNDN);
  mpfr_cosh(tmp, msinh, MPFR_RNDN);
  mpfr_sqr(tmp, tmp, MPFR_RNDN);
  mpfr_mul(*wk, *wk, mcosh, MPFR_RNDN);
  mpfr_div(*wk, *wk, tmp, MPFR_RNDN);

  PetscFunctionReturn(0);
}

/* generates n weights spaced h apart in the positive direction */
/* includes "center" n=0 weight as first entry */
#undef __FUNCT__
#define __FUNCT__ "DExkwk"
PetscErrorCode DExkwk(PetscInt n, mpfr_t h, mpfr_t **xk, mpfr_t **wk)
{
  PetscErrorCode ierr;
  PetscInt k;
  mpfr_t kh, msinh, mcosh, piOver2, tmp;
  PetscReal alp, bet;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  mpfr_inits(kh, msinh, mcosh, piOver2, tmp, NULL);

  // init piOver2
  mpfr_const_pi(piOver2, MPFR_RNDN);
  mpfr_mul_d(piOver2, piOver2, 0.5, MPFR_RNDN); flopCount ++;
  
  
  for(k=0; k<n; ++k) {
    /* calculate wk */
    mpfr_set_d(kh, (double) k, MPFR_RNDN);
    mpfr_mul(kh, kh, h, MPFR_RNDN);
    mpfr_set((*wk)[k], h, MPFR_RNDN);
    mpfr_sinh_cosh(msinh, mcosh, kh, MPFR_RNDN);
    mpfr_mul(msinh, msinh, piOver2, MPFR_RNDN);
    mpfr_mul(mcosh, mcosh, piOver2, MPFR_RNDN);
    mpfr_cosh(tmp, msinh, MPFR_RNDN);
    mpfr_sqr(tmp, tmp, MPFR_RNDN);
    mpfr_mul((*wk)[k], (*wk)[k], mcosh, MPFR_RNDN);
    mpfr_div((*wk)[k], (*wk)[k], tmp, MPFR_RNDN); flopCount += 9;

    /* calculate xk */
    mpfr_set_d((*xk)[k], 1.0, MPFR_RNDZ);
    mpfr_cosh(tmp, msinh, MPFR_RNDN);
    mpfr_div((*xk)[k], (*xk)[k], tmp, MPFR_RNDZ);
    mpfr_exp(tmp, msinh, MPFR_RNDN);
    mpfr_div((*xk)[k], (*xk)[k], tmp, MPFR_RNDZ); flopCount += 4;

    
  }
  mpfr_clears(kh, msinh, mcosh, piOver2, tmp, NULL);

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "inte"
PetscErrorCode inte(mpfr_t *x, mpfr_t *val, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  mpfr_mul(*val, *x, *x, MPFR_RNDN);
  mpfr_mul(*val, *val, *x, MPFR_RNDN); flopCount += 2;

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FindStepSize"
PetscErrorCode FindStepSize(PetscInt prec, PetscInt nPts, mpfr_t *step)
{
  PetscErrorCode ierr;
  const PetscReal tol = 1e-25;
  mpfr_t wk, tmp, h, hl, hr, lwk, rwk, swk, slwk, srwk, error;
  PetscInt cont;
  PetscFunctionBegin;

  if(prec < 1) {
    printf("prec should be positive\n");
    PetscFunctionReturn(1);
  }
    
  mpfr_inits(wk, tmp, h, hl, hr, lwk, rwk, swk, slwk, srwk, error, NULL);

  mpfr_set_d(hl, 0.0, MPFR_RNDN);
  mpfr_set_d(hr, 2.0, MPFR_RNDN);

  /* make hr large enough */
  mpfr_set_d(hr, 1.0, MPFR_RNDN);
  cont = 1;
  while(cont == 1) {
    ierr = DEwk(nPts, hr, &rwk);CHKERRQ(ierr);
    mpfr_log10(srwk, rwk, MPFR_RNDN);
    if(mpfr_get_d(srwk, MPFR_RNDN) < -2*prec)
      cont = 0;
    else
      mpfr_mul_d(hr, hr, 2.0, MPFR_RNDN); // make hr larger if value not small
  }
  //printf("hr: %15.15f\n", mpfr_get_d(hr, MPFR_RNDN));
  /* make hl small enough */
  mpfr_mul_d(hl, hr, 0.5, MPFR_RNDN);
  cont = 1;
  while(cont == 1) {
    ierr = DEwk(nPts, hl, &lwk);CHKERRQ(ierr);
    mpfr_log10(slwk, lwk, MPFR_RNDN);
    if(mpfr_get_d(slwk, MPFR_RNDN) > -2*prec)
      cont = 0;
    else
      mpfr_mul_d(hl, hl, 0.5, MPFR_RNDN); // make hl larger if value not small
  }
  //printf("hl: %15.15f\n", mpfr_get_d(hl, MPFR_RNDN));
  //ierr = DEwk(nPts, hl, &lwk);CHKERRQ(ierr);
  //ierr = DEwk(nPts, hr, &rwk);CHKERRQ(ierr);

  cont = 1;
  while(cont == 1) {
    // average hl and hr
    mpfr_add(h, hl, hr, MPFR_RNDN);
    mpfr_mul_d(h, h, 0.5, MPFR_RNDN);
    
    // calculate weight and size
    ierr = DEwk(nPts, h, &wk);CHKERRQ(ierr);
    ierr = DEwk(nPts, hl, &lwk);CHKERRQ(ierr);
    ierr = DEwk(nPts, hr, &rwk);CHKERRQ(ierr);
    mpfr_log10(swk, wk, MPFR_RNDN);
    mpfr_log10(slwk, lwk, MPFR_RNDN);
    mpfr_log10(srwk, rwk, MPFR_RNDN);
    if( mpfr_get_d(swk, MPFR_RNDN) < -2*prec)
      mpfr_set(hr, h, MPFR_RNDN);
    else
      mpfr_set(hl, h, MPFR_RNDN);
    // calculate maximum possible error
    mpfr_sub(error, slwk, srwk, MPFR_RNDN);
    mpfr_abs(error, error, MPFR_RNDN);
    if(mpfr_get_d(error, MPFR_RNDN) < tol)
       cont = 0;
  }
  
  mpfr_set(*step, h, MPFR_RNDN);
  
  mpfr_clears(wk, tmp, h, hl, hr, lwk, rwk, swk, slwk, srwk, error, NULL);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DEQuad"
PetscErrorCode DEQuad(PetscErrorCode (*f)(mpfr_t*,mpfr_t*,void*), mpfr_t a, mpfr_t b, PetscInt prec, PetscInt nPts, PetscReal *integral, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt k;
  PetscInt nSide;
  mpfr_t alpha, beta;
  mpfr_t leftXk, rightXk;
  mpfr_t wk;
  mpfr_t *wkVals;
  mpfr_t *xkVals;
  mpfr_t lx, rx;
  mpfr_t leftFx, rightFx;
  mpfr_t pi2;
  mpfr_t kh, h;
  mpfr_t sum;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  if(nPts % 2 == 0) {
    nSide = nPts/2; flopCount++;
  }
  if(nPts % 2 == 1) {
    nSide = (nPts-1)/2; flopCount++;
  }

  /* initializations */
  mpfr_inits(alpha, beta, leftXk, rightXk, wk, leftFx, rightFx, pi2, kh, h, sum, lx, rx, NULL);
  wkVals = (mpfr_t*) malloc(sizeof(mpfr_t)*(nSide+1));
  xkVals = (mpfr_t*) malloc(sizeof(mpfr_t)*(nSide+1));
  for(k=0; k<nSide+1; ++k)
    mpfr_inits(xkVals[k], wkVals[k], NULL);

  // init alpha, beta
  mpfr_set(alpha, b, MPFR_RNDN);
  mpfr_sub(alpha, alpha, a, MPFR_RNDN);
  mpfr_mul_d(alpha, alpha, 0.5, MPFR_RNDN);
  mpfr_set(beta, b, MPFR_RNDN);
  mpfr_add(beta, beta, a, MPFR_RNDN);
  mpfr_mul_d(beta, beta, 0.5, MPFR_RNDN); flopCount += 4;
  
  /* determine right step size */
  ierr = FindStepSize(prec, nSide, &h); // not including in flop count
  //mpfr_set_d(h, 0.04, MPFR_RNDN);  


  /* calculate abscissas, weights */
  ierr = DExkwk(nSide+1, h, &xkVals, &wkVals);CHKERRQ(ierr);

  /* center term */
  mpfr_set_d(sum, 0.0, MPFR_RNDN);
  mpfr_set(lx, b, MPFR_RNDN);
  mpfr_add(lx, lx, a, MPFR_RNDN);
  mpfr_mul_d(lx, lx, .5, MPFR_RNDN);
  f(&lx, &sum, ctx);
  mpfr_mul(sum, sum, wkVals[0], MPFR_RNDN);
  mpfr_mul(sum, sum, alpha,     MPFR_RNDN); flopCount += 4;
  
  
  for(k=1; k <= nSide; ++k) {
    // adjust xk to [a,b]
    mpfr_sub_d(lx, xkVals[k], 1.0, MPFR_RNDZ);
    mpfr_mul(lx, lx, alpha, MPFR_RNDU);
    mpfr_add(lx, lx, beta, MPFR_RNDU);
    mpfr_d_sub(rx, 1.0, xkVals[k], MPFR_RNDZ);
    mpfr_mul(rx, rx, alpha, MPFR_RNDD);
    mpfr_add(rx, rx, beta , MPFR_RNDD); flopCount += 6;

    double left = mpfr_get_d(lx, MPFR_RNDN);
    double right = mpfr_get_d(rx, MPFR_RNDN);
    double xk = mpfr_get_d(xkVals[k], MPFR_RNDN);

    
    /* calculate function values */
    f(&lx, &leftFx , ctx);
    f(&rx, &rightFx, ctx);

    /* update sum */
    mpfr_mul(leftFx, leftFx, wkVals[k], MPFR_RNDN);
    mpfr_mul(rightFx, rightFx, wkVals[k], MPFR_RNDN);
    mpfr_mul(leftFx, leftFx, alpha, MPFR_RNDN);
    mpfr_mul(rightFx, rightFx, alpha, MPFR_RNDN);
    //double l1 = mpfr_get_d(leftFx, MPFR_RNDN);
    //double l2 = mpfr_get_d(rightFx, MPFR_RNDN);
    //printf("leftval: %15.15f\n", l1);
    //printf("rightval: %15.15f\n", l2);
    mpfr_add(sum, sum, leftFx, MPFR_RNDN);
    mpfr_add(sum, sum, rightFx, MPFR_RNDN); flopCount += 6;
  }
  PetscReal prog = mpfr_get_d(sum, MPFR_RNDN);
  //printf("sum: %15.15f\n", prog);
  
  *integral = prog;
  
  mpfr_clears(alpha, beta, leftXk, rightXk, wk, leftFx, rightFx, pi2, kh, h, sum, lx, rx, NULL);
  for(k=0; k<nSide+1; ++k)
    mpfr_clears(wkVals[k], xkVals[k], NULL);
  PetscFunctionReturn(0);
}
