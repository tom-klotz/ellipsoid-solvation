#include <math.h>
#include <stdio.h>
#include <mpfr.h>
#include <petsc.h>
#include <stdlib.h>
#include <string.h>

#include "../ellipsoid/ellipsoid.h"

#undef __FUNCT__
#define __FUNCT__ "TestNormalizationMPFR"
PetscErrorCode TestNormalizationMPFR()
{
  PetscErrorCode ierr;
  EllipsoidalSystem e;
  const PetscInt prec = 128;
  PetscReal xA = 3.0;
  PetscReal yB = 2.0;
  PetscReal zC = 1.0;
  mpfr_t a, b, c;
  mpfr_t a2, b2, c2;
  mpfr_t a3, b3, c3;
  mpfr_t temp1, temp2, temp3;
  mpfr_t pi;

  
  PetscFunctionBegin;
  mpfr_set_default_prec(4*prec);

  
  mpfr_inits(a, b, c, NULL);
  mpfr_inits(a2, b2, c2, NULL);
  mpfr_inits(a3, b3, c3, NULL);
  mpfr_inits(temp1, temp2, temp3, pi, NULL);
  
  mpfr_set_d(a, xA, MPFR_RNDN);
  mpfr_set_d(b, yB, MPFR_RNDN);
  mpfr_set_d(c, zC, MPFR_RNDN);


  //init ellipsoidal system
  printf("initEllipsoidalSystem...\n");
  ierr = initEllipsoidalSystem(&e, xA, yB, zC, prec);CHKERRQ(ierr);
  printf("done.\n");

  
  mpfr_mul(a2, a, a, MPFR_RNDN);
  mpfr_mul(a3, a2, a, MPFR_RNDN);
  mpfr_mul(b2, b, b, MPFR_RNDN);
  mpfr_mul(b3, b2, b, MPFR_RNDN);
  mpfr_mul(c2, c, c, MPFR_RNDN);
  mpfr_mul(c3, c2, c, MPFR_RNDN);

  mpfr_t firstTerm, secondTerm, LambdaD, LambdaDprime;
  mpfr_inits(firstTerm, secondTerm, LambdaD, LambdaDprime, NULL); 

  //Dassios (B14)
  //firstTerm = (a*a + b*b + c*c)/3.0;
  mpfr_add(firstTerm, a2, b2, MPFR_RNDN);
  mpfr_add(firstTerm, firstTerm, c2, MPFR_RNDN);
  mpfr_div_d(firstTerm, firstTerm, 3.0, MPFR_RNDN);

  //double secondTerm = sqrt((a*a*a*a - b*b*c*c) + (b*b*b*b - a*a*c*c) + (c*c*c*c - a*a*b*b))/3.0;
  mpfr_mul(temp1, a3, a, MPFR_RNDN);
  mpfr_mul(temp2, b2, c2, MPFR_RNDN);
  mpfr_sub(temp1, temp1, temp2, MPFR_RNDN);
  mpfr_mul(temp2, b3, b, MPFR_RNDN);
  mpfr_mul(temp3, a2, c2, MPFR_RNDN);
  mpfr_sub(temp2, temp2, temp3, MPFR_RNDN);
  mpfr_add(temp1, temp1, temp2, MPFR_RNDN);
  mpfr_mul(temp2, c3, c, MPFR_RNDN);
  mpfr_mul(temp3, a2, b2, MPFR_RNDN);
  mpfr_sub(temp2, temp2, temp3, MPFR_RNDN);
  mpfr_add(temp1, temp1, temp2, MPFR_RNDN);
  mpfr_sqrt(secondTerm, temp1, MPFR_RNDN);

  //double LambdaD = firstTerm + secondTerm;
  mpfr_add(LambdaD, firstTerm, secondTerm, MPFR_RNDN);

  //double LambdaDprime = firstTerm - secondTerm;
  mpfr_sub(LambdaDprime, firstTerm, secondTerm, MPFR_RNDN);

  mpfr_t hx, hy, hz;
  mpfr_inits(hx, hy, hz, NULL);

  //Dassios (B16)-(B20)
  //double hx = sqrt(b*b - c*c);
  mpfr_sub(temp1, b2, c2, MPFR_RNDN);
  mpfr_sqrt(hx, temp1, MPFR_RNDN);

  //double hy = e.k;
  mpfr_set(hy, e.hp_k, MPFR_RNDN);
  //double hz = e.h;
  mpfr_set(hz, e.hp_h, MPFR_RNDN);

  //set pi
  mpfr_const_pi(pi, MPFR_RNDN);

  mpfr_t analytic[9];
  mpfr_t approx[9];
  for(PetscInt k=0; k < 9; ++k)
    mpfr_inits(analytic[k], approx[k], NULL);

  //analytic[0] = 4*PETSC_PI
  mpfr_mul_d(analytic[0], pi, 4.0, MPFR_RNDN);
  printf("analytic[0] = %4.4e\n", mpfr_get_d(analytic[0], MPFR_RNDN));
  
  //analytic[1] = 4*PETSC_PI/3 * hy*hy*hz*hz
  mpfr_div_d(temp1, pi, 3.0, MPFR_RNDN);
  mpfr_mul_d(temp1, temp1, 4.0, MPFR_RNDN);
  mpfr_mul(temp2, hy, hy, MPFR_RNDN);
  mpfr_mul(temp2, temp2, hz, MPFR_RNDN);
  mpfr_mul(temp2, temp2, hz, MPFR_RNDN);
  mpfr_mul(analytic[1], temp1, temp2, MPFR_RNDN);
  printf("analytic[1] = %4.4e\n", mpfr_get_d(analytic[1], MPFR_RNDN));



  ierr = calcNormalizationMPFR(&e, 0, 0, approx+0);CHKERRQ(ierr);
  ierr = calcNormalizationMPFR(&e, 1, 0, approx+1);CHKERRQ(ierr);
  mpfr_t err, l;
  mpfr_inits(err, l, NULL);
  mpfr_set_d(l, 1.0, MPFR_RNDN);


  mpfr_sub(err, approx[1], analytic[1], MPFR_RNDN);
  mpfr_div(err, err, approx[1], MPFR_RNDN);
  mpfr_abs(err, err, MPFR_RNDN);
  mpfr_log10(err, err, MPFR_RNDN);
  printf("the error is %2.2f\n", mpfr_get_d(err, MPFR_RNDN));
  mpfr_clears(err, l, NULL);
  

  for(PetscInt k=0; k < 9; ++k)
    mpfr_clears(analytic[k], approx[k], NULL);
  mpfr_clears(hx, hy, hz, NULL);
  mpfr_clears(temp1, temp2, temp3, pi, NULL);
  mpfr_clears(a3, b3, c3, NULL);
  mpfr_clears(a2, b2, c2, NULL);
  mpfr_clears(a, b, c, NULL);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT __ "main"
PetscErrorCode main(int argc, char **argv)
{
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

  ierr = TestNormalizationMPFR();
  
  ierr = PetscFinalize();

}
