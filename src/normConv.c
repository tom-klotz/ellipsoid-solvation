#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "ellipsoid/ellipsoid.h"
#include "tanhsinh.h"


#undef __FUNCT__
#define __FUNCT__ "NormConstantIntFixed"
PetscErrorCode NormConstantIntFixed(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *normConst)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;

  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };

  PetscReal integrals[4];

  mpfr_t mpfrzero;
  mpfr_t mpfrone;
  mpfr_inits(mpfrzero, mpfrone, NULL);
  mpfr_set_d(mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(mpfrone, 1.0, MPFR_RNDN);
  
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+0, &ctx1);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);

  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]); flopCount += 4;

  mpfr_clears(mpfrzero, mpfrone, NULL);
  
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "NormConstantIntFixedSE"
PetscErrorCode NormConstantIntFixedSE(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *normConst)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;

  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };

  PetscReal integrals[4];

  mpfr_t mpfrzero;
  mpfr_t mpfrone;
  mpfr_inits(mpfrzero, mpfrone, NULL);
  mpfr_set_d(mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(mpfrone, 1.0, MPFR_RNDN);
  
  ierr = SEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+0, &ctx1);

  ierr = SEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);

  ierr = SEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);

  ierr = SEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);

  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]); flopCount += 4;

  mpfr_clears(mpfrzero, mpfrone, NULL);
  
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NormConstantIntFixedERF"
PetscErrorCode NormConstantIntFixedERF(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *normConst)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;

  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };

  PetscReal integrals[4];

  mpfr_t mpfrzero;
  mpfr_t mpfrone;
  mpfr_inits(mpfrzero, mpfrone, NULL);
  mpfr_set_d(mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(mpfrone, 1.0, MPFR_RNDN);
  
  ierr = ERFQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+0, &ctx1);
  //ierr = ERFQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) inte, mpfrzero, mpfrone, prec, nPts, integrals+0, &ctx1);
  //printf("inte: %15.15f\n", integrals[0]);
  ierr = ERFQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);

  ierr = ERFQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);

  ierr = ERFQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);

  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]); flopCount += 4;

  mpfr_clears(mpfrzero, mpfrone, NULL);
  
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "NormPlot"
PetscErrorCode NormPlot()
{
  const PetscInt NUM_SOLUTIONS = 43;
  const PetscInt POINTS_MIN = 4;
  const PetscInt POINTS_STEP = 2;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  PetscReal solOld;
  PetscReal solutions[NUM_SOLUTIONS];
  PetscReal errors   [NUM_SOLUTIONS];
  PetscReal sol2, sol3, sol4;
  PetscLogEvent flopCounts[NUM_SOLUTIONS];
  EllipsoidalSystem e;
  PetscErrorCode ierr;
  PetscInt i;
  PetscEventPerfInfo info;
  PetscInt n = 3;
  PetscInt p = 5;
  PetscFunctionBegin;

  
  initEllipsoidalSystem(&e, a, b, c);

  /* calculate approximate solutions and record flops */
  char text[40] = "%d points";
  char sText[40];
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
    ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = NormConstantIntFixed(&e, n, p, prec, POINTS_MIN + POINTS_STEP*i, solutions+i);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    printf("new norm constant: %15.15f\n", solutions[i]);
  }

  /* calculate "exact" solution */
  ierr = calcNormalization(&e, n, p, &solOld);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solOld);


  /* calculate errors */
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    errors[i] = PetscAbsReal((solOld - solutions[i])/solOld);
    printf("errors[%d] = %15.15f\n", i, errors[i]);
  }


  FILE *fp = fopen("out/normWorkPrec.txt", "w");
  fprintf(fp, "points flops error\n");
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, errors[i]);
  }
  fclose(fp);

  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NormPlotSE"
PetscErrorCode NormPlotSE()
{
  const PetscInt NUM_SOLUTIONS = 60;
  const PetscInt POINTS_MIN = 4;
  const PetscInt POINTS_STEP = 2;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  PetscReal solOld;
  PetscReal solutions[NUM_SOLUTIONS];
  PetscReal errors   [NUM_SOLUTIONS];
  PetscReal sol2, sol3, sol4;
  PetscLogEvent flopCounts[NUM_SOLUTIONS];
  EllipsoidalSystem e;
  PetscErrorCode ierr;
  PetscInt i;
  PetscEventPerfInfo info;
  PetscInt n = 3;
  PetscInt p = 5;
  PetscFunctionBegin;

  
  initEllipsoidalSystem(&e, a, b, c);

  /* calculate approximate solutions and record flops */
  char text[40] = "%d points";
  char sText[40];
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
    ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = NormConstantIntFixedSE(&e, n, p, prec, POINTS_MIN + POINTS_STEP*i, solutions+i);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    printf("new norm constant (N=%d): %15.15f\n", POINTS_MIN+POINTS_STEP*i, solutions[i]);
  }

  /* calculate "exact" solution */
  ierr = calcNormalization(&e, n, p, &solOld);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solOld);


  /* calculate errors */
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    errors[i] = PetscAbsReal((solOld - solutions[i])/solOld);
    printf("errors[%d] = %15.15f\n", i, errors[i]);
  }


  FILE *fp = fopen("out/normWorkPrecSE.txt", "w");
  fprintf(fp, "points flops error\n");
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, errors[i]);
  }
  fclose(fp);

  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "NormPlotERF"
PetscErrorCode NormPlotERF()
{
  const PetscInt NUM_SOLUTIONS = 50;
  const PetscInt POINTS_MIN = 4;
  const PetscInt POINTS_STEP = 2;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  PetscReal solOld;
  PetscReal solutions[NUM_SOLUTIONS];
  PetscReal errors   [NUM_SOLUTIONS];
  PetscReal sol2, sol3, sol4;
  PetscLogEvent flopCounts[NUM_SOLUTIONS];
  EllipsoidalSystem e;
  PetscErrorCode ierr;
  PetscInt i;
  PetscEventPerfInfo info;
  PetscInt n = 5;
  PetscInt p = 3;
  PetscFunctionBegin;

  
  initEllipsoidalSystem(&e, a, b, c);

  /* calculate "exact" solution */
  ierr = calcNormalization(&e, n, p, &solOld);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solOld);

  /* calculate approximate solutions and record flops */
  char text[40] = "%d points";
  char sText[40];
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
    ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = NormConstantIntFixedERF(&e, n, p, prec, POINTS_MIN + POINTS_STEP*i, solutions+i);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    printf("new norm constant (N=%d): %15.15f\n", POINTS_MIN+POINTS_STEP*i, solutions[i]);
  }




  /* calculate errors */
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    errors[i] = PetscAbsReal((solOld - solutions[i])/solOld);
    printf("errors[%d] = %15.15f\n", i, errors[i]);
  }


  FILE *fp = fopen("out/normWorkPrecERF.txt", "w");
  fprintf(fp, "points flops error\n");
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, errors[i]);
  }
  fclose(fp);

  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{

  const PetscInt NUM_DIGITS = 128;
  const PetscInt PRECISION  = 32;
  mpfr_t a, b;
  mpfr_t val;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  const PetscInt nPts = 10;

  /* initializations */
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);
  
  mpfr_set_default_prec(4*NUM_DIGITS);
  mpfr_inits(a, b, val, NULL);
  mpfr_set_d(a, 2.0, MPFR_RNDN);
  mpfr_set_d(b, 3.0, MPFR_RNDN);
  
  /*
  int fac = 1;
  for(int i=0; i<10; ++i) {
    printf("\n########################################################\n");
    ierr = DEQuad(inte, a, b, PRECISION, fac*nPts, &val, NULL);CHKERRQ(ierr);
    fac *= 2;
  }
  */
  
  //ierr = NormPlot();CHKERRQ(ierr);
  ierr = NormPlotSE();CHKERRQ(ierr);
  //ierr = NormPlotERF();CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
