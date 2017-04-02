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
PetscErrorCode NormConstantIntFixed(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *intVals, PetscReal *normConst)
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
  intVals[0] = integrals[0];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);
  intVals[1] = integrals[1];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);
  intVals[2] = integrals[2];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);
  intVals[3] = integrals[3];
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
  PetscReal solExact;
  PetscReal intExact[4];
  PetscReal solutions[NUM_SOLUTIONS];
  PetscReal errors   [NUM_SOLUTIONS];
  PetscReal intSols[4*NUM_SOLUTIONS];
  PetscReal intErrors[4*NUM_SOLUTIONS];
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

  /* calculate "exact" solution */
  ierr = calcNormalization(&e, n, p, &solExact);CHKERRQ(ierr);
  ierr = calcNormalization2(&e, n, p, intExact, &solExact);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solExact);

  /* calculate approximate solutions and record flops */
  char text[40] = "%d points";
  char sText[40];
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
    ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = NormConstantIntFixed(&e, n, p, prec, POINTS_MIN + POINTS_STEP*i, intSols+4*i, solutions+i);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    printf("new norm constant: %15.15f\n", solutions[i]);
  }

  /* calculate errors */
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    errors[i] = PetscAbsReal((solExact - solutions[i])/solExact);
    intErrors[4*i+0] = PetscAbsReal((intExact[0] - intSols[4*i+0])/intExact[0]);
    intErrors[4*i+1] = PetscAbsReal((intExact[1] - intSols[4*i+1])/intExact[1]);
    intErrors[4*i+2] = PetscAbsReal((intExact[2] - intSols[4*i+2])/intExact[2]);
    intErrors[4*i+3] = PetscAbsReal((intExact[3] - intSols[4*i+3])/intExact[3]);
    printf("errors[%d] = %15.15f\n", i, errors[i]);
  }

  FILE *fp1 = fopen("out/normInt1Prec.txt", "w");
  FILE *fp2 = fopen("out/normInt2Prec.txt", "w");
  FILE *fp3 = fopen("out/normInt3Prec.txt", "w");
  FILE *fp4 = fopen("out/normInt4Prec.txt", "w");
  FILE *fp = fopen("out/normWorkPrec.txt", "w");
  fprintf(fp, "points flops error\n");
  fprintf(fp1, "points error\n");
  fprintf(fp2, "points error\n");
  fprintf(fp3, "points error\n");
  fprintf(fp4, "points error\n");
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, errors[i]);
    fprintf(fp1, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+0]);
    fprintf(fp2, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+1]);
    fprintf(fp3, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+2]);
    fprintf(fp4, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+3]);
  }
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  
  PetscFunctionReturn(0);
}


/* 
NormPlot2 - plots convergence of a normalization
            constant for three different ellipsoidal
	    systems to compare convergence
 */


#undef __FUNCT__
#define __FUNCT__ "NormPlot2"
PetscErrorCode NormPlot2()
{
  const PetscInt NUM_SOLUTIONS = 43;
  const PetscInt POINTS_MIN = 4;
  const PetscInt POINTS_STEP = 2;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  PetscReal solExact;
  PetscReal intExact[4];
  PetscReal solutions[NUM_SOLUTIONS];
  PetscReal errors   [NUM_SOLUTIONS];
  PetscReal intSols[4*NUM_SOLUTIONS];
  PetscReal intErrors[4*NUM_SOLUTIONS];
  PetscReal sol2, sol3, sol4;
  PetscLogEvent flopCounts[NUM_SOLUTIONS];
  EllipsoidalSystem e1, e2, e3;
  EllipsoidalSystem *e;
  PetscErrorCode ierr;
  PetscInt i;
  PetscEventPerfInfo info;
  PetscInt n = 3;
  PetscInt p = 5;
  PetscFunctionBegin;

  
  ierr = initEllipsoidalSystem(&e1, 2.5, 2.0, 1.0);CHKERRQ(ierr);
  ierr = initEllipsoidalSystem(&e2, 2.5, 2.0, 1.5);CHKERRQ(ierr);
  ierr = initEllipsoidalSystem(&e3, 2.5, 2.0, 1.9);CHKERRQ(ierr);

  
  for(PetscInt num=0; num<3; ++num) {
    if(num==0)
      e = &e1;
    else if(num==1)
      e = &e2;
    else
      e = &e3;
    /* calculate "exact" solution */
    ierr = calcNormalization(e, n, p, &solExact);CHKERRQ(ierr);
    ierr = calcNormalization2(e, n, p, intExact, &solExact);CHKERRQ(ierr);
    printf("old norm constant: %15.15f\n", solExact);
    
    /* calculate approximate solutions and record flops */
    char text[40] = "%d points";
    char sText[40];
    for(i=0; i<NUM_SOLUTIONS; ++i) {
      sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
      ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
      ierr = NormConstantIntFixed(e, n, p, prec, POINTS_MIN + POINTS_STEP*i, intSols+4*i, solutions+i);CHKERRQ(ierr);
      ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
      printf("new norm constant: %15.15f\n", solutions[i]);
    }
    
    /* calculate errors */
    for(i=0; i<NUM_SOLUTIONS; ++i) {
      errors[i] = PetscAbsReal((solExact - solutions[i])/solExact);
      intErrors[4*i+0] = PetscAbsReal((intExact[0] - intSols[4*i+0])/intExact[0]);
      intErrors[4*i+1] = PetscAbsReal((intExact[1] - intSols[4*i+1])/intExact[1]);
      intErrors[4*i+2] = PetscAbsReal((intExact[2] - intSols[4*i+2])/intExact[2]);
      intErrors[4*i+3] = PetscAbsReal((intExact[3] - intSols[4*i+3])/intExact[3]);
      printf("errors[%d] = %15.15f\n", i, errors[i]);
    }
    FILE *fp1, *fp2, *fp3, *fp4, *fp;
    if(num==0) {
      fp1 = fopen("out/e1normInt1Prec.txt", "w");
      fp2 = fopen("out/e1normInt2Prec.txt", "w");
      fp3 = fopen("out/e1normInt3Prec.txt", "w");
      fp4 = fopen("out/e1normInt4Prec.txt", "w");
      fp = fopen("out/e1normWorkPrec.txt", "w");
    }
    else if(num==1) {
      fp1 = fopen("out/e2normInt1Prec.txt", "w");
      fp2 = fopen("out/e2normInt2Prec.txt", "w");
      fp3 = fopen("out/e2normInt3Prec.txt", "w");
      fp4 = fopen("out/e2normInt4Prec.txt", "w");
      fp = fopen("out/e2normWorkPrec.txt", "w");
    }
    else {
      fp1 = fopen("out/e3normInt1Prec.txt", "w");
      fp2 = fopen("out/e3normInt2Prec.txt", "w");
      fp3 = fopen("out/e3normInt3Prec.txt", "w");
      fp4 = fopen("out/e3normInt4Prec.txt", "w");
      fp = fopen("out/e3normWorkPrec.txt", "w");
    }
    fprintf(fp, "points flops error\n");
    fprintf(fp1, "points error\n");
    fprintf(fp2, "points error\n");
    fprintf(fp3, "points error\n");
    fprintf(fp4, "points error\n");
    for(i=0; i<NUM_SOLUTIONS; ++i) {
      ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
      fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, errors[i]);
      fprintf(fp1, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+0]);
      fprintf(fp2, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+1]);
      fprintf(fp3, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+2]);
      fprintf(fp4, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, intErrors[4*i+3]);
    }
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
  }
  
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
  ierr = NormPlot2();CHKERRQ(ierr);
  //ierr = NormPlotSE();CHKERRQ(ierr);
  //ierr = NormPlotERF();CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
