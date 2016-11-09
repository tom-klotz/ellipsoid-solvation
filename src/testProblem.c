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
#include <gsl/gsl_rng.h>
#include "constants.h"
#include "testProblem.h"





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
  const double cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/PETSC_PI; /* kcal ang/mol */


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
  const double cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/PETSC_PI; /* kcal ang/mol */


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
  ierr = CalcSolidInteriorHarmonicVec(&e, ellPoints, 2, 3, vals); //testing
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
#define __FUNCT__ "CoulombExact"
PetscErrorCode CoulombExact(PetscReal eps1, Vec srcXYZ, Vec srcMag, Vec tarXYZ, Vec sol)
{
  PetscErrorCode ierr;
  PetscReal srcX, srcY, srcZ;
  PetscReal tarX, tarY, tarZ;
  PetscReal dx2, dy2, dz2;
  PetscReal val, mag;
  const PetscScalar *srcXYZArray;
  const PetscScalar *srcMagArray;
  const PetscScalar *tarXYZArray;
  PetscScalar *solArray;
  PetscInt nSrc, nTar;
  PetscInt i, k;
  PetscFunctionBeginUser;

  ierr = VecGetSize(srcMag, &nSrc);CHKERRQ(ierr);
  ierr = VecGetSize(sol, &nTar);CHKERRQ(ierr);

  ierr = VecGetArrayRead(srcXYZ, &srcXYZArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(tarXYZ, &tarXYZArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(srcMag, &srcMagArray);CHKERRQ(ierr);
  ierr = VecGetArray    (sol   , &solArray   );CHKERRQ(ierr);
  for(i=0; i<nTar; ++i) {
    /* get target X,Y,Z */
    tarX = tarXYZArray[3*i+0];
    tarY = tarXYZArray[3*i+1];
    tarZ = tarXYZArray[3*i+2];
    for(k=0; k<nSrc; ++k) {
      /* get source X,Y,Z */
      srcX = srcXYZArray[3*k+0];
      srcY = srcXYZArray[3*k+1];
      srcZ = srcXYZArray[3*k+2];
      mag  = srcMagArray[k];

      dx2 = (tarX-srcX)*(tarX-srcX);
      dy2 = (tarY-srcY)*(tarY-srcY);
      dz2 = (tarZ-srcZ)*(tarZ-srcZ);
      
      val = mag/(eps1*PetscSqrtReal(dx2 + dy2 + dz2));
      solArray[i] += val;
    }
  }
  ierr = VecRestoreArrayRead(srcXYZ, &srcXYZArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(tarXYZ, &tarXYZArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(srcMag, &srcMagArray);CHKERRQ(ierr);
  ierr = VecRestoreArray    (sol   , &solArray   );CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RandomEllipsoidPoints"
PetscErrorCode RandomEllipsoidPoints(PetscReal a, PetscReal b, PetscReal c, Vec xyz)
{
  PetscErrorCode ierr;
  PetscReal r, theta, phi;
  PetscReal x, y, z;
  const PetscReal theta_min = 0;
  const PetscReal theta_max = 2*PETSC_PI;
  const PetscReal phi_min   = 0;
  const PetscReal phi_max   = PETSC_PI;
  const PetscReal r_min     = 0;
  const PetscReal r_max     = 1;
  PetscScalar *xyzArray;
  PetscInt nvals;
  PetscFunctionBegin;

  /* petsc random object */
  PetscRandom rnd1, rnd2;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rnd1);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rnd2);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd1, .1, 1.0);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd2, .5, 1.0);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rnd1);CHKERRQ(ierr);

  ierr = VecGetSize(xyz, &nvals);CHKERRQ(ierr);
  nvals = nvals/3;

  ierr = VecGetArray(xyz, &xyzArray);CHKERRQ(ierr);
  for(PetscInt i=0; i<nvals; ++i) {
    ierr = PetscRandomGetValue(rnd2, &r);CHKERRQ(ierr);
    ierr = PetscRandomGetValue(rnd1, &theta);CHKERRQ(ierr);
    ierr = PetscRandomGetValue(rnd1, &phi);CHKERRQ(ierr);

    //r = r_min + (r_max-r_min)*i/(nvals-1);
    //theta = theta_min + (theta_max-theta_min)*i/(nvals-1);
    //phi = phi_min + (phi_max-phi_min)*i/(nvals-1);
    r = r_min + (r_max-r_min)*r;
    theta = theta_min + (theta_max-theta_min)*theta;
    phi = phi_min + (phi_max-phi_min)*phi;
    
    x = a*r*PetscCosReal(theta)*PetscSinReal(phi);
    y = b*r*PetscSinReal(theta)*PetscSinReal(phi);
    z = c*r*PetscCosReal(theta);

    xyzArray[3*i+0] = x;
    xyzArray[3*i+1] = y;
    xyzArray[3*i+2] = z;
    if(nvals==1) {
      xyzArray[3*i+0] = 0;
      xyzArray[3*i+1] = 0;
      xyzArray[3*i+2] = .1;
    }
  }
  ierr = VecRestoreArray(xyz, &xyzArray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RandomMag"
PetscErrorCode RandomMag(Vec mags)
{
  PetscErrorCode ierr;
  PetscInt nMag, i;
  PetscReal val;
  PetscRandom rnd;
  PetscScalar *magsArray;
  PetscFunctionBegin;

  ierr = VecGetSize(mags, &nMag);CHKERRQ(ierr);
  
  /* petsc random object */
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rnd);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd, -1.0, 1.0);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rnd);CHKERRQ(ierr);
  
  ierr = VecGetArray(mags, &magsArray);CHKERRQ(ierr);
  for(i=0; i<nMag; ++i) {
    ierr = PetscRandomGetValue(rnd, &val);CHKERRQ(ierr);
    magsArray[i] = val;
  }
  ierr = VecRestoreArray(mags, &magsArray);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WorkPrecExample"
PetscErrorCode WorkPrecExample(PetscInt Nmax)
{
  PetscErrorCode ierr;
  
  const PetscInt NUM_SOLUTIONS = 10;
  const PetscInt CHARGE_START = 100000;
  const PetscInt CHARGE_INC   = 50000;
  Vec srcXYZ[NUM_SOLUTIONS], srcMag[NUM_SOLUTIONS], solution[NUM_SOLUTIONS];
  PetscInt i;
  PetscInt chargeNums[NUM_SOLUTIONS];
  PetscLogEvent events[NUM_SOLUTIONS];
  PetscInt charges;

  char tText[30] = "Free energy with %d charges";
  char eText[30];
  PetscInt tempNum;

  EllipsoidalSystem e;
  PetscFunctionBegin;

  const PetscReal eps1 = 4.0;
  const PetscReal eps2 = 80.0;

  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  ierr = initEllipsoidalSystem(&e, a, b, c);CHKERRQ(ierr);
  
  /* create the vector with charge numbers to use */
  for(i=0; i < NUM_SOLUTIONS; ++i) {
    charges = CHARGE_START + CHARGE_INC*i;
    chargeNums[i] = charges;
    /* initialize xyz and solution vectors */
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3*charges, srcXYZ+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, charges, srcMag+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, charges, solution+i);CHKERRQ(ierr);
    /* generate random points and magnitudes */
    ierr = RandomEllipsoidPoints(a, b, c, srcXYZ[i]);CHKERRQ(ierr);
    ierr = RandomMag(srcMag[i]);CHKERRQ(ierr);

    /* register events for flop counting */
    sprintf(eText, tText, charges);
    ierr = PetscLogEventRegister((const char*) eText, 0, events+i);CHKERRQ(ierr);
  }
  

  printf("wow we here\n");
  for(i=0; i < NUM_SOLUTIONS; ++i) {
    printf("calculating olution %d/%d\n", i+1, NUM_SOLUTIONS);
    ierr = PetscLogEventBegin(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = CalcEllipsoidFreeEnergy(&e, eps1, eps2, chargeNums[i], srcXYZ[i], srcMag[i], 1e-5, Nmax, solution[i]);
    ierr = PetscLogEventEnd(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
  }
  //ierr = EasyExample(Nmax, nSrc, nx, xl, xr, ny, yl, yr, zConst);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitGrid"
PetscErrorCode InitGrid(PetscInt nx, PetscReal xl, PetscReal xr, PetscInt ny, PetscReal yl, PetscReal yr, PetscReal zConst, Vec *xyz)
{
  PetscErrorCode ierr;
  PetscInt i, j;
  PetscReal xh, yh;
  PetscInt index;
  PetscScalar *xyzArray;
  PetscFunctionBegin;

  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nx*ny, xyz);CHKERRQ(ierr);

  xh = (xr - xl)/(nx - 1);
  yh = (yr - yl)/(ny - 1);

  ierr = VecGetArray(*xyz, &xyzArray);CHKERRQ(ierr);
  for(i=0; i < nx; ++i) {
    for(j=0; j<ny; ++j) {
      index = 3*(i*ny + j);
      xyzArray[index+0] = xl + xh*i;
      xyzArray[index+1] = yl + yh*j;
      xyzArray[index+2] = zConst;
    }
  }
  ierr = VecRestoreArray(*xyz, &xyzArray);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GridSolution"
/*
  Example is of a some charges
*/
PetscErrorCode GridSolution(PetscInt Nmax, PetscInt nSrc, PetscInt nx, PetscReal xl, PetscReal xr, PetscInt ny, PetscReal yl, PetscReal yr, PetscReal zConst, Vec *xvals, Vec *yvals, Vec *sol)
{
  PetscErrorCode ierr;
  Vec sourceXYZ, sourceMag, tarXYZ, solution;
  PetscReal xh, yh;
  PetscReal xc, yc;
  PetscReal randVal;
  PetscInt ind;
  PetscFunctionBegin;

  const PetscReal eps1 = 4.0;
  const PetscReal eps2 = 4.0;

  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  // initialize XYZ and solution vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nSrc , &sourceXYZ);CHKERRQ(ierr);
  //ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nx*ny, &tarXYZ);   CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , &solution); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSrc   , &sourceMag);CHKERRQ(ierr);

  xh = (xr - xl) / (nx - 1);
  yh = (yr - yl) / (ny - 1);


  ierr = RandomEllipsoidPoints(a, b, c, sourceXYZ);CHKERRQ(ierr);

  // populate source vec with points
  for(PetscInt k=0; k<nSrc; ++k) {
    ind = 3*k;
    ierr = VecSetValue(sourceMag, k  , 1.0      , INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(sourceMag);CHKERRQ(ierr); ierr = VecAssemblyEnd(sourceMag);CHKERRQ(ierr);

  
  ierr = InitGrid(nx, xl, xr, ny, yl, yr, zConst, &tarXYZ);

  // calculate the solvation potential
  ierr = CalcEllipsoidTester( a , b , c,
			      eps1, eps2,
			      nSrc, sourceXYZ, sourceMag,
			      nx*ny, tarXYZ,
			      Nmax, solution); CHKERRQ(ierr);
  
  /*
  FILE *fp = fopen("solvPlot.txt", "w");
  fprintf(fp, "xh: %4.4e\n", xh);
  fprintf(fp, "yh: %4.4e\n", yh);
  const PetscScalar *solutionArray;
  ierr = VecGetArrayRead(solution, &solutionArray);CHKERRQ(ierr);
  PetscInt index = 0;
  for(PetscInt i=0; i<nx; ++i) {
    for(PetscInt j=0; j<ny; ++j) {
      index = i*ny + j;
      fprintf(fp, "%4.4f ", solutionArray[index]);
    }
    fprintf(fp, "\n");
  }
  ierr = VecRestoreArrayRead(solution, &solutionArray);CHKERRQ(ierr);
  */
  if(xvals != NULL) {
    ierr = VecCreateSeq(PETSC_COMM_WORLD, nx, xvals);CHKERRQ(ierr);
    for(PetscInt i=0; i<nx; ++i)
      ierr = VecSetValue(*xvals, i, xl + i*xh, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*xvals);CHKERRQ(ierr); ierr = VecAssemblyEnd(*xvals);CHKERRQ(ierr);
  }
  if(yvals != NULL) {
    ierr = VecCreateSeq(PETSC_COMM_WORLD, nx, yvals);CHKERRQ(ierr);
    for(PetscInt i=0; i<ny; ++i)
      ierr = VecSetValue(*yvals, i, yl + i*yh, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*yvals);CHKERRQ(ierr); ierr = VecAssemblyEnd(*yvals);CHKERRQ(ierr);
  }
  if(sol != NULL) {
    ierr = VecDuplicate(solution, sol);CHKERRQ(ierr);
    ierr = VecCopy(solution, *sol);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&sourceXYZ);CHKERRQ(ierr);
  ierr = VecDestroy(&tarXYZ);CHKERRQ(ierr);
  ierr = VecDestroy(&solution);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GridAnimation"
PetscErrorCode GridAnimation()
{
  const PetscInt NUM_SOLUTIONS = 8;
  PetscErrorCode ierr;
  Vec sourceXYZ, sourceMag, tarXYZ;
  Vec solution[NUM_SOLUTIONS], solutionEx[NUM_SOLUTIONS], errorVec[NUM_SOLUTIONS];
  Vec xOut, yOut, zOut;
  PetscReal xh, yh;
  PetscReal xc, yc;
  PetscReal randVal;
  PetscInt ind;
  PetscFunctionBegin;

  PetscReal eps1 = 4.0;
  PetscReal eps2 = 4.0;
  PetscInt nSrc  = 10;
  PetscInt nx    = 10;
  PetscReal xl   = -4.8;
  PetscReal xr   = 4.8;
  PetscInt ny    = 10;
  PetscReal yl   = -4.2;
  PetscReal yr   =  4.2;
  PetscReal zConst = .1;
  
  
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  // initialize source XYZ, mag
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nSrc , &sourceXYZ);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSrc   , &sourceMag);CHKERRQ(ierr);
  // initialize solution vectors
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , &solution[i]); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , &solutionEx[i]); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , &errorVec[i]); CHKERRQ(ierr);
    
  }


  ierr = RandomEllipsoidPoints(a, b, -.1, sourceXYZ);CHKERRQ(ierr); // on z=0
  
  ierr = RandomMag(sourceMag); CHKERRQ(ierr);
  
  ierr = InitGrid(nx, xl, xr, ny, yl, yr, 0, &tarXYZ); //zConst = 0


  PetscReal err;
  // calculate the solvation potential
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    printf("wow\n");
    ierr = CoulombExact(eps1, sourceXYZ, sourceMag, tarXYZ, solutionEx[i]);
    ierr = CalcEllipsoidTester( a , b , c,
				eps1, eps2,
				nSrc, sourceXYZ, sourceMag,
				nx*ny, tarXYZ,
				i, solution[i]); CHKERRQ(ierr);
    ierr = VecCopy(solutionEx[i], errorVec[i]);CHKERRQ(ierr);
    ierr = VecAXPY(errorVec[i], -1.0, solution[i]);
    ierr = VecNorm(errorVec[i], NORM_2, &err);
    printf("\n\n");
    printf("the error on %d is %15.15f\n", i, err);
    //ierr = VecView(solutionEx[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    //ierr = VecView(solution[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = VecView(errorVec[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    printf("\n\n");
  }

  char fname[20] = "out/sol%d.txt";
  char fnameS[20];
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    sprintf(fnameS, fname, i);
    ierr = WriteToFile((const char*) fnameS, ny, solution[i]);CHKERRQ(ierr);
  }

  ierr = VecCreateSeq(PETSC_COMM_SELF, nx, &xOut);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ny, &yOut);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, 1, &zOut);CHKERRQ(ierr);
  xh = (xr - xl) / (nx - 1);
  yh = (yr - yl) / (ny - 1);
  for(PetscInt i=0; i<nx; ++i) {
    ierr = VecSetValue(xOut, i, xl + xh*i, INSERT_VALUES);CHKERRQ(ierr);
  }
  for(PetscInt i=0; i<ny; ++i) {
    ierr = VecSetValue(yOut, i, yl + yh*i, INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecSetValue(zOut, 0, zConst, INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(xOut);CHKERRQ(ierr); ierr = VecAssemblyEnd(xOut);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(yOut);CHKERRQ(ierr); ierr = VecAssemblyEnd(yOut);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(zOut);CHKERRQ(ierr); ierr = VecAssemblyEnd(zOut);CHKERRQ(ierr);
  
  ierr = WriteToFile((const char*) "out/xVals.txt", 1, xOut);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/yVals.txt", 1, yOut);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/zVals.txt", 1, zOut);CHKERRQ(ierr);

  /* output source xyz values */
  ierr = WriteToFile((const char*) "out/chargeXYZ.txt", 3, sourceXYZ);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/chargeMag.txt", 1, sourceMag);CHKERRQ(ierr);
  FILE *fp = fopen("out/otherinfo.txt", "w");
  fprintf(fp, "%6.6f %6.6f %6.6f\n", a, b, c);
  fclose(fp);
  


  ierr = VecDestroy(&sourceXYZ);CHKERRQ(ierr);
  ierr = VecDestroy(&tarXYZ);CHKERRQ(ierr);
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    ierr = VecDestroy(solution+i);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&sourceMag);CHKERRQ(ierr);
  ierr = VecDestroy(&xOut);CHKERRQ(ierr);
  ierr = VecDestroy(&yOut);CHKERRQ(ierr);
  ierr = VecDestroy(&zOut);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "numChargesPlot"
PetscErrorCode numChargesPlot(PetscInt nMin, PetscInt nMax, PetscInt nStep)
{
  PetscErrorCode ierr;
  PetscInt Nmax = 4;
  //PetscInt nSrc = 20;
  PetscInt nx   = 7;
  PetscReal xl  = -4.6;
  PetscReal xr  = 4.6;
  PetscInt ny   = 7;
  PetscReal yl  = -4.6;
  PetscReal yr  = 4.6;
  PetscReal zConst = .1;
  PetscInt nSource;
  PetscLogStage *chargeStages;
  PetscStageLog stageLog;
  PetscFunctionBegin;

  PetscInt count = 0;
  for(nSource = nMin; nSource <= nMax; nSource += nStep) {
    count++;
  }
  ierr = PetscMalloc1(sizeof(PetscLogStage)*count, &chargeStages);
  char bStr[40] = "Num Charges: %d";
  char indvStr[40];
  count = 0;
  for(nSource = nMin; nSource <= nMax; nSource += nStep) {
    sprintf(indvStr, bStr, nSource);
    ierr = PetscLogStageRegister(indvStr, chargeStages+count);CHKERRQ(ierr);
    count++;
  }
  count = 0;
  for(nSource = nMin; nSource <= nMax; nSource += nStep) {
    printf("wow\n");
    ierr = PetscLogStagePush(chargeStages[count]);CHKERRQ(ierr);
    ierr = GridSolution(Nmax, nSource, nx, xl, xr, ny, yl, yr, zConst, NULL, NULL, NULL);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
    count++;
  }

  ierr = PetscLogGetStageLog(&stageLog);CHKERRQ(ierr);
  count = 0;
  for(nSource = nMin; nSource <= nMax; nSource += nStep) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\nflops = %4.4e\n", nSource, stageLog->stageInfo[chargeStages[count]].perfInfo.flops);CHKERRQ(ierr);
    count++;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "WriteToFile"
PetscErrorCode WriteToFile(const char *fname, PetscInt rowsize, Vec values)
{
  PetscErrorCode ierr;
  PetscInt k;
  PetscInt nPts;
  const PetscScalar *valuesArray;
  PetscFunctionBegin;

  ierr = VecGetSize(values, &nPts);CHKERRQ(ierr);

  ierr = VecGetArrayRead(values, &valuesArray);CHKERRQ(ierr);
  FILE *fp = fopen(fname, "w");
  for(k=0; k < nPts; ++k) {
    fprintf(fp, "%4.4f ", valuesArray[k]);
    if((k+1)%rowsize == 0)
      fprintf(fp, "\n");
  }
  fclose(fp);
  ierr = VecRestoreArrayRead(values, &valuesArray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolutionAnimation"
PetscErrorCode SolutionAnimation()
{
  PetscErrorCode ierr;
  Vec x1, y1, s0, s1, s2, s3;
  PetscFunctionBegin;
  PetscInt nSrc = 10;
  PetscInt nx   = 15;
  PetscReal xl  = -4.6;
  PetscReal xr  = 4.6;
  PetscInt ny   = 15;
  PetscReal yl  = -4.6;
  PetscReal yr  = 4.6;
  PetscReal zConst = .1;

  ierr = GridSolution(0, nSrc, nx, xl, xr, ny, yl, yr, zConst, &x1, &y1, &s0);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/sol0.txt", ny, s0);CHKERRQ(ierr);
  ierr = GridSolution(1, nSrc, nx, xl, xr, ny, yl, yr, zConst, &x1, &y1, &s1);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/sol1.txt", ny, s1);CHKERRQ(ierr);
  ierr = GridSolution(2, nSrc, nx, xl, xr, ny, yl, yr, zConst, &x1, &y1, &s2);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/sol2.txt", ny, s2);CHKERRQ(ierr);
  ierr = GridSolution(3, nSrc, nx, xl, xr, ny, yl, yr, zConst, &x1, &y1, &s3);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/sol3.txt", ny, s3);CHKERRQ(ierr);

  ierr = WriteToFile((const char*) "out/xVals.txt", 1, x1);CHKERRQ(ierr);
  ierr = WriteToFile((const char*) "out/yVals.txt", 1, y1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main( int argc, char **argv )
{

  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;
  //initialize Petsc
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);  

  
  //testSphericalCoordinates();
  //runTest1();
  //testLegendre();
  //runTest4();
  //test4sphere();
  //RunArg();
  //ierr = RunArgTester(); CHKERRQ(ierr);
  //PetscErrorCode GridSolution(PetscInt Nmax, PetscInt nSrc, PetscInt nx, PetscReal xl, PetscReal xr, PetscInt ny, PetscReal yl, PetscReal yr, PetscReal zConst)
  
  PetscInt Nmax = 4;
  PetscInt nSrc = 25;
  PetscInt nx   = 15;
  PetscReal xl  = -5.23;
  PetscReal xr  = 5.23;
  PetscInt ny   = 15;
  PetscReal yl  = -4.45;
  PetscReal yr  = 4.45;
  PetscReal zConst = .1;
  PetscReal eps1 = 4.0;
  PetscReal eps2 = 80.0;
  
  //ierr = GridSolution(Nmax, nSrc, nx, xl, xr, ny, yl, yr, zConst, NULL, NULL);CHKERRQ(ierr);
  //ierr = SolutionAnimation(); //<---- total crap
  ierr = GridAnimation();CHKERRQ(ierr);
  //ierr = WorkPrecExample(Nmax);

  //ierr = numChargesPlot(100, 500, 100);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
