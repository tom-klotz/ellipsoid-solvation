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
#include "constants.h"
#include "testProblem.h"


//compares ellipsoidal and spherical
//has ellipsoidal distribution of charges just inside ellipsoidal boundary
//computes the solution at two point, one inside smallest enclosing sphere where spherical diverges
//other solution point is outside of sphere and ellipse to compare convergence
#undef __FUNCT__
#define __FUNCT__ "EllVsSphConvergence"
PetscErrorCode EllVsSphConvergence()
{
  PetscErrorCode ierr;
  //semi-axes of ellipsoid
  PetscReal a = 3.0;
  PetscReal b = 2.5;
  PetscReal c = 1.2;
  //sphere radius
  PetscReal sphRad = 3;
  //dielecric permitivities both set to 1
  const PetscReal eps1 = 1.0;
  const PetscReal eps2 = 1.0;
  //charges placed on ellipsoid with sem-axes smaller by epsilon
  PetscReal epsilon = .5;
  //interior charge XYZ and magnitude
  const PetscInt SQRT_NUM_CHARGES = 3;
  Vec chargeXYZ;
  Vec chargeMag;
  //maximum expansion order for test
  const PetscInt MAX_N = 20;
  //solutions vector
  const PetscInt NUM_SOL_PTS = 2;
  Vec ellSolutions[MAX_N], sphSolutions[MAX_N];
  Vec ellError[MAX_N], sphError[MAX_N];
  Vec exactSolution;
  Vec solXYZ;
  //solution point inside of sphere
  PetscReal sol1X = .1;
  PetscReal sol1Y = .1;
  PetscReal sol1Z = 1.5;
  //solution point outside of sphere
  PetscReal sol2X = 2.1;
  PetscReal sol2Y = 2.6;
  PetscReal sol2Z = 1;
  //pointer for accessing petsc vectors
  PetscScalar *vecPt;
  PetscFunctionBegin;

  //create, set charge magnitudes to 1
  ierr = VecCreateSeq(PETSC_COMM_SELF, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, &chargeMag);CHKERRQ(ierr);
  ierr = VecSet(chargeMag, 1.0);CHKERRQ(ierr);

  //create, set charge XYZ values to be on a-eps,b-eps,c-eps ellipsoid
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, &chargeXYZ);CHKERRQ(ierr);
  ierr = VecGetArray(chargeXYZ, &vecPt);CHKERRQ(ierr);
  PetscReal pu, pv;
  for(PetscInt i=0; i<SQRT_NUM_CHARGES; ++i) {
    for(PetscInt j=0; j<SQRT_NUM_CHARGES; ++j) {
      pu = i*2.0*PETSC_PI/(SQRT_NUM_CHARGES-1);
      pv = j*PETSC_PI/(SQRT_NUM_CHARGES-1);
      PetscInt ind = i*SQRT_NUM_CHARGES + j;
      //there is a slight issue with on-axes coordinate transforms which is what the +.01 is for
      vecPt[3*ind + 0] = (a-epsilon)*PetscCosReal(pu)*PetscSinReal(pv) +.01;
      vecPt[3*ind + 1] = (b-epsilon)*PetscSinReal(pu)*PetscSinReal(pv) + .01;
      vecPt[3*ind + 2] = (c-epsilon)*PetscCosReal(pv) + .01;
      printf("pu: %3.3f pv: %3.3f\n", pu, pv);
    }
  }
  ierr = VecRestoreArray(chargeXYZ, &vecPt);CHKERRQ(ierr);

  //set solution points vector
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*NUM_SOL_PTS, &solXYZ);CHKERRQ(ierr);
  ierr = VecGetArray(solXYZ, &vecPt);CHKERRQ(ierr);
  vecPt[0] = sol1X; vecPt[1] = sol1Y; vecPt[2] = sol1Z;
  vecPt[3] = sol2X; vecPt[4] = sol2Y; vecPt[5] = sol2Z;
  ierr = VecRestoreArray(solXYZ, &vecPt);CHKERRQ(ierr);
  
  //caulculate exact solution
  ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, &exactSolution);CHKERRQ(ierr);
  ierr = CoulombExact(eps1, chargeXYZ, chargeMag, solXYZ, exactSolution);CHKERRQ(ierr);

  //initialize solutions and error vectors
  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, ellSolutions+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, sphSolutions+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, ellError+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, sphError+k);CHKERRQ(ierr);
  }

  //compute solutions for expansions up to order MAX_N
  for(PetscInt k=1; k<=MAX_N; ++k) {
    printf("N = %d\n", k);
    printf("calculating ellipsoidal solution...\n");
    ierr = CalcEllipsoidCoulombPotential(a, b, c, eps1, eps2, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, chargeXYZ, chargeMag, NUM_SOL_PTS, solXYZ, k, ellSolutions[k-1]);CHKERRQ(ierr);
    printf("calculating spherical solution...\n");
    ierr = CalcSphericalCoulombPotential(1.0, eps1, eps2, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, chargeXYZ, chargeMag, NUM_SOL_PTS, solXYZ, k, sphSolutions[k-1]);CHKERRQ(ierr);
  }

  //compute errors
  printf("Computing errors...\n");
  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecCopy(exactSolution, ellError[k]);CHKERRQ(ierr);
    ierr = VecCopy(exactSolution, sphError[k]);CHKERRQ(ierr);
    ierr = VecAXPY(ellError[k], -1.0, ellSolutions[k]);CHKERRQ(ierr);
    ierr = VecAXPY(sphError[k], -1.0, sphSolutions[k]);CHKERRQ(ierr);
    ierr = VecAbs(ellError[k]);CHKERRQ(ierr);
    ierr = VecAbs(sphError[k]);CHKERRQ(ierr);
  }

  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecGetArrayRead(ellError[k], &vecPt);
    printf("Error[%d] = %4.4e\n", k, vecPt[0]);
    ierr = VecRestoreArrayRead(ellError[k], &vecPt);
  }
  printf("Exact: \n");
  ierr = VecView(exactSolution, PETSC_VIEWER_STDOUT_SELF);
  printf("Ell: \n");
  ierr = VecView(ellSolutions[MAX_N-1], PETSC_VIEWER_STDOUT_SELF);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "runTest1"
void runTest1() {
  
  const int nCharges = 1;

  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0, 64);

  //grid dimensions
  int nx = 10;
  int ny = 10;
  int nz = 2;
  int nPoints = nx*ny*nz;
  double xa = -5;
  double xb = 5;
  double ya = -5;
  double yb = 5;
  double za = -3;
  double zb = 3;
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
	printf("rs[3] = %15.15f\n", rs[index].x3);
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


  FILE *fp = fopen("ellArgData.txt", "r");
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
  initEllipsoidalSystem(&e, a, b, c, 64);
  
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

  printf("\n");
  printf("Free energy: %15.15f\n", freeEnergy*.5*cf);
  
  free(xyz);
  free(chargeValues);
  free(solPoints);
  free(solution);
  free(chargePoints);
  PetscFunctionReturn(0);
}

/* DOESN'T DO ANYTHING YET */
#undef __FUNCT__
#define __FUNCT__ "RunArgWorkPrec"
PetscErrorCode RunArgWorkPrec()
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


  FILE *fp = fopen("ellArgData.txt", "r");
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
  initEllipsoidalSystem(&e, a, b, c, 64);
  
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
  initEllipsoidalSystem(&e, a, b, c, 64);
  
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
#define __FUNCT__ "CoulombExactZeroInside"
PetscErrorCode CoulombExactZeroInside(PetscReal a, PetscReal b, PetscReal c, PetscReal eps1, Vec srcXYZ, Vec srcMag, Vec tarXYZ, Vec sol)
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

  ierr = VecZeroEntries(sol);CHKERRQ(ierr);
  
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
    if((tarX*tarX)/(a*a) + (tarY*tarY)/(b*b) + (tarZ*tarZ)/(c*c) <= 1)
      solArray[i] = 0;
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
  ierr = PetscRandomSetInterval(rnd2, .9, 1.0);CHKERRQ(ierr);
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
      xyzArray[3*i+0] = 2.02;
      xyzArray[3*i+1] = -0.4;
      xyzArray[3*i+2] = -0.463;
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
#define __FUNCT__ "SphereLimitExample"
PetscErrorCode SphereLimitExample()
{
  PetscErrorCode ierr;
  const PetscInt NUM_DELTAS = 14;
  const PetscReal DELTAS[14] = {1.3, 1.0, .5, .25, .15, .05, .01, .005, .001, .0005, .0001, .00005, .00001, .000005};
  const PetscInt  Nmax = 10;
  const PetscReal eps1 = 4.0;
  const PetscReal eps2 = 80.0;
  Vec srcXYZ[NUM_DELTAS], srcMag[NUM_DELTAS], solution[NUM_DELTAS];
  PetscScalar *srcXYZArray, *srcMagArray;
  PetscReal a, b, c;
  PetscReal delta;
  PetscReal *freeE;
  PetscInt i;
  PetscLogEvent events[NUM_DELTAS];
  char tText[40] = "Free energy with delta=%8.8f";
  char eText[40];
  EllipsoidalSystem ells[NUM_DELTAS];
  PetscFunctionBegin;

  ierr = PetscMalloc1(sizeof(PetscReal)*NUM_DELTAS, &freeE);

  for(i=0; i<NUM_DELTAS; ++i) {
    /* get size of delta and set a,b,c */
    delta = DELTAS[i];
    a = 1.0 + delta;
    b = 1.0 + (delta/5.0);
    c = 1.0 + (delta/10.0);
    /* init ellipsoidal system */
    ierr = initEllipsoidalSystem(ells+i, a, b, c, 64);CHKERRQ(ierr);
    /* initialize xyz and solution vectors */
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3, srcXYZ+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, 1, srcMag+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, 1, solution+i);CHKERRQ(ierr);
    ierr = VecZeroEntries(srcXYZ[i]);CHKERRQ(ierr);
    ierr = VecGetArray(srcMag[i], &srcMagArray);CHKERRQ(ierr);
    srcMagArray[0] = 1.0; // set only charge to have magnitude 1.0
    ierr = VecRestoreArray(srcMag[i], &srcMagArray);CHKERRQ(ierr);
    /* register events for flop counting */
    sprintf(eText, tText, delta);
    ierr = PetscLogEventRegister((const char*) eText, 0, events+i);CHKERRQ(ierr); 
  }
  /* solutions */
  for(i=0; i<NUM_DELTAS; ++i) {
    printf("calculating solution with delta=%8.8f\n", DELTAS[i]);
    ierr = PetscLogEventBegin(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = CalcEllipsoidFreeEnergy(ells+i, eps1, eps2, 1, srcXYZ[i], srcMag[i], 1e-5, Nmax, solution[i], freeE+i);CHKERRQ(ierr);
    printf("freeE: %15.15f\n", freeE[i]);
    
    ierr = PetscLogEventEnd(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
  }
  FILE *fp = fopen("out/sphConv.txt", "w");
  fprintf(fp, "delta error");
  PetscReal err;
  PetscReal exact = 1./eps2 - 1./eps1;
  for(i=0; i<NUM_DELTAS; ++i) {
    err = PetscAbsReal((exact - freeE[i])/exact);
    fprintf(fp, "%4.4e %4.4e\n", DELTAS[i], err);
    ierr = VecDestroy(srcXYZ+i);CHKERRQ(ierr);
    ierr = VecDestroy(srcMag+i);CHKERRQ(ierr);
    ierr = VecDestroy(solution+i);CHKERRQ(ierr);
  }
  fclose(fp);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ChargeFlopsExample"
PetscErrorCode ChargeFlopsExample(PetscInt Nmax)
{
  PetscErrorCode ierr;
  
  const PetscInt NUM_SOLUTIONS = 10;
  const PetscInt CHARGE_START = 100;
  const PetscInt CHARGE_INC   = 50;
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
  const PetscReal eps2 = 4.0;

  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  ierr = initEllipsoidalSystem(&e, a, b, c, 64);CHKERRQ(ierr);
  
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
  

  PetscReal freeE;
  for(i=0; i < NUM_SOLUTIONS; ++i) {
    printf("calculating solution %d/%d\n", i+1, NUM_SOLUTIONS);
    ierr = PetscLogEventBegin(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = CalcEllipsoidFreeEnergy(&e, eps1, eps2, chargeNums[i], srcXYZ[i], srcMag[i], 1e-5, Nmax, solution[i], &freeE);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
  }
  //ierr = EasyExample(Nmax, nSrc, nx, xl, xr, ny, yl, yr, zConst);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "WorkPrecExample"
PetscErrorCode WorkPrecExample()
{
  PetscErrorCode ierr;
  
  const PetscInt NUM_SOLUTIONS = 45;
  const PetscInt EXACT_NUM   = NUM_SOLUTIONS+10;
  const PetscInt NUM_CHARGES = 5;
  Vec srcXYZ, srcMag;
  Vec potential[NUM_SOLUTIONS];
  Vec potExact;
  Vec solution;
  PetscScalar *solutionArray;
  PetscReal exact;
  PetscInt i;
  PetscLogEvent events[NUM_SOLUTIONS];
  const PetscInt Nmin = 0;
  const PetscInt Nstep = 1;
  PetscInt N;
  char tText[30] = "Free energy with Nmax = %d";
  char eText[30];
  PetscInt tempNum;

  EllipsoidalSystem e;
  PetscFunctionBegin;

  const PetscReal eps1 = 4.0;
  const PetscReal eps2 = 80.0;

  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  ierr = initEllipsoidalSystem(&e, a, b, c, 64);CHKERRQ(ierr);
  
  /* create charge vectors and initialize with random values */
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*NUM_CHARGES, &srcXYZ);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,   NUM_CHARGES, &srcMag);CHKERRQ(ierr);
  ierr = RandomEllipsoidPoints(a, b, c, srcXYZ);CHKERRQ(ierr);
  ierr = RandomMag(srcMag);CHKERRQ(ierr);
  ierr = VecSet(srcMag, 1.0);CHKERRQ(ierr);
  //ierr = VecScale(srcMag, 1./NUM_CHARGES);CHKERRQ(ierr); //scale charge magnitudes so total energy is more normalized
  
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    /* register events for flop counting */
    N = Nmin + Nstep*i;
    sprintf(eText, tText, N);
    ierr = PetscLogEventRegister((const char*) eText, 0, events+i);CHKERRQ(ierr);
    /* create potential vectors */
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_CHARGES, potential+i);CHKERRQ(ierr);
  }
  /* create free energy solution vector */
  ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOLUTIONS, &solution);CHKERRQ(ierr);
  
  ierr = VecGetArray(solution, &solutionArray);CHKERRQ(ierr);
  /* approximate solutions */
  for(i=0; i < NUM_SOLUTIONS; ++i) {
    N = Nmin + Nstep*i;
    printf("calculating solution %d/%d\n", i+1, NUM_SOLUTIONS);
    printf("NNNNN = %d\n", N);
    ierr = PetscLogEventBegin(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = CalcEllipsoidFreeEnergy(&e, eps1, eps2, NUM_CHARGES, srcXYZ, srcMag, 1e-16, N, potential[i], solutionArray+i);CHKERRQ(ierr);
    //ierr = VecView(solution, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = VecView(potential[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(events[i], 0, 0, 0, 0);CHKERRQ(ierr);
  }

  
  /* exact(ish) solution */
  ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_CHARGES, &potExact);CHKERRQ(ierr);
  printf("calculating exact solution\n");
  ierr = CalcEllipsoidFreeEnergy(&e, eps1, eps2, NUM_CHARGES, srcXYZ, srcMag, 1e-16, EXACT_NUM, potExact, &exact);CHKERRQ(ierr);


  FILE *fp = fopen("out/workprec.txt", "w");
  fprintf(fp, "work   precision\n");
  PetscReal errVal;
  PetscEventPerfInfo info;
  //* calculate error for approximate solutions */
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    errVal = PetscAbsReal((solutionArray[i] - exact)/(exact));
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, events[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%4.4e %4.4e\n", info.flops, errVal);
    printf("exact: %15.15f\n", solutionArray[i]);
  }
  fclose(fp);
  ierr = VecRestoreArray(solution, &solutionArray);CHKERRQ(ierr);

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
#define __FUNCT__ "SubtractCoulombInside"
PetscErrorCode SubtractCoulombInside(PetscReal a, PetscReal b, PetscReal c, Vec tarXYZ, Vec solution, Vec solutionEx, Vec ncSolution)
{
  PetscErrorCode ierr;
  const PetscScalar *tarXYZArray, *solutionArray, *solutionExArray;
  PetscScalar *ncSolutionArray;
  PetscInt nTar;
  PetscInt i;
  PetscReal x, y, z;
  PetscFunctionBegin;
  
  ierr = VecGetSize(solution, &nTar);CHKERRQ(ierr);

  ierr = VecGetArrayRead(tarXYZ, &tarXYZArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(solution, &solutionArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(solutionEx, &solutionExArray);CHKERRQ(ierr);
  ierr = VecGetArray(ncSolution, &ncSolutionArray);CHKERRQ(ierr);

  for(i=0; i<nTar; ++i) {
    x = tarXYZArray[3*i+0];
    y = tarXYZArray[3*i+1];
    z = tarXYZArray[3*i+2];

    if( ((x*x)/(a*a)) + ((y*y)/(b*b)) + ((z*z)/(c*c)) <= 1) {
      ncSolutionArray[i] = solutionArray[i];
    }
    else {
      ncSolutionArray[i] = solutionArray[i] - solutionExArray[i];
    }
  }

  ierr = VecRestoreArrayRead(tarXYZ, &tarXYZArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(solution, &solutionArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(solutionEx, &solutionExArray);CHKERRQ(ierr);
  ierr = VecRestoreArray(ncSolution, &ncSolutionArray);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridAnimation"
PetscErrorCode GridAnimation()
{
  const PetscInt NUM_SOLUTIONS = 8;
  PetscErrorCode ierr;
  Vec sourceXYZ, sourceMag, tarXYZ;
  Vec solution[NUM_SOLUTIONS], solutionEx[NUM_SOLUTIONS], ncSolution[NUM_SOLUTIONS];
  Vec errorVec[NUM_SOLUTIONS];
  Vec xOut, yOut, zOut;
  PetscReal xh, yh;
  PetscReal xc, yc;
  PetscReal randVal;
  PetscInt ind;
  PetscFunctionBegin;

  PetscReal eps1 = 2.0;
  PetscReal eps2 = 20.0;
  PetscInt nSrc  = 5;
  PetscInt nx    = 15;
  PetscReal xl   = -5.13;
  PetscReal xr   = 5.13;
  PetscInt ny    = 15;
  PetscReal yl   = -4.02;
  PetscReal yr   =  4.02;
  PetscReal zConst = .1;


  xh = (xr-xl)/ny;
  yh = (yr-yl)/ny;
  
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;

  // initialize source XYZ, mag
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nSrc , &sourceXYZ);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSrc   , &sourceMag);CHKERRQ(ierr);
  // initialize solution vectors
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , solution+i); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , solutionEx+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , ncSolution+i);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nx*ny  , errorVec+i);CHKERRQ(ierr);
    
  }


  ierr = RandomEllipsoidPoints(a, b, c, sourceXYZ);CHKERRQ(ierr); // on z=0
  
  //ierr = RandomMag(sourceMag); CHKERRQ(ierr);
  ierr = VecSet(sourceMag, 1.0);CHKERRQ(ierr);
  
  ierr = InitGrid(nx, xl, xr, ny, yl, yr, zConst, &tarXYZ);CHKERRQ(ierr); //zConst = 0


  //Vec tester;
  //ierr = VecDuplicate(tarXYZ, &tester);
  PetscReal err;
  // calculate the solvation potential
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    printf("wow\n");
    ierr = CoulombExact(eps1, sourceXYZ, sourceMag, tarXYZ, solutionEx[i]);
    //ierr = VecCopy(tarXYZ, tester);
    printf("eps: %15.15f\n", eps1);
    ierr = CalcEllipsoidTester( a , b , c,
				eps1, eps2,
				nSrc, sourceXYZ, sourceMag,
				nx*ny, tarXYZ,
				i, solution[i]); CHKERRQ(ierr);
    ierr = SubtractCoulombInside(a, b, c, tarXYZ, solution[i], solutionEx[i], ncSolution[i]);CHKERRQ(ierr);
    PetscInt zero = 0;
    PetscReal sol0, solEx0;
    ierr = VecGetValues(solution[i], 1, &zero, &sol0);CHKERRQ(ierr);
    ierr = VecGetValues(solutionEx[i], 1, &zero, &solEx0);CHKERRQ(ierr);
    //ierr = VecAXPY(tester, -1.0, tarXYZ);CHKERRQ(ierr);
    //ierr = VecView(solutionEx[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    //ierr = VecView(solution[i], PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = VecCopy(solutionEx[i], errorVec[i]);CHKERRQ(ierr);
    ierr = VecAXPY(errorVec[i], -1.0, solution[i]);
    ierr = VecNorm(errorVec[i], NORM_2, &err);
    printf("\n");
    printf("the error on %d is %15.15f\n", i, err*xh*yh);
  }


  /* subtract coulomb from solution inside */


  
  char fname[20] = "out/sol%d.txt";
  char fnameEx[20] = "out/solEx%d.txt";
  char fnameErr[20] = "out/solErr%d.txt";
  char fnameS[20];
  char fnameExS[20];
  char fnameErrS[20];
  for(PetscInt i=0; i < NUM_SOLUTIONS; ++i) {
    sprintf(fnameS, fname, i);
    sprintf(fnameExS, fnameEx, i);
    sprintf(fnameErrS, fnameErr, i);
    ierr = WriteToFile((const char*) fnameS, ny, ncSolution[i]);CHKERRQ(ierr);
    ierr = WriteToFile((const char*) fnameExS, ny, solutionEx[i]);CHKERRQ(ierr);
    ierr = WriteToFile((const char*) fnameErrS, ny, errorVec[i]);CHKERRQ(ierr);
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
    ierr = VecDestroy(solutionEx+i);CHKERRQ(ierr);
    ierr = VecDestroy(ncSolution+i);CHKERRQ(ierr);
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
#define __FUNCT__ "testLame"
PetscErrorCode testLame()
{
  PetscErrorCode ierr;
  PetscInt n, p;
  PetscInt Nmax;
  PetscReal Fnp;
  PetscReal Enp;
  PetscReal a, b, c;
  PetscReal cval;
  PetscInt signm, signn;
  EllipsoidalSystem e;
  PetscFunctionBegin;

  cval = .3;
  signm = -1;
  signn = 1;
  
  a = 3.0; b = 2.0; c = 1.0;
  initEllipsoidalSystem(&e, a, b, c, 64);
  FILE *fp = fopen("newellout.txt", "w");
  Nmax = 6;
  for(n=0; n < Nmax; ++n) {
    for(p=0; p < 2*n+1; ++p) {
      ierr = calcLame(&e, n, p, cval, signm, signn, &Enp);CHKERRQ(ierr);
      fprintf(fp, "%d %d %15.15f\n", n, p, Enp);
    }
  }
  fclose(fp);

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


  //ierr = EllVsSphConvergence();CHKERRQ(ierr);
  
  //testSphericalCoordinates();
  //runTest1();
  //testLegendre();
  //runTest4();
  //test4sphere();
  //RunArg();
  //ierr = RunArgTester(); CHKERRQ(ierr);
  //PetscErrorCode GridSolution(PetscInt Nmax, PetscInt nSrc, PetscInt nx, PetscReal xl, PetscReal xr, PetscInt ny, PetscReal yl, PetscReal yr, PetscReal zConst)
  /*
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
  */
  //ierr = GridSolution(Nmax, nSrc, nx, xl, xr, ny, yl, yr, zConst, NULL, NULL);CHKERRQ(ierr);
  //ierr = SolutionAnimation(); //<---- total crap
  //runTest1();
  //testLame();


//  ierr = runTest1();CHKERRQ(ierr);
  //ierr = RunArg();
  //ierr = GridAnimation();CHKERRQ(ierr);
  //ierr = ChargeFlopsExample(10);CHKERRQ(ierr);
  ierr = WorkPrecExample();CHKERRQ(ierr);
  //ierr = SphereLimitExample();CHKERRQ(ierr);




  //ierr = numChargesPlot(100, 500, 100);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
