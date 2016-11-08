#ifndef __ellSolv
#define __ellSolv


typedef struct Problem {
  Point *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
  EllipsoidalSystem *e;
} Problem;

PetscErrorCode CalcEllipsoidSolvationPotential(PetscReal, PetscReal, PetscReal, PetscReal, PetscReal, PetscInt, Vec, Vec, PetscInt, Vec, PetscInt, Vec);
PetscErrorCode CalcSolidInteriorHarmonic(EllipsoidalSystem*, PetscReal, PetscReal, PetscReal, PetscInt, PetscInt, PetscReal*);
PetscErrorCode CalcSolidInteriorHarmonicVec(EllipsoidalSystem*, Vec, PetscInt, PetscInt, Vec);
PetscErrorCode CalcCoulombEllCoefs(EllipsoidalSystem*, PetscInt, Vec, Vec, PetscInt, Vec*);
PetscErrorCode CalcSolidExteriorHarmonic(EllipsoidalSystem*, PetscReal, PetscReal, PetscReal, PetscInt, PetscInt, PetscReal*);
PetscErrorCode CalcSolidExteriorHarmonicVec(EllipsoidalSystem*, Vec, PetscInt, PetscInt, Vec);
PetscErrorCode CalcReactAndExtCoefsFromCoulomb(EllipsoidalSystem*, PetscReal, PetscReal, PetscInt, Vec, Vec, Vec);
PetscErrorCode CalcEllipsoidTester(PetscReal a, PetscReal b, PetscReal c, PetscReal eps1, PetscReal eps2, PetscInt nSrc, Vec srcXYZ, Vec srcMag, PetscInt nTar, Vec tarXYZ, PetscInt Nmax, Vec tarSol);
PetscErrorCode CalcEllipsoidFreeEnergy(EllipsoidalSystem *e, PetscReal eps1, PetscReal eps2, PetscInt nSrc, Vec srcXYZ, Vec srcMag, PetscReal tol, PetscInt Nmax, Vec tarSol);

double calcEnp(EllipsoidalSystem*, Point*, int, int);
double calcGnp(EllipsoidalSystem*, Point*, double*, int, int, int);
double calcFnp(EllipsoidalSystem*, Point*, int, int);
void calcBnpAndCnpFromGnp(Problem*, int, int, double, double*, double*);
double calcCoulomb(Problem*, int, Point*);
double calcCoulombEllipsoidalGrid(Problem*, int, Point*, int, double*);

#endif
