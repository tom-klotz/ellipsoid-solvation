#ifndef __sphSolv
#define __sphSolv

typedef struct SProblem {
  SPoint *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
  double b;  //radius of molecule
} SProblem;


void calcEnm(SPoint*, double*, int, int, int, double*, double*);
void calcBnmFromEnm(int, int, double, double, SProblem*, double*, double*);
void calcCnmFromEnm(int, int, double, double, SProblem*, double*, double*);
double calcCoulombSpherical(SProblem*, int, SPoint*);
double calcCoulombSphericalGrid(SProblem*, int, SPoint*, int, double*);



#endif
