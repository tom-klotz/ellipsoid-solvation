#ifndef __ellSolv
#define __ellSolv


//spherical





//ellipsoidal
typedef struct Problem {
  Point *positions;
  double *charges;
  int nCharges;

  double e1; //permitivity inside
  double e2; //permitivity outside
  EllipsoidalSystem *e;
} Problem;

double calcEnp(EllipsoidalSystem*, Point*, int, int);
double calcGnp(EllipsoidalSystem*, Point*, double*, int, int, int);
double calcFnp(EllipsoidalSystem*, Point*, int, int);
void calcBnpAndCnpFromGnp(Problem*, int, int, double, double*, double*);
double calcCoulomb(Problem*, int, Point*);
double calcCoulombEllipsoidalGrid(Problem*, int, Point*, int, double*);

#endif
