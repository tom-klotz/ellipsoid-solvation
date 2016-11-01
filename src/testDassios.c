#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include "ellipsoid.h"

int main() {
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  ellipsoidInitToOrderN(&e, 10);
  int n = 6;
  //double *coefs = (double*) malloc(sizeof(double)*(n/2+1));
  //printf("the value is: %15.15f\n", calcLame2(&e, 3, 1, .5));
  printf("the value is: %15.15f\n", calcLame (&e, 1, 2, 2.1, 1, 1));
  //printf("the value is: %15.15f\n", calcLame2(&e, 3, 1, 2.1));
}
