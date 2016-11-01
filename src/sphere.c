#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sphere.h"

void cartesianToSpherical(SPoint *point) {
  
  double r, theta, phi;

  r = sqrt((point->x1)*(point->x1) + (point->x2)*(point->x2) + (point->x3)*(point->x3));  
  theta = acos(point->x3/r);
  phi = atan2((point->x2),(point->x1));

  point->x1 = r;
  point->x2 = 0;
  point->x3 = phi;

  if(r > 0) point->x2 = theta;

  point->type = 's';
}

void sphericalToCartesian(SPoint *point) {
  double r = point->x1;
  double theta = point->x2;
  double phi = point->x3;

  double x = r*sin(theta)*cos(phi);
  double y = r*sin(theta)*sin(phi);
  double z = r*cos(theta);

  point->x1 = x;
  point->x2 = y;
  point->x3 = z;
}
