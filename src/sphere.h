#ifndef __sphere
#define __sphere

typedef struct SPoint {
  double x1;
  double x2;
  double x3;
  char type;
} SPoint;

void cartesianToSpherical(SPoint*);
void sphericalToCartesian(SPoint*);

#endif
