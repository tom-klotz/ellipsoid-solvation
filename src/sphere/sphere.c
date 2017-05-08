#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sphere.h"
#include <petsc.h>
#undef __FUNCT__
#define __FUNCT__ "CartesianToSpherical"
PetscErrorCode CartesianToSpherical(Vec xyz, Vec sph)
{
  PetscErrorCode ierr;
  PetscScalar *xyzPt, *sphPt;
  PetscInt size;
  PetscReal x, y, z, r, theta, phi;
  PetscFunctionBegin;
  //get size of xyz vec
  ierr = VecGetSize(xyz, &size);CHKERRQ(ierr);
  
  ierr = VecGetArrayRead(xyz, &xyzPt);CHKERRQ(ierr);
  ierr = VecGetArray(sph, &sphPt);CHKERRQ(ierr);
  for(PetscInt k=0; k<size/3; ++k) {
    //get xyz points
    x = xyzPt[3*k+0];
    y = xyzPt[3*k+1];
    z = xyzPt[3*k+2];
    //compute transformation
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
    //set values in sph vec
    sphPt[3*k+0] = r;
    sphPt[3*k+1] = 0;
    sphPt[3*k+2] = phi;
    if(r>0) sphPt[3*k+1] = theta;
    
  }
  ierr = VecRestoreArray(sph, &sphPt);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(xyz, &xyzPt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
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
