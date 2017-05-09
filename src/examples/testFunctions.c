#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "../ellipsoid/ellipsoid.h"
#include "../ellipsoid/ellSolv.h"
#include "../sphere/sphere.h"
#include "../sphere/sphSolv.h"
#include "../constants.h"

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
