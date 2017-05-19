#ifndef __MPFReigs
#define __MPFReigs

PetscErrorCode MatVecMPFR(PetscInt size_d, mpfr_t *Mat, mpfr_t *vec, mpfr_t *result);
PetscErrorCode eigsMPFR(PetscInt size_d, mpfr_t *mat, mpfr_t *result);
PetscErrorCode MatViewMPFR(PetscInt m, PetscInt n, mpfr_t *mat);
PetscErrorCode VecViewMPFR(PetscInt size_d, mpfr_t *vec);
PetscErrorCode OrthogMap(PetscInt size_d, mpfr_t *x, mpfr_t *v1, mpfr_t *lambda,  mpfr_t *result);
PetscErrorCode Vec2NormMPFR(PetscInt size_d, mpfr_t *vec, mpfr_t *norm);
#endif
