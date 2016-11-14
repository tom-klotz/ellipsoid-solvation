#ifndef __tanhsinh
#define __tanhsinh


PetscErrorCode DExkwk(PetscInt n, mpfr_t h, mpfr_t **xk, mpfr_t **wk);
PetscErrorCode DEQuad(PetscErrorCode (*f)(mpfr_t*,mpfr_t*,void*), mpfr_t a, mpfr_t b, PetscInt prec, PetscInt nPts, PetscReal *integral, void *ctx);
PetscErrorCode inte(mpfr_t *x, mpfr_t *val, void *ctx);


#endif
