#include <math.h>
#include <stdio.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include <petsc.h>
#include <time.h>
#include "ellipsoid.h"


extern void dgeev_(char *jobvl, char *jobvr, int *N, double *A, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

extern void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *B, int *ldb, int *info);

extern void dgetrf_(int *M, int *N, double *A, int *lda, int *ipiv, int *info);

static double max(double a, double b) {
  if( a > b)
    return a;
  return b;
}

static double min(double a, double b) {
  if( a > b )
    return b;
  return a;
}

void matTranspose(double *A, int n) {
  double *cpy = (double*) malloc(sizeof(double)*n*n);
  memcpy(cpy, A, sizeof(double)*n*n);
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      A[i*n+j] = cpy[j*n+i];
    }
  }
  free(cpy);
}


#undef __FUNCT__
#define __FUNCT__ "initEllipsoidalSystem"
PetscErrorCode initEllipsoidalSystem(struct EllipsoidalSystem *s, double a, double b, double c)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  //set pointers to NULL and maxN to 0
  s->Dconsts = NULL;
  s->Rconsts = NULL;
  s->normConstants = NULL;
  s->DmaxN = 0;
  s->RmaxN = 0;

  //default init Romain constants to order 40
  int N = 80;
  
  s->a = a;
  s->b = b;
  s->c = c;

  s->precision = 32;

  mpfr_t temp; 
  mpfr_init(temp);
  
  s->h2 = a*a-b*b;
  s->h = sqrt(s->h2);
  s->k2 = a*a-c*c;
  s->k = sqrt(s->k2); flopCount += 8;

  mpfr_inits2(4*s->precision, s->hp_h2, s->hp_h, s->hp_k2, s->hp_k, NULL);
  //calculates high precision h2
  mpfr_set_d(s->hp_h2, b, MPFR_RNDN);
  mpfr_mul_d(s->hp_h2, s->hp_h2, b, MPFR_RNDN);
  mpfr_set_d(temp, a, MPFR_RNDN);
  mpfr_mul_d(temp, temp, a, MPFR_RNDN);
  mpfr_sub(s->hp_h2, temp, s->hp_h2, MPFR_RNDN); flopCount += 3;

  //calculates high precision h
  mpfr_sqrt(s->hp_h, s->hp_h2, MPFR_RNDN); flopCount += 1;

  //calculates high precision k2
  mpfr_set_d(s->hp_k2, c, MPFR_RNDN);
  mpfr_mul_d(s->hp_k2, s->hp_k2, c, MPFR_RNDN);
  mpfr_sub(s->hp_k2, temp, s->hp_k2, MPFR_RNDN); flopCount += 2;
  
  //calculates high precision k
  mpfr_sqrt(s->hp_k, s->hp_k2, MPFR_RNDN); flopCount += 1;


  //initializes all the variables for integrate.
  mpfr_inits2(4*s->precision, s->alpha, s->beta, s->h_step, s->sum, s->osum, s->psum, s->yk, s->wk, s->lx, s->rx, s->tmp, s->maxTerm, s->curTerm, s->pi2, s->kh, s->msinh, s->mcosh, s->lval, s->rval, NULL);

  //initialize other temp variables
  mpfr_inits2(4*s->precision, s->temp1, s->temp2, s->temp3, s->tempa, s->tempb, s->endpt, NULL);

  //constants
  mpfr_inits2(4*s->precision, s->mpfrzero, s->mpfrone, NULL);
  mpfr_set_d(s->mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(s->mpfrone , 1.0, MPFR_RNDN);


  ierr = initRomainConstsToOrderN(s, N);CHKERRQ(ierr);

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void ellipsoidInitToOrderN(EllipsoidalSystem *s, int N)
{
  s->Dconsts = (double***) malloc(sizeof(double**)*(N+1));
  s->normConstants = (double**) malloc(sizeof(double*)*(N+1));
  s->tVals  = (char**) malloc(sizeof(char*)*(N+1));
  s->tpVals = (int**)  malloc(sizeof(int* )*(N+1));
  for(int n=0; n<N+1; ++n) {
    s->Dconsts[n] = (double**) malloc(sizeof(double*)*(2*n+1));
    s->normConstants[n] = (double*) malloc(sizeof(double)*(2*n+1));
    s->tVals[n]  = (char*) malloc(sizeof(char)*(2*n+1));
    s->tpVals[n] = (int *) malloc(sizeof(int) *(2*n+1));
  }
  
  //loop over n
  for(int n=0; n<=N; ++n) {
    getCoefsK(s, n, s->Dconsts[n]);
    if(n != 0) {
      getCoefsL(s, n, s->Dconsts[n] + (n/2)+1);
      getCoefsM(s, n, s->Dconsts[n] + (n/2)+1 + (n+1)/2);
      if(n != 1) {
	getCoefsN(s, n, s->Dconsts[n] + (n/2)+1 + (n+1)/2 + (n+1)/2);
      }
    }
    
    
    //loop over p
    for(int p=0; p<2*n+1; ++p) {
      s->tVals [n][p] = getLameTypeT (n, p);
      s->tpVals[n][p] = getLameTypeTp(n, p);
    }
  }
  
}

void printMatrix(double *M, int matSize) {
  for(int n=0; n<matSize; ++n) {
    for(int m=0; m<matSize; ++m) {
      printf("%4.4f ", M[n*matSize+m]);
    }
    printf("\n");
  }
}  


void eigs(double *M, int matSize, double *real, double *imag) {
  //transpose M for lapack
  matTranspose(M, matSize);


  int lwork = 16*matSize;
  int info;
  double *work, *wr, *wi;
  //vr     = (double*) malloc(sizeof(double)*matSize*matSize);
  wr     = (double*) malloc(sizeof(double)*matSize);
  wi     = (double*) malloc(sizeof(double)*matSize);
  work   = (double*) malloc(sizeof(double)*lwork);


  //compute eigenvalues with lapack
  char complefteig, comprighteig;
  complefteig = 'N';  //don't compute left eigenvectors
  comprighteig = 'N'; //don't compute right eigenvectors

  
  dgeev_( &complefteig, &comprighteig, &matSize, M, &matSize, wr, wi, NULL, &matSize, NULL, &matSize, work, &lwork, &info );
  
  matTranspose(M, matSize);
  
  for(int i=0; i<matSize; ++i)
    real[i] = wr[i];
  if(imag != NULL) {
    for(int i=0; i<matSize; ++i)
      imag[i] = wi[i];
  }
}



#undef __FUNCT__
#define __FUNCT__ "getCoefsK"
//////////////
//FROM DASSIOS
//////////////
PetscErrorCode getCoefsK(EllipsoidalSystem *s, int n, double **coefs) {
  PetscErrorCode ierr;
  PetscInt flopCount;
  
  PetscFunctionBegin;
  //r=n/2 if even, (n-1)/2 if odd
  int r = n/2;
  //matrix size for class K is r+1
  int matSize = r+1;
  //constants
  double alpha = s->h2 + s->k2;
  double beta  = s->h2 * s->k2;

  //initialize matrix and tridiagonals
  double *M = (double*) malloc(sizeof(double)*matSize*matSize);
  double *d = (double*) malloc(sizeof(double)*matSize);
  double *f = (double*) malloc(sizeof(double)*(matSize-1));
  double *g = (double*) malloc(sizeof(double)*(matSize-1));

  flopCount = 2;
  //construct the diagonal (size matSize)
  for(int i=1; i<=matSize; ++i)
    d[i-1] = (n - 2*i + 2)*(n - 2*i + 2);
  flopCount += 7*matSize;
  //construct above and below diagonal (size matSize-1)
  for(int i=1; i<matSize; ++i) {
    g[i-1] = (2./alpha)*i*(2*n - 2*i + 1);
    f[i-1] = -(beta/alpha)*(n - 2*i + 1)*(n - 2*i + 2);
  }
  flopCount += 17*matSize;
  
  //fill diagonal with d
  for(int k=0; k < matSize; ++k)
    M[k*matSize + k]   = d[k];
  //fill above diagonal with g, below with f
  for(int k=0; k < matSize-1; ++k) {
    M[(k)*matSize + k+1] = g[k];
    M[(k+1)*matSize + k] = f[k];
  }

  double *vals = (double*) malloc(sizeof(double)*matSize);
  eigs(M, matSize, vals, NULL);
  
  //for(int i=0; i<matSize; ++i)
  //printf("Kn=%d vals[%d] = %15.15f\n", n, i, vals[i]);

  double pj;
  //loop over each root
  for(int p=0; p<matSize; ++p) {
    coefs[p] = (double*) malloc(sizeof(double)*(matSize));
    pj = vals[p];
    coefs[p][0] = 1.0;
    if(matSize != 1)
      coefs[p][1] = alpha*(pj - n*n)/(2*(2*n - 1));
    for(int k=2; k<matSize; ++k) {
      coefs[p][k] = (alpha*(pj - (n - 2*k + 2)*(n - 2*k + 2)) * coefs[p][k-1]
		     + beta*(n - 2*k + 4)*(n - 2*k + 3) * coefs[p][k-2])
	/(2*k*(2*n - 2*k + 1));
    }
  }
  flopCount += 7*matSize + matSize*PetscMax(0,(matSize-2))*27;
  
  //for(int i=0; i<matSize; ++i)
  //printf("%5.5f\n", vals[i]);
  ierr = PetscLogFlops(flopCount);
  PetscFunctionReturn(0);
}

void getCoefsL(EllipsoidalSystem *s, int n, double **coefs) {
  //r=n/2 if even, (n-1)/2 if odd
  int r = (n+1)/2;
  //matrix size for class K is r+1
  int matSize = r;
  //constants
  double alpha = s->h2 + s->k2;
  double beta  = s->h2 * s->k2;

  //initialize matrix and tridiagonals
  double *M = (double*) malloc(sizeof(double)*matSize*matSize);
  double *d = (double*) malloc(sizeof(double)*matSize);
  double *f = (double*) malloc(sizeof(double)*(matSize-1));
  double *g = (double*) malloc(sizeof(double)*(matSize-1));

  //construct the diagonal (size matSize)
  for(int i=1; i<=matSize; ++i)
    d[i-1] = (n - 2*i + 1)*(n - 2*i + 1) + (s->k2/alpha)*(2*n - 4*i + 3);
  //construct above and below diagonal (size matSize-1)
  for(int i=1; i<matSize; ++i) {
    g[i-1] = (2./alpha)*i*(2*n - 2*i + 1);
    f[i-1] = -(beta/alpha)*(n - 2*i + 1)*(n - 2*i);
  }
  
  //fill diagonal with d
  for(int k=0; k < matSize; ++k)
    M[k*matSize + k]   = d[k];
  //fill above diagonal with g, below with f
  for(int k=0; k < matSize-1; ++k) {
    M[(k)*matSize + k+1] = g[k];
    M[(k+1)*matSize + k] = f[k];
  }


  double *vals = (double*) malloc(sizeof(double)*matSize);
  eigs(M, matSize, vals, NULL);
  
  for(int i=0; i<matSize; ++i)
    printf("n=%d vals[%d] = %15.15f\n", n, i, vals[i]);
  //printf("%15.15f\n", (1/alpha)*(2*alpha + 3*s->k2 - 2*sqrt((alpha+s->k2)*(alpha+s->k2) - 5*beta)));
  
  double pj;
  //loop over each root
  for(int p=0; p<matSize; ++p) {
    coefs[p] = (double*) malloc(sizeof(double)*(matSize));
    pj = vals[p];
    coefs[p][0] = 1.0;
    if(matSize != 1)
      coefs[p][1] = (alpha*(pj - (n-1)*(n-1)) - (2*n - 1)*(2*n - 1)*s->k2)/(2*(2*n-1));
    for(int k=2; k<matSize; ++k) {
      coefs[p][k] = ((alpha*(pj - (n - 2*k + 3)*(n - 2*k + 3)) - (2*n - 4*k + 7)*s->k2) * coefs[p][k-1]
		     + beta*(n - 2*k + 5)*(n - 2*k + 4)*coefs[p][k-2])/(2*(k - 1)*(2*n - 2*k + 3));
    }
  } 
    
}


void getCoefsM(EllipsoidalSystem *s, int n, double **coefs) {
  int wow = 1;
  wow = wow/2;
}
void getCoefsN(EllipsoidalSystem *s, int n, double **coefs) {
  int wow = 1;
  wow = wow/2;
}


#undef __FUNCT__
#define __FUNCT__ "calcLame2"
PetscErrorCode calcLame2(EllipsoidalSystem *s, int n, int p, double l, double *sol)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  char t = getLameTypeT(n, p);
  *sol = 0;
  if (t == 'K') {
    int r = n/2;
    for(int k=0; k <= r; ++k) {
      *sol += s->Dconsts[n][p][k] * pow(l,n - 2*k); flopCount += 2;
    }
  }
  if (t == 'L') {
    int r = (n+1)/2;
    for(int k=0; k <= r-1; ++k) {
      *sol += s->Dconsts[n][p][k] * pow(l, n - 1 - 2*k); flopCount += 2;
    }
    *sol *= sqrt(fabs(l*l - s->h2)); flopCount += 3;
  }


  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ellipsoidToCartesian"
PetscErrorCode ellipsoidToCartesian(EllipsoidalSystem *s, Point *p)
{
  PetscErrorCode ierr;


  PetscFunctionBegin;
  if(p->type != 'c') {
    double h2 = s->h2;
    double k2 = s->k2;
  
    double l = p->x1;
    double m = p->x2;
    double n = p->x3;

    //Romain (7)
    p->x1 = sqrt((l * l * m * m * n * n)/(h2 * k2));
    p->x2 = sqrt(((l * l - h2) * (m * m - h2) * (h2-n * n)) / (h2 * (k2 - h2)));
    p->x3 = sqrt(((l * l - k2) * (k2 - m * m) * (k2 - n * n)) / (k2 * (k2 - h2)));
    

    if(l * m * n < 0)
      p->x1 = -p->x1;
    if(l * n < 0)
      p->x2 = -p->x2;
    if(l * m < 0)
      p->x3 = -p->x3;
    
    p->type = 'c';
  }

  ierr = PetscLogFlops(36);
  PetscFunctionReturn(0);
} 


void cartesianToEllipsoidal2(struct EllipsoidalSystem *e, struct Point *p) {

  if(p->type != 'e') {
    double h2 = e->h2;
    double k2 = e->k2;
    double x  = p->x1;
    double y  = p->x2;
    double z  = p->x3;
    double a  = -1.0*(x*x + y*y + z*z + h2 + k2);
    double b  = (h2 + k2)*x*x + k2*y*y + h2*z*z + h2*k2;
    double c  = -h2*k2*x*x;

    double Areal  = -(c/2.0) + (a*b/6.0) - (a*a*a/27.0);
    double Aimag  = sqrt(-1.0*((a*a*a*c + b*b*b - a*a*b*b)/27.0 + (9.0*c*c - 6.0*a*b*c + a*a*b*b)/36.0));
    
    double r = sqrt(Areal*Areal + Aimag*Aimag);
    double Atheta = atan2(Aimag, Areal);

    //Not sure what I was doing with Broot1, Broot2. Commented out
    double Aroot1, Aroot2;//, Broot1, Broot2;
    double ArootReal, ArootImag; //, BrootReal, BrootImag;
    int Ayes = 0;

    for(int k=0; k<3; ++k) {
      Aroot1 = cbrt(r)*cos((Atheta + 2.0*k*PETSC_PI)/3.0);
      Aroot2 = cbrt(r)*sin((Atheta + 2.0*k*PETSC_PI)/3.0);
      //what is this jibber jabber
      //Broot1 = cbrt(r)*cos((Btheta + 2.0*k*PETSC_PI)/3.0);
      //Broot2 = cbrt(r)*sin((Btheta + 2.0*k*PETSC_PI)/3.0);
      if(Aroot1 >= 0 && Aroot2 >= 0) {
	ArootReal = Aroot1;
	ArootImag = Aroot2;
	Ayes = 1;
      }
      //Not sure what this is for
      //
      //if(Broot1 >= 0 && Broot2 >= 0) {
      //BrootReal = Broot1;
      //BrootImag = Broot2;
      //}
      //printf("A,Ai,B,Bi: %15.15f, %15.15f, %15.15f, %15.15f\n", Aroot1, Aroot2, Broot1, Broot2);
    }
    if(Ayes == 0)
      printf("\n\n\n\n COULD NOT FIND PROPER ROOT FOR A");
  
    double kap1 = -a/3.0 + 2.0*ArootReal;
    double kap2 = -a/3.0 - ArootReal - sqrt(3.0)*ArootImag;
    double kap3 = -a/3.0 - ArootReal + sqrt(3.0)*ArootImag;
    
    double newX1 = sqrt(kap1);
    double newX2 = sqrt(kap2);
    double newX3 = sqrt(kap3);
    
    //p->x1 = sqrt(kap1);
    //p->x2 = sqrt(kap3);
    //p->x3 = sqrt(kap2);
    
    if(x*y*z < 0)
      newX1 = - newX1;
    if(x*y < 0)
      newX2 = -newX2;
    if(x*z < 0)
      newX3 = -newX3;

    p->x1 = newX1;
    p->x2 = newX2;
    p->x3 = newX3;
  }
  p->type = 'e';
}

void cartesianToEllipsoidal(struct EllipsoidalSystem *s, struct Point *p)
{
  if(p->type != 'e') {
    double h2 = s->h2;
    double k2 = s->k2;
    double x = p->x1;
    double y = p->x2;
    double z = p->x3;
    //Use Romain (6) for approximate values
    double a1 = -(x*x + y*y + z*z + h2 + k2);
    double a2 = x*x*(h2 + k2) + y*y*k2 + z*z*h2 + h2*k2;
    double a3 = -x*x*h2*k2;
    double Q = (a1*a1 - 3.*a2)/9.;
    double R = (9.*a1*a2 - 27.*a3 - 2.*a1*a1*a1)/54.;
    double cth = R/sqrt(Q*Q*Q);
    double th = acos(cth);
    
    /*
      #    Wikipedia http://en.wikipedia.org/wiki/List_of_trigonometric_identities                                             
      #     cos x = 4 cos^3 x/3 - 3 cos x/3                                                                                   
      #     4 x^3 - 3x - v = 0  ==> {a = 4, b = 0, c = -3, d = -v}                                                            
      #     Discriminant = 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2                                                         
      #                  = 0 - 0 + 0 + 432 - 432*v^2 = 432*(1 - v^2) > 0 so 3 real roots                                      
      #     r1 = -1/3a (cbrt((27a^2d + sqrt(729a^4d^2 + 108a^3c^3))/2) + cbrt((27a^2d - sqrt(729a^4d^2 + 108a^3c^3))/2))      
      #        = -1/12 (cbrt((432d + 432 sqrt(d^2 - 1))/2) + cbrt((432d - 432 sqrt(d^2 - 1))/2))                              
      #        = -cbrt(216)/12 (cbrt((d + sqrt(v^2 - 1))) + cbrt((d - sqrt(v^2 - 1))))                                        
      #        = cbrt(216)/(12) (cbrt(cth - i sth) + cbrt(cth + i sth))      
    */
    double lambda[3] = { cos(th/3.0),
			 cos((th+4.0*PETSC_PI)/3.0),
			 cos((th+2.0*PETSC_PI)/3.0) };
    
    double val3 = 2*sqrt(Q)*lambda[2] - a1/3.0;
    if(val3 < 0 && fabs(val3) < 1e-10)
      val3 = 0;
    else if(val3 < 0) {
      printf("\n\nTHERE WAS AN ISSUE WITH THE TRANSFORM and we got %8.8e\n\n", 2*sqrt(Q)*lambda[2] -a1/3.0);
      val3 = 0;
    }
    //for(int i=0; i<3; ++i) printf("lambda[%d] = %15.15f\n", i, lambda[i]);
    p->x1 = sqrt(2*sqrt(Q)*lambda[0] - a1/3);
    p->x2 = sqrt(2*sqrt(Q)*lambda[1] - a1/3);
    p->x3 = sqrt(val3);

    //printf("p->x1: %15.15f\n", p->x1);
    //printf("p->x2: %15.15f\n", p->x2);
    //printf("p->x3: %15.15f\n", p->x3);
    /*
      # Improve the estimate                                                                                                  
      #   Get $\min(\lambda^2_i, |\lambda^2 - h^2|, |\lambda^2_i - k^2|)$                                                     
      #   Rewrite ellipsoid cubic so that smallest value is root                                                              
      #   Solve the equation using Newton                                                                                     
      #   Transform back
    */
    if(x*y*z < 0)
      p->x1 = -p->x1;
    if(x*y < 0)
      p->x2 = -p->x2;
    if(x*z < 0)
      p->x3 = -p->x3;
    
    p->type = 'e';
  }
}


#undef __FUNCT__
#define __FUNCT__ "CartesianToEllipsoidalVec"
PetscErrorCode CartesianToEllipsoidalVec(EllipsoidalSystem *e, Vec xyzP, Vec ellP)
{
  PetscErrorCode ierr;
  PetscInt  nPts;
  PetscReal h2, k2;
  PetscReal x, y, z;
  PetscReal a1, a2, a3;
  PetscReal Q, R, cth, th;
  PetscReal val3;
  PetscReal lambda, mu, nu;
  const PetscScalar *xyzPArray;
  PetscScalar *ellPArray;
  PetscFunctionBegin;

  h2 = e->h2;
  k2 = e->k2;
  ierr = VecGetSize(xyzP, &nPts);CHKERRQ(ierr);
  nPts = nPts/3;
  
  ierr = VecGetArray    (ellP, &ellPArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(xyzP, &xyzPArray);CHKERRQ(ierr);
  for(PetscInt i=0; i < nPts; ++i) {

    x = xyzPArray[3*i+0];
    y = xyzPArray[3*i+1];
    z = xyzPArray[3*i+2];
    
    //Use Romain (6) for approximate values
    a1 = -(x*x + y*y + z*z + h2 + k2);
    a2 = x*x*(h2 + k2) + y*y*k2 + z*z*h2 + h2*k2;
    a3 = -x*x*h2*k2;
    Q = (a1*a1 - 3.*a2)/9.;
    R = (9.*a1*a2 - 27.*a3 - 2.*a1*a1*a1)/54.;
    cth = R/sqrt(Q*Q*Q);
    th = acos(cth);
    PetscReal lam[3] = { cos(th/3.0),
			    cos((th+4.0*PETSC_PI)/3.0),
			    cos((th+2.0*PETSC_PI)/3.0) };
    
    val3 = 2*sqrt(Q)*lam[2] - a1/3.0;
    if(val3 < 0 && fabs(val3) < 1e-10) {
      printf("WOOOOW\n\n");
      val3 = 0;
    }
    else if(val3 < 0) {
      printf("\n\nTHERE WAS AN ISSUE WITH THE TRANSFORM and we got %8.8e\n\n", 2*sqrt(Q)*lam[2] -a1/3.0);
      val3 = 0;
    }
    
    lambda = sqrt(2*sqrt(Q)*lam[0] - a1/3);
    mu     = sqrt(2*sqrt(Q)*lam[1] - a1/3);
    nu     = sqrt(val3);
    
    if(x*y*z < 0)
      lambda = -lambda;
    if(x*y < 0)
      mu = -mu;
    if(x*z < 0)
      nu = -nu;

    ellPArray[3*i+0] = lambda;
    ellPArray[3*i+1] = mu;
    ellPArray[3*i+2] = nu;
    
  }

  ierr = VecRestoreArrayRead(xyzP, &xyzPArray);CHKERRQ(ierr);
  ierr = VecRestoreArray    (ellP, &ellPArray);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


char getLameTypeT(int n, int p)
{
  int r = n/2;
  if(p < r+1)
    return 'K';
  else if(p < (n-r) + (r+1))
    return 'L';
  else if(p < (n-r) + (n-r) + (r+1))
    return 'M';
  else if(p < 2*n+1)
    return 'N';
  else
    printf("getLameType passed invalid n,p = (%d,%d)\n", n, p);
  return 'E';
}

int getLameTypeTp(int n, int p)
{
  int r = n/2;
  if(p < r+1)
    return p;
  else if(p < (n-r) + (r+1))
    return p - (r+1);
  else if(p < (n-r) + (n-r) + (r+1))
    return p - (n-r) - (r+1);
  else if(p < 2*n+1)
    return p - (n-r) - (n-r) - (r+1);
  else
    printf("getLameTypeTP passed invalid n,p = (%d,%d)\n", n, p);
  return 0;
}

void getLameCoefficientMatrix2(struct EllipsoidalSystem *s, char t, int n, int *mat_size, double *out) {

   //For a given Lame type and ordr n, return the tridiagonal matrix
  //the size of the generated coefficient matrix is stored to mat_size
  //The eigenvectors of this matrix give the coefficients for the polynomaial P in the Lame function definition
  //The expressions come from Romain Annex 3
  //OPT: We could memorize this
  double alpha = s->h2;
  double beta = s->k2 - s->h2;
  double gamma = alpha - beta;
  int r = n/2;
  int size_d;
  double *d, *g, *f;
  if(t == 'K') {
    size_d = r+1;
  }
  else if(t == 'L' || t == 'M') {
    size_d = n-r;
  }
  else if(t == 'N') {
    size_d = r;
  }
  d = (double*) malloc(sizeof(double)*size_d);
  g = (double*) malloc(sizeof(double)*(size_d-1));
  f = (double*) malloc(sizeof(double)*(size_d-1));
  
  if (t == 'K') {
    //g missing last item
    for(int k = 0; k < r; ++k)
      g[k] = -(2*k + 2)*(2*k + 1)*beta;
    if(n%2) { //n is odd
      for(int k=0; k < r+1; ++k)
	d[k] = ((2*r + 1)*(2*r + 2) - 4*k*k)*alpha + (2*k + 1)*(2*k + 1)*beta;
      for(int k=1; k < r+1; ++k)
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) + 1); 
    }
    else { //n is even
      for(int k=0; k < r+1; ++k)
	d[k] = 2*r*(2*r + 1)*alpha - 4*k*k*gamma;
      for(int k=1; k < r+1; ++k)
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) - 1);
    }
  }
  else if (t == 'L') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 3)*beta;
      //printf("g[%d] = %15.15f\n", k, g[k]);
      //printf("n, r, %d, %d\n", n, r);
      //printf("beta: %15.15f\n", beta);
    }
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k)
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 1)*(2*k + 1)*gamma;
      for(int k=1; k < n-r; ++k)
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1);
    }
    else { //n is even
      for(int k=0; k < n-r; ++k)
	d[k] = (2*r*(2*r + 1) - (2*k + 1)*(2*k + 1))*alpha + (2*k + 2)*(2*k + 2)*beta;
      for(int k=1; k < n-r; ++k)
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1);
    }
  }
  else if (t == 'M') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k)
      g[k] = -(2*k + 2)*(2*k + 1)*beta;
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k)
	d[k] = ((2*r + 1)*(2*r + 2) - (2*k + 1)*(2*k + 1))*alpha + 4*k*k*beta;
      for(int k=1; k < n-r; ++k)
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1);
    }
    else { //n is even
      for(int k=0; k < n-r; ++k)
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 1)*(2*k + 1)*gamma;
      for(int k=1; k < n-r; ++k)
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1);
    }
  }
  else if (t == 'N') {
    //g missing last item
    for(int k = 0; k < r-1; ++k)
      g[k] = -(2*k + 2)*(2*k + 3)*beta;
    if(n%2) { //n is odd
      for(int k=0; k < r; ++k)
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 2)*(2*k + 2)*gamma;
      for(int k=1; k < r; ++k)
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k +3);
    }
    else { //n is even
      for(int k=0; k < r; ++k)
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 2)*(2*k + 2)*alpha + (2*k + 1)*(2*k + 1)*beta;
      for(int k=1; k < r; ++k)
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1);
    }
  }


  double *M = (double*) calloc(sizeof(double), size_d*size_d);
  //fill diagonal with d
  for(int k=0; k < size_d; ++k)
    M[k*size_d + k]   = d[k];
  //fill above diagonal with g
  for(int k=0; k < size_d-1; ++k)
    M[(k)*size_d + k+1] = g[k];
  //fill below diagonal with f
  for(int k=0; k < size_d-1; ++k)
    M[(k+1)*size_d + k] = f[k];
    
  
  //transpose M for lapack
  matTranspose(M, size_d);
  

  int lwork = 16*size_d;
  int info;
  double *vr, *work, *wr, *wi;
  vr     = (double*) malloc(sizeof(double)*size_d*size_d);
  wr     = (double*) malloc(sizeof(double)*size_d);
  wi     = (double*) malloc(sizeof(double)*size_d);
  work   = (double*) malloc(sizeof(double)*lwork);

  //compute eigenvalues with lapack
  char complefteig, comprighteig;
  complefteig = 'N'; //don't compute left eigenvectors
  comprighteig = 'V'; //compute right eigenvectors
  

  dgeev_( &complefteig, &comprighteig, &size_d, M, &size_d, wr, wi, NULL, &size_d, vr, &size_d, work, &lwork, &info );
  

  //matTranspose(vr, size_d);  
  //printf("hot\n");
  //printf("size_d: %d\n", size_d);
  free(wr);
  free(wi);
  free(d);
  free(f);
  free(g);
  free(work);
  free(M);
  //printf("amazing\n");
  //printf("diggity\n");
  *mat_size = size_d;

  //return vr;
}

/*
#undef __FUNCT__
#define __FUNCT__ "MatVecMult"
PetscErrorCode MatVecMult(PetscInt m, PetscInt n, double **mat, double *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for(PetscInt i=0; i<m; ++i) {
    for(PetscInt j=0; j<n; ++j) {
      *mat[i*n+j] = 

    }
  }

  PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "getLameCoefficientMatrixSymmetric"
PetscErrorCode getLameCoefficientMatrixSymmetric(struct EllipsoidalSystem *s, char t, int n, int* const mat_size, double **mat)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  //For a given Lame type and ordr n, return the tridiagonal matrix
  //the size of the generated coefficient matrix is stored to mat_size
  //The eigenvectors of this matrix give the coefficients for the polynomaial P in the Lame function definition
  //The expressions come from Romain Annex 3
  //OPT: We could memorize this
  double alpha = s->h2;
  double beta = s->k2 - s->h2;
  double gamma = alpha - beta; flopCount += 2;
  int r = n/2;
  int size_d;
  double *d, *g, *f, *sigma;
  if(t == 'K') {
    size_d = r+1;
  }
  else if(t == 'L' || t == 'M') {
    size_d = n-r;
  }
  else if(t == 'N') {
    size_d = r;
  }
  d = (double*) malloc(sizeof(double)*size_d);
  g = (double*) malloc(sizeof(double)*(size_d-1));
  f = (double*) malloc(sizeof(double)*(size_d-1));
  sigma = (double*) malloc(sizeof(double)*size_d);
    
  if (t == 'K') {
    //g missing last item
    for(int k = 0; k < r; ++k) {
      g[k] = -(2*k + 2)*(2*k + 1)*beta; flopCount++;
    }
    if(n%2) { //n is odd
      for(int k=0; k < r+1; ++k) {
	d[k] = ((2*r + 1)*(2*r + 2) - 4*k*k)*alpha + (2*k + 1)*(2*k + 1)*beta; flopCount += 6;
      }
      for(int k=1; k < r+1; ++k) {
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) + 1); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < r+1; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - 4*k*k*gamma; flopCount += 5;
      }
      for(int k=1; k < r+1; ++k) {
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) - 1); flopCount += 2;
      }
    }
  }
  else if (t == 'L') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 3)*beta; flopCount += 2;
      //printf("g[%d] = %15.15f\n", k, g[k]);
      //printf("n, r, %d, %d\n", n, r);
      //printf("beta: %15.15f\n", beta);
    }
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k) {
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 1)*(2*k + 1)*gamma; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1); flopCount += 3;
      }
    }
    else { //n is even
      for(int k=0; k < n-r; ++k) {
	d[k] = (2*r*(2*r + 1) - (2*k + 1)*(2*k + 1))*alpha + (2*k + 2)*(2*k + 2)*beta; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
  }
  else if (t == 'M') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 1)*beta; flopCount += 2;
    }
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k) {
	d[k] = ((2*r + 1)*(2*r + 2) - (2*k + 1)*(2*k + 1))*alpha + 4*k*k*beta; flopCount += 5;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < n-r; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 1)*(2*k + 1)*gamma; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
  }
  else if (t == 'N') {
    //g missing last item
    for(int k = 0; k < r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 3)*beta; flopCount += 2;
    }
    if(n%2) { //n is odd
      for(int k=0; k < r; ++k) {
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 2)*(2*k + 2)*gamma; flopCount += 3;
      }
      for(int k=1; k < r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k +3); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < r; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 2)*(2*k + 2)*alpha + (2*k + 1)*(2*k + 1)*beta;
	flopCount += 4;
      }
      for(int k=1; k < r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 3;
      }
    }
  }


  sigma[0] = 1;
  for(PetscInt k=1; k<size_d; ++k) {
    sigma[k] = sqrt(g[k-1]/f[k-1])*sigma[k-1];
  }

  double *M = (double*) calloc(sizeof(double), size_d*size_d);
  //fill diagonal with d
  for(int k=0; k < size_d; ++k)
    M[k*size_d + k]   = d[k];
  //fill above diagonal with g
  for(int k=0; k < size_d-1; ++k)
    M[(k)*size_d + k+1] = g[k];
  //fill below diagonal with f
  for(int k=0; k < size_d-1; ++k)
    M[(k+1)*size_d + k] = f[k];
    
  //scale entries by values of sigma to get symmetric matrix
  for(PetscInt i=0; i<size_d; ++i) {
    for(PetscInt j=0; j<size_d; ++j) {
      M[i*size_d+j] = M[i*size_d+j]*sigma[i]*(1./sigma[j]);
    }
  }
  if(n==5) {
    printf("\nSYMETRIC MATRIX:\n");
    for(PetscInt i=0; i<size_d; ++i) {
      for(PetscInt j=0; j<size_d; ++j) {
	printf("%4.4f ", M[i*size_d+j]);
      }
      printf("\n");
    }
  }
  
  //transpose M for lapack
  matTranspose(M, size_d);
  

  int lwork = 16*size_d;
  int info;
  double *work, *wr, *wi;
  *mat    = (double*) malloc(sizeof(double)*size_d*size_d);
  wr     = (double*) malloc(sizeof(double)*size_d);
  wi     = (double*) malloc(sizeof(double)*size_d);
  work   = (double*) malloc(sizeof(double)*lwork);

  //compute eigenvalues with lapack
  char complefteig, comprighteig;
  complefteig = 'N'; //don't compute left eigenvectors
  comprighteig = 'V'; //compute right eigenvectors
  
  /* #################################################### */
  /* ############## NEED TO ADD FLOP COUNT HERE ######### */
  dgeev_( &complefteig, &comprighteig, &size_d, M, &size_d, wr, wi, NULL, &size_d, *mat, &size_d, work, &lwork, &info );
  /* #################################################### */
  /* #################################################### */
  
  //matTranspose(vr, size_d);  
  //printf("hot\n");
  //printf("size_d: %d\n", size_d);
  free(wr);
  free(wi);
  free(d);
  free(f);
  free(g);
  free(work);
  free(M);
  //printf("amazing\n");
  //printf("diggity\n");
  if(mat_size != NULL)
    *mat_size = size_d;
  //for(int k=0; k<(m+1); ++k)
  //b[k] = b[k]/(b[(m+1)-1]/pow(-e->h2,(m+1)-1));


  //transform back to get eigenvectors
  for(PetscInt i=0; i<size_d; ++i) {
    for(PetscInt j=0; j<size_d; ++j) {
      (*mat)[i*size_d+j] = (*mat)[i*size_d+j]*1./(sigma[j]);
    }
  }
  
  for(int i=0; i<size_d; ++i) {
    for(int k=0; k<size_d; ++k) {
      (*mat)[i*size_d+k] = (*mat)[i*size_d+k]/((*mat)[(i*size_d) + size_d-1]/pow(-s->h2,(size_d)-1));
      flopCount += 3;
    }
  }
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  //return vr;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "getLameCoefficientMatrix"
PetscErrorCode getLameCoefficientMatrix(struct EllipsoidalSystem *s, char t, int n, int* const mat_size, double **mat)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  //For a given Lame type and ordr n, return the tridiagonal matrix
  //the size of the generated coefficient matrix is stored to mat_size
  //The eigenvectors of this matrix give the coefficients for the polynomaial P in the Lame function definition
  //The expressions come from Romain Annex 3
  //OPT: We could memorize this
  double alpha = s->h2;
  double beta = s->k2 - s->h2;
  double gamma = alpha - beta; flopCount += 2;
  int r = n/2;
  int size_d;
  double *d, *g, *f;
  if(t == 'K') {
    size_d = r+1;
  }
  else if(t == 'L' || t == 'M') {
    size_d = n-r;
  }
  else if(t == 'N') {
    size_d = r;
  }
  d = (double*) malloc(sizeof(double)*size_d);
  g = (double*) malloc(sizeof(double)*(size_d-1));
  f = (double*) malloc(sizeof(double)*(size_d-1));
  
  if (t == 'K') {
    //g missing last item
    for(int k = 0; k < r; ++k) {
      g[k] = -(2*k + 2)*(2*k + 1)*beta; flopCount++;
    }
    if(n%2) { //n is odd
      for(int k=0; k < r+1; ++k) {
	d[k] = ((2*r + 1)*(2*r + 2) - 4*k*k)*alpha + (2*k + 1)*(2*k + 1)*beta; flopCount += 6;
      }
      for(int k=1; k < r+1; ++k) {
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) + 1); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < r+1; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - 4*k*k*gamma; flopCount += 5;
      }
      for(int k=1; k < r+1; ++k) {
	f[k-1] = -alpha*(2*(r - k) + 2)*(2*(r + k) - 1); flopCount += 2;
      }
    }
  }
  else if (t == 'L') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 3)*beta; flopCount += 2;
      //printf("g[%d] = %15.15f\n", k, g[k]);
      //printf("n, r, %d, %d\n", n, r);
      //printf("beta: %15.15f\n", beta);
    }
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k) {
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 1)*(2*k + 1)*gamma; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1); flopCount += 3;
      }
    }
    else { //n is even
      for(int k=0; k < n-r; ++k) {
	d[k] = (2*r*(2*r + 1) - (2*k + 1)*(2*k + 1))*alpha + (2*k + 2)*(2*k + 2)*beta; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
  }
  else if (t == 'M') {
    //g missing last item
    for(int k = 0; k < n-r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 1)*beta; flopCount += 2;
    }
    if(n%2) { //n is odd
      for(int k=0; k < n-r; ++k) {
	d[k] = ((2*r + 1)*(2*r + 2) - (2*k + 1)*(2*k + 1))*alpha + 4*k*k*beta; flopCount += 5;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k + 2)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < n-r; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 1)*(2*k + 1)*gamma; flopCount += 3;
      }
      for(int k=1; k < n-r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 2;
      }
    }
  }
  else if (t == 'N') {
    //g missing last item
    for(int k = 0; k < r-1; ++k) {
      g[k] = -(2*k + 2)*(2*k + 3)*beta; flopCount += 2;
    }
    if(n%2) { //n is odd
      for(int k=0; k < r; ++k) {
	d[k] = (2*r + 1)*(2*r + 2)*alpha - (2*k + 2)*(2*k + 2)*gamma; flopCount += 3;
      }
      for(int k=1; k < r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k +3); flopCount += 2;
      }
    }
    else { //n is even
      for(int k=0; k < r; ++k) {
	d[k] = 2*r*(2*r + 1)*alpha - (2*k + 2)*(2*k + 2)*alpha + (2*k + 1)*(2*k + 1)*beta;
	flopCount += 4;
      }
      for(int k=1; k < r; ++k) {
	f[k-1] = -alpha*(2*r - 2*k)*(2*r + 2*k + 1); flopCount += 3;
      }
    }
  }


  double *M = (double*) calloc(sizeof(double), size_d*size_d);
  //fill diagonal with d
  for(int k=0; k < size_d; ++k)
    M[k*size_d + k]   = d[k];
  //fill above diagonal with g
  for(int k=0; k < size_d-1; ++k)
    M[(k)*size_d + k+1] = g[k];
  //fill below diagonal with f
  for(int k=0; k < size_d-1; ++k)
    M[(k+1)*size_d + k] = f[k];
    
  
  //transpose M for lapack
  matTranspose(M, size_d);
  

  int lwork = 16*size_d;
  int info;
  double *work, *wr, *wi;
  *mat    = (double*) malloc(sizeof(double)*size_d*size_d);
  wr     = (double*) malloc(sizeof(double)*size_d);
  wi     = (double*) malloc(sizeof(double)*size_d);
  work   = (double*) malloc(sizeof(double)*lwork);

  //compute eigenvalues with lapack
  char complefteig, comprighteig;
  complefteig = 'N'; //don't compute left eigenvectors
  comprighteig = 'V'; //compute right eigenvectors
  
  /* #################################################### */
  /* ############## NEED TO ADD FLOP COUNT HERE ######### */
  dgeev_( &complefteig, &comprighteig, &size_d, M, &size_d, wr, wi, NULL, &size_d, *mat, &size_d, work, &lwork, &info );
  /* #################################################### */
  /* #################################################### */
  
  //matTranspose(vr, size_d);  
  //printf("hot\n");
  //printf("size_d: %d\n", size_d);
  free(wr);
  free(wi);
  free(d);
  free(f);
  free(g);
  free(work);
  free(M);
  //printf("amazing\n");
  //printf("diggity\n");
  if(mat_size != NULL)
    *mat_size = size_d;
  //for(int k=0; k<(m+1); ++k)
  //b[k] = b[k]/(b[(m+1)-1]/pow(-e->h2,(m+1)-1));
  for(int i=0; i<size_d; ++i) {
    for(int k=0; k<size_d; ++k) {
      (*mat)[i*size_d+k] = (*mat)[i*size_d+k]/((*mat)[(i*size_d) + size_d-1]/pow(-s->h2,(size_d)-1));
      flopCount += 3;
    }
  }
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  //return vr;
  PetscFunctionReturn(0);
}

double *computeLameCoefficients(EllipsoidalSystem *e, int n, int p, int *vecsize) {
  char t = getLameTypeT(n, p);
  int tp = getLameTypeTp(n, p);
  int bsize;
  double *B;
  getLameCoefficientMatrix(e, t, n, &bsize, &B);
  //printf("B SIZE: %d\n", bsize);
  double *ret = (double*) malloc(sizeof(double)*bsize);
  memcpy(ret, B+tp*bsize, sizeof(double)*bsize);
  free(B);
  *vecsize = bsize;
  return ret;
}


#undef __FUNCT__
#define __FUNCT__ "initRomainConstsToOrderN"
PetscErrorCode initRomainConstsToOrderN(EllipsoidalSystem *e, int N)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  
  if(e->Rconsts != NULL)
    free(e->Rconsts);
  e->Rconsts = (double***) malloc(sizeof(double**)*(N+1));
  //int* mat_size; *mat_size = 1;
  for(int n=0; n<=N; ++n) {
    e->Rconsts[n] = (double**) malloc(sizeof(double*)*4);
    //e->Rconsts[n][0] = getLameCoefficientMatrix(e, 'K', n, NULL);
    ierr = getLameCoefficientMatrixSymmetric(e, 'K', n, NULL, e->Rconsts[n]+0);CHKERRQ(ierr);
    //printf("wowze\n");
    if(n != 0) {
      //e->Rconsts[n][1] = getLameCoefficientMatrix(e, 'L', n, NULL);
      ierr = getLameCoefficientMatrixSymmetric(e, 'L', n, NULL, e->Rconsts[n]+1);CHKERRQ(ierr);
      //e->Rconsts[n][2] = getLameCoefficientMatrix(e, 'M', n, NULL);
      ierr = getLameCoefficientMatrixSymmetric(e, 'M', n, NULL, e->Rconsts[n]+2);CHKERRQ(ierr);
      if(n != 1) {
	//e->Rconsts[n][3] = getLameCoefficientMatrix(e, 'N', n, NULL);
	ierr = getLameCoefficientMatrixSymmetric(e, 'N', n, NULL, e->Rconsts[n]+3);CHKERRQ(ierr);
      }
    }
  }

  ierr = PetscLogFlops(0);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "calcLame"
PetscErrorCode calcLame(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *Enp)
{
  PetscErrorCode ierr;
  PetscInt countFlops;
  PetscFunctionBegin;
  countFlops = 0;
  
  double signh, signk;
  if(signm*l >= 0)
    signh = 1;
  else
    signh = -1;
  if(signn*l >= 0)
    signk = 1;
  else
    signk = -1;
  char t = getLameTypeT(n, p);
  int tp = getLameTypeTp(n, p);


  //double *B = getLameCoefficientMatrix(e, t, n, &bsize);
  //printf("%c\n", t);

  int r = n/2;
  double l2 = l*l;
  double psi;
  int m;
  int type;
  //This comes from Romain Table II
  if(t == 'K') {
    m = r;
    psi = pow(l, n-2*r); countFlops++;
    type = 0;
  }
  else if(t == 'L') {
    m = n-r-1;
    psi = pow(l, 1-n+2*r)*signh*sqrt(fabs(l2-e->h2)); countFlops += 3;
    type = 1;
  }
  else if(t == 'M') {
    m = n-r-1;
    psi = pow(l, 1-n+2*r)*signk*sqrt(fabs(l2-e->k2)); countFlops += 3;
    type = 2;
  }
  else if(t == 'N') {
    m = r-1;
    psi = pow(l, n-2*r)*signh*signk*sqrt(fabs((l2-e->k2)*(l2-e->h2))); countFlops += 8;
    type = 3;
  }
  else
    printf("Invalid Lame Type\n");
  //if(isnan(psi))
  //printf("psi is nan and s: %15.15f\n", l);
  //Romain (30) \sum^m_{j=0} b_j (1 - \frac{\lambda^2}{h^2}) = \sum^m_{j=0} b_j \Lambda^j
  double Lambda_Romain = (1.0 - (l2/e->h2)); 

  //b = B[tp]
  //double *b = (double*) malloc(sizeof(double)*(m+1));
  double *b = e->Rconsts[n][type]+tp*(m+1); countFlops += 4;
  //for(int i=0; i<=m; ++i) {
  //memcpy(b, e->Rconsts[n][type]+tp*(m+1), sizeof(double)*(m+1));
    //}
  //memcpy(b, B+tp*(m+1), sizeof(double)*(m+1));
  
  //Normalize b so that the highest power of lambda has coefficient unity
  //Romain p242 and Dassios (D9,D10 and their gamma in B16-B20)
  //for(int k=0; k<(m+1); ++k)
  //b[k] = b[k]/(b[(m+1)-1]/pow(-e->h2,(m+1)-1));

  double P = b[m];
  for(int j=m-1; j>-1; --j) {
    P = P*Lambda_Romain + b[j];
    countFlops += 2;
  }
  //free(B);
  //free(b);
  //printf("psi: %15.15f\n", psi);
  //printf("p: %d\n", p);
  
  *Enp = psi*P; countFlops++;
  //printf("countFlops: %d\n", countFlops);
  
  ierr = PetscLogFlops(countFlops);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "calcLameDerivative"
PetscErrorCode calcLameDerivative(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *sol)
{
  PetscErrorCode ierr;
  PetscInt logFlops;
  PetscFunctionBegin;
  logFlops = 0;
  //Evaluate \frac{dE^p_n(l)}{dl}
  //The signs of \mu and \nu are necessary in order to determine the sign of \psi
  double signh, signk;
  if(signm*l >= 0)
    signh = 1;
  else
    signh = -1;
  if(signn*l >= 0)
    signk = 1;
  else
    signk = -1;
  
  char t = getLameTypeT(n, p);
  int tp = getLameTypeTp(n, p);
  
  int bsize;
  double *B;
  ierr = getLameCoefficientMatrix(e, t, n, &bsize, &B);CHKERRQ(ierr);
  
  int r = n/2;
  double l2 = l*l; logFlops ++;
  
  int m;
  double psi, psider;
  int type;
  //This comes from Romain Table II
  if(t == 'K') {
    m      = r;
    psi    = pow(l, n-2*r);
    psider = (n-2*r)*pow(l, n-2*r-1);
    type = 0;
    logFlops += 3;
  }
  else if(t == 'L') {
    m = n-r-1;
    psi = pow(l, 1-n+2*r)*signh*sqrt(fabs(l2 - e->h2)); logFlops += 4;
    psider = (1-n+2*r)*pow(l, -n+2*r)*signh*sqrt(fabs(l2 - e->h2)) + pow(l, 2-n+2*r)*signh/sqrt(fabs(l2 - e->h2)); logFlops += 12;
    type = 1;
  }
  else if(t == 'M') {
    m = n-r-1;
    psi    = pow(l, 1-n+2*r)*signk*sqrt(fabs(l2 - e->k2)); logFlops += 4;
    psider = (1-n+2*r)*pow(l, -n+2*r)*signk*sqrt(fabs(l2 - e->k2)) + pow(l, 2-n+2*r)*signk/sqrt(fabs(l2 - e->k2)); logFlops += 12;
    type = 2;
  }
  else if(t == 'N') {
    m = r-1;
    psi    = pow(l, n-2*r)*signh*signk*sqrt(fabs((l2 - e->k2)*(l2 - e->h2))); logFlops += 8;
    psider = ((n-2*r)*pow(l, n-2*r-1)*signh*signk*sqrt(fabs((l2 - e->k2)*(l2 - e->h2))) +
	      pow(l, n-2*r+1)*signh*signk*sqrt(fabs((l2 - e->k2)/(l2 - e->h2))) +
	      pow(l, n-2*r+1)*signh*signk*sqrt(fabs((l2 - e->h2)/(l2 - e->k2)))); logFlops += 24;
    type = 3;
  }
  else
    printf("Invalid Lame type\n");
  
  //Romain (30) \sum^m_{j=0} b_j (1 - \frac{\lambda^2}{h^2}) = \sum^m_{j=0} b_j \Lambda^j
  double Lambda_Romain = (1.0 - (l2/e->h2)); logFlops += 2; //Romain bottom of p.252
  //b = B[tp]
  double *b = (double*) malloc(sizeof(double)*(m+1));
  memcpy(b, e->Rconsts[n][type]+tp*(m+1), sizeof(double)*(m+1));
  double P = b[m];
  for(int j=m-1; j>-1; --j) {
    P = P*Lambda_Romain + b[j];
    logFlops += 2;
  }
  //P' = \sum^{m-1}_{k=0} c_k \Lambda^k where c_k = b_{k+1} (-2 (k+1) \frac{\lambda}{h^2}) 
  double Pder = b[m]*(-2*m*l/e->h2); logFlops += 4;
  free(b);
  free(B);

  *sol = psider*P + psi*Pder;
  
  ierr = PetscLogFlops(logFlops);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "integrand"
PetscErrorCode integrand(mpfr_t *x, mpfr_t *val, FuncInfo *ctx)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  mpfr_t *temp = &((*ctx).e->temp1);
  mpfr_d_div(*temp, 1.0, *x, MPFR_RNDN); flopCount ++;
  double s = mpfr_get_d(*temp, MPFR_RNDZ);
  double E;
  ierr = calcLame((*ctx).e, (*ctx).n, (*ctx).p, s, ctx->signm, ctx->signn, &E);CHKERRQ(ierr);
  mpfr_mul(*temp, *x, *x, MPFR_RNDN);
  mpfr_set(*val, *temp, MPFR_RNDN);
  mpfr_mul_d(*temp, *temp, (*ctx).e->k2, MPFR_RNDN);
  mpfr_mul_d(*val,  *val,  (*ctx).e->h2, MPFR_RNDN);
  mpfr_d_sub(*temp, 1.0, *temp, MPFR_RNDZ);
  mpfr_d_sub(*val,  1.0, *val,  MPFR_RNDZ);
  mpfr_sqrt(*temp, *temp, MPFR_RNDN);
  mpfr_sqrt(*val,  *val,  MPFR_RNDN);
  mpfr_mul(*val, *temp, *val, MPFR_RNDN);
  mpfr_mul_d(*val, *val, E, MPFR_RNDN);
  mpfr_mul_d(*val, *val, E, MPFR_RNDN); flopCount += 10;
  
  mpfr_d_div(*val, 1.0, *val, MPFR_RNDN); flopCount++;
  //printf("\n\n\nVAL: %15.15f\n\n\n", mpfr_get_d(*val, MPFR_RNDN));
  
  ierr = PetscLogFlops(flopCount);
  PetscFunctionReturn(0);
}

/*
  void integrand(mpfr_t *x , mpfr_t *val, FuncInfo *ctx)
  {
  double x2 = x*x;
  double s;
  int n = (*ctx).n;
  int p = (*ctx).p;
  
  if(x==0.0)
  s = 1e6;
  else
  s = 1.0/x;
  double E = calcLame((*ctx).e, n, p, s);
  return 1.0 / (E*E*sqrt(1.0 - (*ctx).e->k2*x2)*sqrt(1.0 - (*ctx).e->h2*x2));
  }
*/

#undef __FUNCT__
#define __FUNCT__ "calcI"
PetscErrorCode calcI(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *sol)
{
  PetscErrorCode ierr;
  PetscInt logFlops;
  PetscFunctionBegin;
  logFlops = 0;
  if(l < 0) {
    l *= -1.0;
    logFlops ++;
  }
  //int digits = 14;
  FuncInfo ctx = { e, n, p, signm, signn};
  mpfr_t *a = &(e->tempa);
  mpfr_t *b = &(e->tempb);
  mpfr_set_d(*a, 0.0, MPFR_RNDN);
  mpfr_set_d(*b, 1.0, MPFR_RNDN);
  mpfr_div_d(*b, *b, l, MPFR_RNDN); logFlops++;
  //printf("a: %15.15f\nb: %15.15f\n", mpfr_get_d(*a, MPFR_RNDN), mpfr_get_d(*b, MPFR_RNDN));
  int err = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*))integrand, e, *a, *b, 14, sol, &ctx);
  if(err)
    printf("integral failed\n");
  //printf("the integral is: %15.15f\n", *sol);

  ierr = PetscLogFlops(logFlops);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "calcIDerivative"
PetscErrorCode calcIDerivative(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *Ideriv)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  
  //Evaluate \frac{dI^p_n(l)}{dl} = \frac{-1}{(E^p_n(l))^2 \sqrt{l^2 - k^2}\sqrt{l^2 - h^2}}
  double l2 = l*l;
  double k2 = e->k2;
  double h2 = e->h2;
  double E;
  calcLame(e, n, p, l, signm, signn, &E);
  *Ideriv = -1.0 / (E*E*sqrt(l2 - k2)*sqrt(l2 - h2));

  ierr = PetscLogFlops(9);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "normFunction1"
PetscErrorCode normFunction1(mpfr_t *x, mpfr_t *val, FuncInfo2 *ctx)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  int n = (*ctx).n;
  int p = (*ctx).p;
  int numeratorType = (*ctx).numeratorType;
  int denomSign = (*ctx).denomSign;
  //double h2 = (*ctx).e->h2;
  //double k2 = (*ctx).e->k2;
  EllipsoidalSystem *e = (*ctx).e;
  double top;
  calcLame((*ctx).e, n, p, mpfr_get_d(*x, MPFR_RNDN), 1, 1, &top);
  mpfr_t *temp = &((*ctx).e->temp1);
  mpfr_mul(*temp, *x, *x, MPFR_RNDN);
  mpfr_d_sub(*val, e->k2, *temp, MPFR_RNDN);
  mpfr_sub_d(*temp, *temp, e->h2, MPFR_RNDN);
  mpfr_mul(*val, *temp, *val, MPFR_RNDN);
  mpfr_mul_d(*val, *val, (double) denomSign, MPFR_RNDN); flopCount += 5;
  if(mpfr_get_d(*val, MPFR_RNDN) < 0 && fabs(mpfr_get_d(*val, MPFR_RNDN)) < tol) {
    mpfr_mul_d(*val, *val, -1.0, MPFR_RNDN); flopCount++;
  }
  mpfr_sqrt(*temp, *val, MPFR_RNDN);
  mpfr_set_d(*val, top, MPFR_RNDN);
  mpfr_mul_d(*val, *val, top, MPFR_RNDN); flopCount += 2;
  flopCount += 3;
  if(numeratorType == 1) {
    mpfr_mul(*val, *val, *x, MPFR_RNDN);
    mpfr_mul(*val, *val, *x, MPFR_RNDN); flopCount += 2;
  }
  //printf("val before is %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));
  //printf("denom before is %8.8e\n", mpfr_get_d(*temp, MPFR_RNDN));
  mpfr_div(*val, *val, *temp, MPFR_RNDN); flopCount++;
  //printf("val is %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void normFunction2(mpfr_t *x, mpfr_t *val, FuncInfo3 *ctx)
{
  //printf("x: %8.8e\n", mpfr_get_d(*x, MPFR_RNDN));

  //double h2 = (*ctx).e->h2;
  //double k2 = (*ctx).e->k2;
  mpfr_t *temp = &((*ctx).e->temp1);
  EllipsoidalSystem *e = ctx->e;
  mpfr_set(*temp, e->hp_k2, MPFR_RNDN);
  mpfr_sub(*temp, *temp, e->hp_h2, MPFR_RNDN);
  mpfr_mul(*val, *x, e->hp_h2, MPFR_RNDN);
  mpfr_add(*val, *val, *temp, MPFR_RNDN);
  mpfr_sqrt(*val, *val, MPFR_RNDN);
  
  //printf("val 1: %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));

  mpfr_abs(*temp, *x, MPFR_RNDN);
  //printf("absx(actual): %8.8e\n", mpfr_get_d(*temp, MPFR_RNDN));
  mpfr_sqrt(*temp, *temp, MPFR_RNDN);
  //printf("before multiply val 2: %8.8e\n", mpfr_get_d(*temp, MPFR_RNDN));
  //printf("absx: %8.8e\n", mpfr_get_d(*temp, MPFR_RNDN));
  mpfr_mul(*val, *val, *temp, MPFR_RNDN);

  //printf("val 2: %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));

  mpfr_set_d(*temp, 1.0, MPFR_RNDN);
  mpfr_sub(*temp, *temp, *x, MPFR_RNDN);
  mpfr_sqrt(*temp, *temp, MPFR_RNDN);
  //printf("before multiply val 3: %8.8e\n", mpfr_get_d(*temp, MPFR_RNDN));
  mpfr_mul(*val, *val, *temp, MPFR_RNDN);
  //printf("val 3(denom): %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));
  if((*ctx).numeratorType == 0) {
    //printf("numerator is 1.0\n");
    mpfr_d_div(*val, 1.0, *val, MPFR_RNDN);
  }
  else
  {
    mpfr_div(*val, *x, *val, MPFR_RNDN);
    //printf("numerator is x\n");
  }
  //printf("val: %8.8e\n", mpfr_get_d(*val, MPFR_RNDN));
  
}


#undef __FUNCT__
#define __FUNCT__ "calcNormalization"
PetscErrorCode calcNormalization(EllipsoidalSystem *e, int n, int p, double *normConst)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  //double h = e->h;
  //double h2 = e->h2;
  //double k = e->k;
  //double k2 = e->k2;
  //The four integrals making up Romain (51)
  //structs containing f unction info
  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };
  int err[4];
  double integrals[4];

  mpfr_t *mpfrzero = &(e->mpfrzero);
  mpfr_t *mpfrone  = &(e->mpfrone);
  
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, e->hp_h, e->hp_k, 16, integrals, &ctx1);

  err[1] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, e->hp_h, e->hp_k, 16, integrals+1, &ctx2);

  err[2] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, *mpfrzero, e->hp_h, 16, integrals+2, &ctx3);

  err[3] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, *mpfrzero, e->hp_h, 16, integrals+3, &ctx4);

  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]);


  ierr = PetscLogFlops(4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
  /*
  for(int i=0; i<4; ++i) {
    //printf("integrals[i] = %8.8e\n", i, integrals[i]);
    if(err[i])
      printf("Integral %d failed in calcNormalization\n", i);
  }

  //The basic elliptic integrals Romain (53)
  //He has an error in the equation, the prefactor should be h/2
  //structs containing function info
  FuncInfo3 ctx20 = { .e = e, .numeratorType = 0, .denomSign = -1 };
  FuncInfo3 ctx21 = { .e = e, .numeratorType = 1, .denomSign = -1 };
  FuncInfo3 ctx30 = { .e = e, .numeratorType = 0, .denomSign = 1 };
  FuncInfo3 ctx31 = { .e = e, .numeratorType = 1, .denomSign = 1 };
  double integrals2[4];

  mpfr_t *endpt = &(e->endpt);
  mpfr_set(*endpt, e->hp_k2, MPFR_RNDN);
  mpfr_div(*endpt, *endpt, e->hp_h2, MPFR_RNDN);
  mpfr_d_sub(*endpt, 1.0, *endpt, MPFR_RNDN);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *endpt, 14, integrals2, &ctx20);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *endpt, 14, integrals2+1, &ctx21);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *mpfrone, 14, integrals2+2, &ctx30);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *mpfrone, 14, integrals2+3, &ctx31);

  //Multiply by prefactor h/2
  for(int i=0; i<4; ++i) {
    //integrals2[i] = .5*h*integrals2[i];
  }


  for(int i=0; i<4; ++i) {
    if(err[i])
      printf("Integral2 %d failed in calcNormalization\n", i);
  }
  
  //Solve system for coefficients
  double *M = (double*) malloc(sizeof(double)*4);
  M[0] = integrals2[0]; //integrals2[0]; //transposed for lapack
  M[1] = integrals2[2]; //
  M[2] = integrals2[1]; //
  M[3] = integrals2[3]; //

  double *sol = (double*) malloc(sizeof(double)*2);
  sol[0] = integrals[0];
  sol[1] = integrals[2];

  //variables to store solutions
  double alpha, beta, A, B;

  char tchar = 'N';
  int two = 2;
  int one = 1;

  //ipiv pivot (output for dgetrf, input for dgetrs)
  int *ipiv = (int*) malloc(sizeof(int)*2);
  int info;

  //first solve
  dgetrf_(&two, &two, M, &two, ipiv, &info);
  dgetrs_(&tchar, &two, &one, M, &two, ipiv, sol, &two, &info);
  alpha = sol[0];
  beta  = sol[1];

  sol[0] = integrals[1];
  sol[1] = integrals[3];

  //second solve
  dgetrs_(&tchar, &two, &one, M, &two, ipiv, sol, &two, &info);
  A = sol[0];
  B = sol[1];
  free(sol);
  free(M);
  free(ipiv);
  //Romain (54)
  //The factor 8 is from Dassios (see discussion in Deng), since we integrate over the whole ellipsoid rather than just one octant
  if(fabs(alpha*B - beta*A) == 0) {
    printf("\n\nbad\n\n\n\n");
    printf("alpha*B: %16.16e\n", alpha*B);
    printf("beta*A: %16.16e\n", beta*A);
    return 8 * (PETSC_PI/2.0)*1e-14;
  }
  return 8 * (PETSC_PI/2.0)*(alpha*B - beta*A);
  */
}

#undef __FUNCT__
#define __FUNCT__ "calcNormalization2"
PetscErrorCode calcNormalization2(EllipsoidalSystem *e, PetscInt n, PetscInt p, double *intVals, double *normConst)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  //double h = e->h;
  //double h2 = e->h2;
  //double k = e->k;
  //double k2 = e->k2;
  //The four integrals making up Romain (51)
  //structs containing f unction info
  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };
  int err[4];
  double integrals[4];

  mpfr_t *mpfrzero = &(e->mpfrzero);
  mpfr_t *mpfrone  = &(e->mpfrone);
  
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, e->hp_h, e->hp_k, 16, integrals, &ctx1);
  intVals[0] = integrals[0];
  err[1] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, e->hp_h, e->hp_k, 16, integrals+1, &ctx2);
  intVals[1] = integrals[1];
  err[2] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, *mpfrzero, e->hp_h, 16, integrals+2, &ctx3);
  intVals[2] = integrals[2];
  err[3] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e, *mpfrzero, e->hp_h, 16, integrals+3, &ctx4);
  intVals[3] = integrals[3];
  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]);


  ierr = PetscLogFlops(4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
  /*
  for(int i=0; i<4; ++i) {
    //printf("integrals[i] = %8.8e\n", i, integrals[i]);
    if(err[i])
      printf("Integral %d failed in calcNormalization\n", i);
  }

  //The basic elliptic integrals Romain (53)
  //He has an error in the equation, the prefactor should be h/2
  //structs containing function info
  FuncInfo3 ctx20 = { .e = e, .numeratorType = 0, .denomSign = -1 };
  FuncInfo3 ctx21 = { .e = e, .numeratorType = 1, .denomSign = -1 };
  FuncInfo3 ctx30 = { .e = e, .numeratorType = 0, .denomSign = 1 };
  FuncInfo3 ctx31 = { .e = e, .numeratorType = 1, .denomSign = 1 };
  double integrals2[4];

  mpfr_t *endpt = &(e->endpt);
  mpfr_set(*endpt, e->hp_k2, MPFR_RNDN);
  mpfr_div(*endpt, *endpt, e->hp_h2, MPFR_RNDN);
  mpfr_d_sub(*endpt, 1.0, *endpt, MPFR_RNDN);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *endpt, 14, integrals2, &ctx20);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *endpt, 14, integrals2+1, &ctx21);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *mpfrone, 14, integrals2+2, &ctx30);
  err[0] = integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction2, e, *mpfrzero, *mpfrone, 14, integrals2+3, &ctx31);

  //Multiply by prefactor h/2
  for(int i=0; i<4; ++i) {
    //integrals2[i] = .5*h*integrals2[i];
  }


  for(int i=0; i<4; ++i) {
    if(err[i])
      printf("Integral2 %d failed in calcNormalization\n", i);
  }
  
  //Solve system for coefficients
  double *M = (double*) malloc(sizeof(double)*4);
  M[0] = integrals2[0]; //integrals2[0]; //transposed for lapack
  M[1] = integrals2[2]; //
  M[2] = integrals2[1]; //
  M[3] = integrals2[3]; //

  double *sol = (double*) malloc(sizeof(double)*2);
  sol[0] = integrals[0];
  sol[1] = integrals[2];

  //variables to store solutions
  double alpha, beta, A, B;

  char tchar = 'N';
  int two = 2;
  int one = 1;

  //ipiv pivot (output for dgetrf, input for dgetrs)
  int *ipiv = (int*) malloc(sizeof(int)*2);
  int info;

  //first solve
  dgetrf_(&two, &two, M, &two, ipiv, &info);
  dgetrs_(&tchar, &two, &one, M, &two, ipiv, sol, &two, &info);
  alpha = sol[0];
  beta  = sol[1];

  sol[0] = integrals[1];
  sol[1] = integrals[3];

  //second solve
  dgetrs_(&tchar, &two, &one, M, &two, ipiv, sol, &two, &info);
  A = sol[0];
  B = sol[1];
  free(sol);
  free(M);
  free(ipiv);
  //Romain (54)
  //The factor 8 is from Dassios (see discussion in Deng), since we integrate over the whole ellipsoid rather than just one octant
  if(fabs(alpha*B - beta*A) == 0) {
    printf("\n\nbad\n\n\n\n");
    printf("alpha*B: %16.16e\n", alpha*B);
    printf("beta*A: %16.16e\n", beta*A);
    return 8 * (PETSC_PI/2.0)*1e-14;
  }
  return 8 * (PETSC_PI/2.0)*(alpha*B - beta*A);
  */
}

//calculate the eigenvalues of the surface operator
#undef __FUNCT__
#define __FUNCT__ "calcSurfaceOperatorEigenvalues"
PetscErrorCode calcSurfaceOperatorEigenvalues(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *sol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  double a = e->a;
  double b = e->b;
  double c = e->c;
  double Inp;
  ierr = calcI(e, n, p, l, signm, signn, &Inp);CHKERRQ(ierr);
  double Enp;
  ierr = calcLame(e, n, p, l, signm, signn, &Enp);CHKERRQ(ierr);
  double EnpDer;
  calcLameDerivative(e, n, p, l, signm, signn, &EnpDer);

  *sol = (2.0*a*b*c*(EnpDer/a)*Inp*Enp - 1)/2;
  
  ierr = PetscLogFlops(9);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "integrateMPFR"
PetscErrorCode integrateMPFR(PetscErrorCode (*f)(mpfr_t *,mpfr_t*,void*), EllipsoidalSystem *e, mpfr_t a, mpfr_t b, int digits, double *integral, void *ctx) 
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  //printf("can we get out of here?\n");
  const int  safetyFactor = 2;    /* Calculate abcissa until 2*p digits */
  int        p      = e->precision;   /* Digits of precision in the evaluation */
  int        l      = 0;    /* Level of refinement, h = 2^{-l} */
  mpfr_t     *alpha = &(e->alpha);    /* Half-width of the integration interval */
  mpfr_t     *beta  = &(e->beta);                /* Center of the integration interval */
  mpfr_t     *h     = &(e->h_step);                   /* Step size, length between x_k */
  mpfr_t     *osum  = &(e->osum);                /* Integral on last level */
  mpfr_t     *psum  = &(e->psum);                /* Integral on the level before the last level */
  mpfr_t     *sum   = &(e->sum);                 /* Integral on current level */
  mpfr_t     *yk    = &(e->yk);                  /* Quadrature point 1 - x_k on reference domain [-1, 1] */
  mpfr_t     *lx    = &(e->lx);
  mpfr_t     *rx    = &(e->rx);              /* Quadrature points to the left and right of 0 on the real domain [a, b] */
  mpfr_t     *wk    = &(e->wk);                  /* Quadrature weight at x_k */
  mpfr_t     *lval  = &(e->lval);
  mpfr_t     *rval  = &(e->rval);          /* Terms in the quadature *sum to the left and right of 0 */
  mpfr_t     *pi2   = &(e->pi2);
  mpfr_t     *kh    = &(e->kh);
  mpfr_t     *msinh = &(e->msinh);
  mpfr_t     *mcosh = &(e->mcosh);
  mpfr_t   *maxTerm = &(e->maxTerm);
  mpfr_t   *curTerm = &(e->curTerm);
  mpfr_t     *tmp   = &(e->tmp);
  int        d;                   /* Digits of precision in the integral */
  int counter = 1;

  if(digits <= 0) {
    printf("Please give a positive number of significant digits\n");
    return 1;
  }
  digits = 16;
  /* Create high precision storage */


  /* Initialization */
  mpfr_set(*alpha, b, MPFR_RNDN);
  mpfr_sub(*alpha, *alpha, a, MPFR_RNDN);
  mpfr_mul_d(*alpha, *alpha, 0.5, MPFR_RNDN);
  mpfr_set(*beta, b, MPFR_RNDN);
  mpfr_add(*beta, *beta, a, MPFR_RNDN);
  mpfr_mul_d(*beta, *beta, 0.5, MPFR_RNDN);
  mpfr_set_d(*osum,  0.0, MPFR_RNDN);
  mpfr_set_d(*psum,  0.0, MPFR_RNDN);
  mpfr_set_d(*h,     1.0, MPFR_RNDN);
  mpfr_const_pi(*pi2, MPFR_RNDN);
  mpfr_mul_d(*pi2, *pi2, 0.5, MPFR_RNDN); flopCount += 5;

  /* Center term */
  //mpfr_set_d(*lx, 0.5*(b+a), MPFR_RNDN);
  mpfr_set(*lx, b, MPFR_RNDN);
  mpfr_add(*lx, *lx, a, MPFR_RNDN);
  mpfr_mul_d(*lx, *lx, .5, MPFR_RNDN);
  f(lx, sum, ctx);
  mpfr_mul(*sum, *sum, *alpha, MPFR_RNDN);
  mpfr_mul(*sum, *sum, *pi2, MPFR_RNDN); flopCount += 4;


  //DELETE LATER
  int maxL = 6;
  double *SUMS = (double*) malloc(sizeof(double)*maxL);
  int *insideSums = (int*) malloc(sizeof(int)*maxL);


  /* */
  do {
    double d1, d2, d3, d4;
    int  k = 1;

    ++l;

    mpfr_set_d(*maxTerm, 0.0, MPFR_RNDN);
    /* PetscPrintf(PETSC_COMM_SELF, "LEVEL %D sum: %15.15f\n", l, sum); */
    //printf("LEVEL %d sem: %15.15f\n", l, sum);
    /* At each level of refinement, h --> h/2 and sum --> sum/2 */
    mpfr_set(*psum, *osum, MPFR_RNDN);
    mpfr_set(*osum, *sum, MPFR_RNDN);
    mpfr_mul_d(*h,  *h,  0.5, MPFR_RNDN);
    mpfr_mul_d(*sum, *sum, 0.5, MPFR_RNDN); flopCount += 2;

    do {
      //printf("doot\n");
      //if(mpfr_number_p(sum) == 0)
	//printf("at iterations %d sum is %8.8e a: %8.8f b: %8.8f\n", k, mpfr_get_d(sum, MPFR_RNDN), a, b);
      mpfr_set_si(*kh, k, MPFR_RNDN);
      mpfr_mul(*kh, *kh, *h, MPFR_RNDN);
      /* Weight */
      mpfr_set(*wk, *h, MPFR_RNDN);
      mpfr_sinh_cosh(*msinh, *mcosh, *kh, MPFR_RNDN);
      mpfr_mul(*msinh, *msinh, *pi2, MPFR_RNDN);
      mpfr_mul(*mcosh, *mcosh, *pi2, MPFR_RNDN);
      mpfr_cosh(*tmp, *msinh, MPFR_RNDN);
      mpfr_sqr(*tmp, *tmp, MPFR_RNDN);
      mpfr_mul(*wk, *wk, *mcosh, MPFR_RNDN);
      mpfr_div(*wk, *wk, *tmp, MPFR_RNDN); flopCount += 8;
      /* Abscissa */
      mpfr_set_d(*yk, 1.0, MPFR_RNDZ);
      mpfr_cosh(*tmp, *msinh, MPFR_RNDN);
      mpfr_div(*yk, *yk, *tmp, MPFR_RNDZ);
      mpfr_exp(*tmp, *msinh, MPFR_RNDN);
      mpfr_div(*yk, *yk, *tmp, MPFR_RNDZ); flopCount += 4;
      /* Quadrature points */
      mpfr_sub_d(*lx, *yk, 1.0, MPFR_RNDZ);
      mpfr_mul(*lx, *lx, *alpha, MPFR_RNDU);
      mpfr_add(*lx, *lx, *beta, MPFR_RNDU);
      mpfr_d_sub(*rx, 1.0, *yk, MPFR_RNDZ);
      mpfr_mul(*rx, *rx, *alpha, MPFR_RNDD);
      mpfr_add(*rx, *rx, *beta, MPFR_RNDD); flopCount += 4;
      //printf("poot\n");
      /* Evaluation */
      f(lx, lval, ctx);
      //printf("lval: %8.8e\n", mpfr_get_d(lval, MPFR_RNDN));
      f(rx, rval, ctx);
      //printf("rval: %8.8e\n", mpfr_get_d(rval, MPFR_RNDN));
      //check if function evaluates to inf or nan
      int update = 1;
      if(mpfr_number_p(*lval) == 0) {
	update = 0;
	//printf("lval is nan/inf *lx: %15.15e\n", mpfr_get_d(*lx, MPFR_RNDN));
      }
      if(mpfr_number_p(*rval) == 0) {
	update = 0;
	//printf("rval is nan/inf\n");
      }
      /* Update if function not inf or nan*/
      if(update == 1) {
	mpfr_mul(*tmp, *wk, *alpha, MPFR_RNDN);
	mpfr_mul(*tmp, *tmp, *lval, MPFR_RNDN);
	mpfr_add(*sum, *sum, *tmp, MPFR_RNDN);
	mpfr_abs(*tmp, *tmp, MPFR_RNDN);
	mpfr_max(*maxTerm, *maxTerm, *tmp, MPFR_RNDN);
	mpfr_set(*curTerm, *tmp, MPFR_RNDN);
	mpfr_mul(*tmp, *wk, *alpha, MPFR_RNDN);
	mpfr_mul(*tmp, *tmp, *rval, MPFR_RNDN);
	mpfr_add(*sum, *sum, *tmp, MPFR_RNDN);
	mpfr_abs(*tmp, *tmp, MPFR_RNDN);
	mpfr_max(*maxTerm, *maxTerm, *tmp, MPFR_RNDN);
	mpfr_max(*curTerm, *curTerm, *tmp, MPFR_RNDN); flopCount += 8;
	counter += 2;
      }
      /* if (l == 1) printf("k is %d and sum is %15.15f and *wk is %15.15f\n", k, sum, *wk); */
      ++k;
      /* Only need to evaluate every other point on refined levels */
      if (l != 1) ++k;
      mpfr_log10(*tmp, *wk, MPFR_RNDN); flopCount++;
      mpfr_abs(*tmp, *tmp, MPFR_RNDN);
    } while (mpfr_get_d(*tmp, MPFR_RNDN) < safetyFactor*p); /* Only need to evaluate sum until weights are < 32 digits\
							      of precision */
    SUMS[l-1] = mpfr_get_d(*sum, MPFR_RNDN);
    insideSums[l-1] = counter;
    mpfr_sub(*tmp, *sum, *osum, MPFR_RNDN);
    mpfr_abs(*tmp, *tmp, MPFR_RNDN);
    mpfr_log10(*tmp, *tmp, MPFR_RNDN);
    d1 = mpfr_get_d(*tmp, MPFR_RNDN);
    mpfr_sub(*tmp, *sum, *psum, MPFR_RNDN);
    mpfr_abs(*tmp, *tmp, MPFR_RNDN);
    mpfr_log10(*tmp, *tmp, MPFR_RNDN);
    d2 = mpfr_get_d(*tmp, MPFR_RNDN);
    mpfr_log10(*tmp, *maxTerm, MPFR_RNDN);
    d3 = mpfr_get_d(*tmp, MPFR_RNDN) - p;
    mpfr_log10(*tmp, *curTerm, MPFR_RNDN);
    d4 = mpfr_get_d(*tmp, MPFR_RNDN);
    d  = (int) fabs(min(0, max(max(max((d1*d1)/d2, 2*d1), d3), d4))); flopCount += 9;

  } while (d < digits && l < maxL);

  *integral = mpfr_get_d(*sum, MPFR_RNDN);
  /* Cleanup */
  //mpfr_clears(*alpha, *beta, h, sum, *osum, *psum, *yk, *wk, *lx, *rx, *tmp, *maxTerm, *curTerm, *pi2, *kh, *msinh, *mcosh, lval, rval, NULL);
  //printf("yes we can\n");
  //printf("levels: %d\n", l);
  /*
  if(l == maxL) {
    for(int i=0; i<l; ++i) {
      printf("SUMS[%d] = %15.15f\n", i, SUMS[i]);
      printf("inside total = %d\n", insideSums[i]);
    }
  }
  */
  //printf("total is: %d\n", counter);
  free(SUMS);

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "integrateMidpoint"
PetscErrorCode integrateMidpoint(PetscErrorCode (*f)(mpfr_t *,mpfr_t*,void*), mpfr_t a, mpfr_t b, int n, double *integral, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;
  
  int digs = 32;

  mpfr_t delT, sum, fval, xval;
  mpfr_init2(delT, 4*digs);
  mpfr_init2(sum,  4*digs);
  mpfr_init2(fval, 4*digs);
  mpfr_init2(xval, 4*digs);


  mpfr_set_d(sum, 0.0, MPFR_RNDN);

  //calculate delT
  mpfr_sub(delT, b, a, MPFR_RNDN);
  mpfr_div_ui(delT, delT, n, MPFR_RNDN);
  //start xval at a+.5*delT
  mpfr_mul_d(xval, delT, 0.5, MPFR_RNDN);
  mpfr_add(xval, a, xval, MPFR_RNDN); flopCount += 4;

  for(int i=0; i < n; ++i) {

    *integral = mpfr_get_d(sum, MPFR_RNDN);
    
    //calculate function value
    f(&xval, &fval, ctx);
    
    //increment sum
    mpfr_add(sum, sum, fval, MPFR_RNDN);
    
    //increment xval
    mpfr_add(xval, xval, delT, MPFR_RNDN); flopCount += 2;

  }
  mpfr_mul(sum, sum, delT, MPFR_RNDN); flopCount ++;
  
  *integral = mpfr_get_d(sum, MPFR_RNDN);

  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
