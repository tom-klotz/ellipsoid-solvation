#ifndef __ellipsoid
#define __ellipsoid

static double tol = 1e-12;

//structs
typedef struct EllipsoidalSystem {
  double a;
  double b;
  double c;

  double ***Dconsts;
  double **normConstants;
  int DmaxN;
  double ***Rconsts;
  int RmaxN;
  char **tVals;
  int **tpVals;

  double h2, h;
  mpfr_t hp_h2, hp_h;
  double k2, k;
  mpfr_t hp_k2, hp_k;

  //precision used for mpfr calculations
  int precision;

  //mpfr variables for integration function
  mpfr_t     alpha;               /* Half-width of the integration interval */
  mpfr_t     beta;                /* Center of the integration interval */
  mpfr_t     h_step;                   /* Step size, length between x_k */
  mpfr_t     osum;                /* Integral on last level */
  mpfr_t     psum;                /* Integral on the level before the last level */
  mpfr_t     sum;                 /* Integral on current level */
  mpfr_t     yk;                  /* Quadrature point 1 - x_k on reference domain [-1, 1] */
  mpfr_t     lx, rx;              /* Quadrature points to the left and right of 0 on the real domain [a, b] */
  mpfr_t     wk;                  /* Quadrature weight at x_k */
  mpfr_t     lval, rval;          /* Terms in the quadature sum to the left and right of 0 */
  mpfr_t     pi2;
  mpfr_t     kh;
  mpfr_t     msinh;
  mpfr_t     mcosh;
  mpfr_t     maxTerm;
  mpfr_t     curTerm;
  mpfr_t     tmp;
  
  //temporary mpfr variables for use
  mpfr_t temp1;
  mpfr_t temp2;
  mpfr_t temp3;
  mpfr_t tempa;
  mpfr_t tempb;

  //other used mpfr variables
  mpfr_t endpt;
  
  //mpfr constants
  mpfr_t mpfrzero;
  mpfr_t mpfrone;

} EllipsoidalSystem;
typedef struct Point {
  double x1;
  double x2;
  double x3;
  
  char type;
} Point;
typedef struct FuncInfo {
  EllipsoidalSystem *e;
  int n;
  int p;
  int signm;
  int signn;
} FuncInfo;
typedef struct FuncInfo2 {
  EllipsoidalSystem *e;
  int n;
  int p;
  int numeratorType;
  int denomSign;
} FuncInfo2;
typedef struct FuncInfo3 {
  EllipsoidalSystem *e;
  int numeratorType;
  int denomSign;
} FuncInfo3;
typedef struct FuncInfo4 {
  EllipsoidalSystem *e;
  double a2;
  double b2;
  double c2;
  double botVar;
} FuncInfo4;


//petsc functions
PetscErrorCode CartesianToEllipsoidalVec(EllipsoidalSystem*, Vec, Vec);
PetscErrorCode ellipsoidToCartesian(struct EllipsoidalSystem *s, struct Point *p);
PetscErrorCode getCoefsK(EllipsoidalSystem *s, int n, double **coefs);
PetscErrorCode calcLame2(EllipsoidalSystem *s, int n, int p, double l, double *sol);
PetscErrorCode calcI(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *sol);
PetscErrorCode calcSurfaceOperatorEigenvalues(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn, double *sol);
PetscErrorCode calcLame(struct EllipsoidalSystem *s, int n, int p, double l, int signm, int signn, double *Enp);

//functions
void initEllipsoidalSystem(struct EllipsoidalSystem *s, double a, double b, double c);
void initRomainConstsToOrderN(EllipsoidalSystem *e, int N);
void ellipsoidInitToOrderN(struct EllipsoidalSystem *s, int N);
//void getCoefsK(EllipsoidalSystem *s, int n, double **coefs);
void getCoefsL(EllipsoidalSystem *s, int n, double **coefs);
void getCoefsM(EllipsoidalSystem *s, int n, double **coefs);
void getCoefsN(EllipsoidalSystem *s, int n, double **coefs);
//double calcLame2(EllipsoidalSystem *s, int n, int p, double l);
//void ellipsoidToCartesian(struct EllipsoidalSystem *s, struct Point *p);
void cartesianToEllipsoidal(struct EllipsoidalSystem *s, struct Point *p);
void cartesianToEllipsoidal2(struct EllipsoidalSystem *e, struct Point *p);
char getLameTypeT(int n, int p);
int getLameTypeTp(int n, int p);
double *getLameCoefficientMatrix(struct EllipsoidalSystem *s, char t, int n, int *mat_size);
double *computeLameCoefficients(struct EllipsoidalSystem *s, int n, int p, int *vecsize);
//double calcLame(struct EllipsoidalSystem *s, int n, int p, double l, int signm, int signn);
double calcLameDerivative(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn);
void integrand(mpfr_t *x, mpfr_t *val, FuncInfo *ctx);
void normFunction1(mpfr_t *x, mpfr_t *val, FuncInfo2 *ctx);
void normFunction2(mpfr_t *x, mpfr_t *val, FuncInfo3 *ctx);
//double calcI(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn);
double calcIDerivative(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn);
double calcNormalization(EllipsoidalSystem *e, int n, int p);
//double calcSurfaceOperatorEigenvalues(EllipsoidalSystem *e, int n, int p, double l, int signm, int signn);
int integrate(void (*f)(mpfr_t*, mpfr_t*, void*), double a, double b, int digits, double *integral, void *ctx);
int integrateMPFR(void (*f)(mpfr_t *,mpfr_t*,void*), EllipsoidalSystem *e, mpfr_t a, mpfr_t b, int digits, double *integral, void *ctx);
int integrateMidpoint(void (*f)(mpfr_t *,mpfr_t*,void*), mpfr_t a, mpfr_t b, int digits, double *integral, void *ctx);
#endif
