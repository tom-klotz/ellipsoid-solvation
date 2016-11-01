//contains old versions of functions/prototypes, mostly for reference. Try not to use this file directly


//old, was in speedtester.c
int main(int argc, char **argv)
{
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3, 2, 1);
  
  /*
  Point x = {.x1 = 0.1, .x2 = 0.2, .x3 = 0.3};
  cartesianToEllipsoidal(&e, &x);
  printf("1: %15.15f\n2: %15.15f\n3: %15.15f\n", x.x1, x.x2, x.x3);
  int num = 20;
  FILE *fp = fopen("plotlame.txt", "w");
  for(int n=0; n<8; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      for(int i=0; i<1; ++i) {
	fprintf(fp, "%d %d %15.15f\n", n, p, calcLame(&e, n, p, .35+((double)i)/n));
      }
    }
  }
  */
  printf("calcNorm: %15.15f\n", calcNormalization(&e, 1, 2));
  return 1;
}

void testSphericalCoordinates() {

  const int numRandPoints = 20;

  SPoint testPoints[numRandPoints];
  SPoint testPointsCpy[numRandPoints];
  
  srand(time(NULL));
  for(int i=0; i < numRandPoints; ++i) {
    testPoints[i].x1 = (double)rand()/(double)RAND_MAX;
    testPoints[i].x2 = (double)rand()/(double)RAND_MAX;
    testPoints[i].x3 = (double)rand()/(double)RAND_MAX;
    testPointsCpy[i].x1 = testPoints[i].x1;
    testPointsCpy[i].x2 = testPoints[i].x2;
    testPointsCpy[i].x3 = testPoints[i].x3;
  }
  
  //convert to spherical and back to cartesian and print error
  for(int i=0; i < numRandPoints; ++i) {
    cartesianToSpherical(testPoints+i);
    sphericalToCartesian(testPoints+i);

    printf("error[%d].x1 = %8.8e\n", i, fabs(testPoints[i].x1-testPointsCpy[i].x1));
    printf("error[%d].x2 = %8.8e\n", i, fabs(testPoints[i].x2-testPointsCpy[i].x2));
    printf("error[%d].x3 = %8.8e\n", i, fabs(testPoints[i].x3-testPointsCpy[i].x3));
  }
 
}


double calcRegion2(SProblem *problem, int maxDeg, SPoint *r) {
  double val = 0;
  double Enm, Cnm;
  double e1 = problem->e1;
  double e2 = problem->e2;
  double rad = r->x1;
  double theta = r->x2;
  double phi = r->x3;
  double radpow = 1;
  double lterm;
  double rterm, iterm;
  for(int n=0; n<=maxDeg; ++n) {
    radpow *= rad;
    for(int m=-n; m<=n; ++m) {
      calcEnm(problem->positions, problem->charges, problem->nCharges, n, m, &rterm, &iterm);
      //Cnm = calcCnmFromEnm(m, n, rterm, problem);
      Cnm = 0;
      if(fabs(Enm-Cnm) > 1e-14)
	printf("\n\n\nWHATWHATWHAT\n\n\n");
      lterm = gsl_sf_legendre_Plm(n, m, cos(theta)); //boost::math::legendre_p(n, m, cos(theta));
      //if(m%2 == 1)
      //lterm *= -1;
      val += Cnm*lterm/radpow;
    }
    printf("n: %d\nval: %15.15f\n\n", n, val);

  }
  return val;
}

//output leg is vector of length nz*(2*n+1)
void legendre2(int n, int nz, double z, double *leg)
{
  double fac = 1.0;
  double sqz2   = sqrt(1.0 - z*z);
  double hsqz2  = 0.5*sqz2;
  double ihsqz2 = z/hsqz2;
  for(int i=2; i<=n; ++i)
    fac *= i;
  double hsqz2ton = 1.0;
  for(int i=0; i<n; ++i)
    hsqz2ton *= hsqz2;
  int pre = n % 2 ? -1 : 1;

  
  if(sqz2 == 0) {
    //printf("WOW!\n");
    for(int mr = 0; mr < 2*n+1; ++mr)
      leg[mr] = gsl_sf_legendre_Plm(n, mr-n, z); //boost::math::legendre_p(n, mr-n, z);
  }
  else if(n==0) {
    leg[0] = 1.0;
  }
  else if(n==1) {
    leg[0] = -hsqz2;
    leg[1] = z;
    leg[2] = sqz2;
  }
  else {
    leg[0] = (1.0 - 2.0*fabs(n - 2.0*floor(n/2.))) * hsqz2ton/fac;
    leg[1] = -leg[0]*n*ihsqz2;
    for(int mr = 1; mr < 2*n; ++mr)
      leg[mr+1] = (mr - n)*ihsqz2*leg[mr] - (2.0*n - mr + 1)*mr*leg[mr-1];
  }
  for(int m=0; m<2*n; ++m, pre = -pre) leg[m] *= pre;
}


void speedTester(Problem *problem, int maxDeg, Point *r)
{
  double Gnp, Fnp, Enp;
  printf("gnp starting\n");
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      Gnp = calcGnp(problem->e, problem->positions, problem->charges, problem->nCharges, n, p);
      //Gnp = calcI(problem->e, n, p, .67);
    }
  }
  printf("gnp done\n");

  printf("fnp starting\n");  
  for(int n=0; n<=maxDeg; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      //calcI(problem->e, n, p, .67);
      //Fnp = (2*n + 1) * calcEnp(problem->e, problem->positions, n, p) * calcI(problem->e, n, p, problem->positions->x1);
      Fnp = calcFnp(problem->e, r, n, p);
    }
  }
  printf("fnp done\n");
}
