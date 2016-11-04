#include <math.h>
#include <stdio.h>
#include <mpfr.h>
#include <petsc.h>
#include <stdlib.h>
#include <string.h>

#include "ellipsoid/ellipsoid.h"

int testCoordinateTransform()
{
  //Verify that the ellipsoidal to Cartesian transformation is a bijection
  //the tolerance is 10 significant figures
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  
  double relErrorX, relErrorY, relErrorZ;

  //loop over octants
  for(int sx=-1; sx<=1; sx+=2) { //sx = [-1,1]
    for(int sy=-1; sy<=1; sy+=2) { //sy = [-1,1]
      for(int sz=-1; sz<=1; sz+=2) { //sz = [-1,1]
	//loop over brick in cartesian space
	for(double x=4; x<5+tol; x+=1.5) {
	  for(double y=4; y<5+tol; y+=1.5) {
	    for(double z=4; z<5+tol; z+=1.5) {
	      Point temp = { .x1 = sx*x, .x2 = sy*y, .x3 = sz*z };
	      cartesianToEllipsoidal(&e, &temp);
	      ellipsoidToCartesian(&e, &temp);
	      relErrorX = fabs(sx*x - temp.x1)/x;
	      relErrorY = fabs(sy*y - temp.x2)/y;
	      relErrorZ = fabs(sz*z - temp.x3)/z;
	      //check error
	      if( relErrorX > tol ||
		  relErrorY > tol ||
		  relErrorZ > tol) {
		printf("testCoordinateTransform failed\n");
		printf("relative error(x): %16.16f\n", relErrorX);
		printf("relative error(y): %16.16f\n", relErrorY);
		printf("relative error(z): %16.16f\n", relErrorZ);
		return 0;
	      }
	    }
	  }
	}
      }
    }
  }
  return 1;
}

int testLameType()
{

  static char formula1[18] = "K^0_1 L^0_1 M^0_1";
  static char formula2[30] = "K^0_2 K^1_2 L^0_2 M^0_2 N^0_2";
  static char formula3[66] = "K^0_5 K^1_5 K^2_5 L^0_5 L^1_5 L^2_5 M^0_5 M^1_5 M^2_5 N^0_5 N^1_5";

  char **checks = (char**) malloc(sizeof(char*)*3);
  checks[0] = (char*) malloc(sizeof(char)*17);
  checks[1] = (char*) malloc(sizeof(char)*29);
  checks[2] = (char*) malloc(sizeof(char)*65);

  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  int index[3] = { 1, 2, 5 };
  int k;
  int n;
  for(int i=0; i<3; ++i) {
    n = index[i];
    for(int p=0; p < 2*n + 1; ++p) {
      k = p*6;
      checks[i][k]   = getLameTypeT(n, p);
      checks[i][k+1] = '^';
      checks[i][k+2] = getLameTypeTp(n, p) + '0';
      checks[i][k+3] = '_';
      checks[i][k+4] = n + '0';
      if( p != 2*n )
	checks[i][k+5] = ' ';
    }
  }
  //check first string
  for(int i=0; i < 17; ++i) {
    if(checks[0][i] != formula1[i]) {
      return 0;
    }
  }
  //check second string
  for(int i=0; i < 29; ++i) {
    if(checks[1][i] != formula2[i]) {
      return 0;
    }
  }
  //check first string
  for(int i=0; i < 65; ++i) {
    if(checks[2][i] != formula3[i]) {
      return 0;
    }
  }
  free(checks[0]);
  free(checks[1]);
  free(checks[2]);
  free(checks);
  return 1;
}

int testI()
{
  //Using results from Dassios, test the computation of I^p_n for n = 0, 1, 2
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  double a = e.a;
  double b = e.b;
  double c = e.c;
  double h = e.h;
  double k = e.k;
  double l = a;
  double rootSum = sqrt((a*a*a*a - b*b*c*c) + (b*b*b*b - a*a*c*c) + (c*c*c*c - a*a*b*b));
  double DassiosLambda = (1/3.)*(a*a + b*b + c*c) + (1/3.)*rootSum;
  double DassiosLambdaPrime = (1/3.)*(a*a + b*b + c*c) - (1/3.)*rootSum;

  double errors[18];
  int derror[18] = { 1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 8, 9, 10, 11, 12, 13 };
  int ecount = 0;
  //Dassios D1 error
  double lhs = 3 * (DassiosLambda + DassiosLambdaPrime);
  double rhs = 2 * (a*a + b*b + c*c);
  errors[ecount] = fabs(lhs - rhs); ecount++;
  //Dassios D2 error
  lhs = 3 * DassiosLambda * DassiosLambdaPrime;
  rhs = a*a*b*b + a*a*c*c + b*b*c*c;
  errors[ecount] = fabs(lhs - rhs); ecount ++;
  //Dassios D3 error
  lhs = ((-1)*(b*b - c*c)*(DassiosLambda - a*a)) + (k*k*(DassiosLambda - b*b)) + (-h*h*(DassiosLambda - c*c));
  rhs = ((-1)*(b*b - c*c)*(DassiosLambdaPrime - a*a)) + (k*k*(DassiosLambdaPrime - b*b)) + (-h*h*(DassiosLambdaPrime - c*c));
  errors[ecount] = fabs(lhs-0); ecount++;
  errors[ecount ] = fabs(rhs-0); ecount++;
  //Dassios D4 error
  lhs = (-a*a*(b*b - c*c)*(DassiosLambda - a*a)) + (b*b*k*k*(DassiosLambda - b*b)) + (-c*c*h*h*(DassiosLambda - c*c));
  rhs = (-a*a*(b*b - c*c)*(DassiosLambdaPrime - a*a)) + (b*b*k*k*(DassiosLambdaPrime - b*b)) + (-c*c*h*h*(DassiosLambdaPrime - c*c));
  double analytical = h*h*k*k*(b*b - c*c);
  errors[ecount] = fabs(lhs - analytical); ecount++;
  errors[ecount] = fabs(rhs - analytical); ecount++;
  //Dassios D5 error
  double alpha[3] = { a, b, c };
  lhs = 0;
  rhs = 0;
  for(int i=0; i<3; ++i) {
    lhs += alpha[i]*alpha[i]/(alpha[i]*alpha[i] - DassiosLambda);
    rhs += alpha[i]*alpha[i]/(alpha[i]*alpha[i] - DassiosLambdaPrime);
  }
  errors[ecount] = fabs(lhs - 3); ecount++;
  errors[ecount] = fabs(rhs - 3); ecount++;
  //Dassios D6 error
  double anal1 = h*h*k*k*(b*b - c*c);
  double num1  = 3*(b*b - c*c)*(DassiosLambda - a*a)*(DassiosLambdaPrime - a*a);
  double anal2 = -h*h*k*k*(b*b - c*c);
  double num2  = 3*k*k*(DassiosLambda - b*b)*(DassiosLambdaPrime - b*b);
  double anal3 = h*h*k*k*(b*b - c*c);
  double num3  = 3*h*h*(DassiosLambda - c*c)*(DassiosLambdaPrime - c*c);
  errors[ecount] = fabs(num1 - anal1); ecount++;
  errors[ecount] = fabs(num2 - anal2); ecount++;
  errors[ecount] = fabs(num3 - anal3); ecount++;
  
  //Calculate I^p_n
  int nmax = 2;
  int Isize = 0;
  for(int n=0; n < nmax+1; ++n) {
    for(int p=0; p < 2*n+1; ++p) {
      Isize++;
    }
  }
  double *Ival = (double*) malloc(sizeof(double)*Isize);
  int i=0;
  for(int n=0; n < nmax+1; ++n) {
    for(int p=0; p < 2*n+1; ++p) {
      calcI(&e, n, p, l, 1, 1, Ival+i); //Ival[i] = calcI(&e, n, p, l, 1, 1);
      i++;
    }
  }
  
  //Dassios D7: their h_3 = our h and their h_2 = our k
  double analyticalSumTest1 = 1.0/(l*sqrt(l*l - h*h)*sqrt(l*l - k*k));
  double numericalSumTest1 = 0;
  for(int k=1; k<4; ++k) {
    numericalSumTest1 += Ival[k];
  }
  errors[ecount] = fabs(analyticalSumTest1 - numericalSumTest1); ecount++;
  
  //Dassios D8: their alpha1 = our a, their alpha2 = our b, their alpha3 = our c
  double analyticalSumTest2 = Ival[0] - (l*l - a*a)/(l*sqrt(l*l - h*h)*sqrt(l*l - k*k));
  double numericalSumTest2 = a*a*Ival[1] + b*b*Ival[2] + c*c*Ival[3];
  errors[ecount] = fabs(analyticalSumTest2 - numericalSumTest2); ecount++;
  //Dassios D9: their Lambda = our DassiosLambda, their LambdaPrime = our DassiosLambdaPrime
  double analyticalSumTest3 = 1.0/(2.0*(DassiosLambda - a*a + l*l)*l*sqrt(l*l - h*h)*sqrt(l*l - k*k)) - (1.0/2.0)*(Ival[1]/(DassiosLambda - a*a) + Ival[2]/(DassiosLambda - b*b) + Ival[3]/(DassiosLambda - c*c));
  double numericalSumTest3 = Ival[4];
  errors[ecount] = fabs(analyticalSumTest3 - numericalSumTest3); ecount++;
  //Dassios D10: their Lambda = our DassiosLambda, their LambdaPrime = our DassiosLambdaPrime
  double analyticalSumTest4 = 1.0/(2.0*(DassiosLambdaPrime - a*a + l*l)*l*sqrt(l*l - h*h)*sqrt(l*l - k*k)) - (1.0/2.0)*(Ival[1]/(DassiosLambdaPrime - a*a) + Ival[2]/(DassiosLambdaPrime - b*b) + Ival[3]/(DassiosLambdaPrime - c*c));
  double numericalSumTest4 = Ival[5];
  errors[ecount] = fabs(analyticalSumTest4 - numericalSumTest4); ecount++;
  //Dassios D11: their h1*h1 = our (b*b - c*c)
  double analyticalSumTest5 = (1.0/(h*h)) * (Ival[2] - Ival[1]);
  double numericalSumTest5 = Ival[6];
  errors[ecount] = fabs(analyticalSumTest5 - numericalSumTest5); ecount++;
  //Dassios D12: their h1*h1 = our (b*b - c*c)
  double analyticalSumTest6 = (1.0/(k*k)) * (Ival[3] - Ival[1]);
  double numericalSumTest6 = Ival[7];
  errors[ecount] = fabs(analyticalSumTest6 - numericalSumTest6); ecount++;
  //Dassios D13: their h1*h1 = our (b*b - c*c)
  double analyticalSumTest7 = (1.0/(b*b - c*c)) * (Ival[3] - Ival[2]);
  double numericalSumTest7 = Ival[8];
  errors[ecount] = fabs(analyticalSumTest7 - numericalSumTest7); ecount++;
  
  //checking errors
  for(int i=0; i < 18; ++i) {
    if(errors[i] > tol) {
      printf("Dassios D%d failed with error %8.8e\n", derror[i], errors[i]);
      free(Ival);
      return 0;
    }
  }
  
  free(Ival);
  return 1;
}

int testNormalization()
{
  //Check calculated \gamma^p_n against analytic results for n = 0, 1, 2
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, a, b, c);
  //Dassios (B14)
  double firstTerm = (a*a + b*b + c*c)/3.0;
  double secondTerm = sqrt((a*a*a*a - b*b*c*c) + (b*b*b*b - a*a*c*c) + (c*c*c*c - a*a*b*b))/3.0;
  double LambdaD = firstTerm + secondTerm;
  double LambdaDprime = firstTerm - secondTerm;
  //Dassios (B16)-(B20)
  double hx = sqrt(b*b - c*c);
  double hy = e.k;
  double hz = e.h;
  double analytic[9] = { 4*M_PI,
			 4*M_PI/3 * hy*hy*hz*hz,
			 4*M_PI/3 * hx*hx*hz*hz,
			 4*M_PI/3 * hx*hx*hy*hy,
			 -8*M_PI/5 * (LambdaD - LambdaDprime)*(LambdaD - a*a)*(LambdaD - b*b)*(LambdaD - c*c),
			 8*M_PI/5 * (LambdaD - LambdaDprime)*(LambdaDprime - a*a)*(LambdaDprime - b*b)*(LambdaDprime - c*c),
			 4*M_PI/15 * hx*hx*hy*hy*hz*hz*hz*hz,
			 4*M_PI/15 * hx*hx*hy*hy*hy*hy*hz*hz,
			 4*M_PI/15 * hx*hx*hx*hx*hy*hy*hz*hz };
  double value, estimate, error;
  int counter = 0;
  for(int n=0; n<3; ++n) {
    for(int p=0; p<2*n+1; ++p) {
      value = analytic[counter];
      calcNormalization(&e, n, p, &estimate);
      error = estimate - value;
      if(fabs(error) > tol) {
	printf("testnormalization failed iteration %d at (n,p) = (%d,%d) with error %8.8e\n", counter, n, p, fabs(error));
	return 0;
      }
      counter++;
    }
  }
  return 1;

}

void testfunc(mpfr_t *x, mpfr_t *val, FuncInfo4 *ctx)
{
  double a2 = (*ctx).a2;
  double b2 = (*ctx).b2;
  double c2 = (*ctx).c2;
  double botVar = (*ctx).botVar;
  
  mpfr_t temp, s;
  mpfr_init2(temp, mpfr_get_prec(*x));
  mpfr_init2(s   , mpfr_get_prec(*x));
  //change of variables
  //s = (1-x)/x
  mpfr_d_sub(s, 1.0, *x, MPFR_RNDZ);
  mpfr_div(s, s, *x, MPFR_RNDN);

  mpfr_add_d(temp, s, a2, MPFR_RNDN);
  mpfr_add_d(*val, s, b2, MPFR_RNDN);
  mpfr_mul(*val, *val, temp, MPFR_RNDN);
  mpfr_add_d(temp, s, c2, MPFR_RNDN);
  mpfr_mul(*val, *val, temp, MPFR_RNDN);
  mpfr_sqrt(*val, *val, MPFR_RNDN);
  
  mpfr_add_d(temp, s, botVar, MPFR_RNDN);
  mpfr_mul(*val, *val, temp, MPFR_RNDN);
  
  mpfr_d_div(*val, 1.0, *val, MPFR_RNDN);

  //multiply by dt = ds * (1+t)^2
  mpfr_add_d(temp, s, 1.0, MPFR_RNDN);
  mpfr_mul(temp, temp, temp, MPFR_RNDN);
  mpfr_mul(*val, *val, temp, MPFR_RNDN);
  //mpfr_clears(temp, s, NULL);
}

int testSurfaceOperatorEigenvalues()
{
  EllipsoidalSystem e;
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  initEllipsoidalSystem(&e, a, b, c);
  double l = a;
  int n = 1;
  double analytic[3];
  double integral;


  mpfr_t mpfrzero, mpfrone;
  mpfr_inits(mpfrzero, mpfrone, NULL);
  mpfr_set_d(mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(mpfrone, 1.0, MPFR_RNDN);

  //Integrals are from Ritter, need to transform integrals to finite interval
  FuncInfo4 ctx = { &e, a*a, b*b, c*c, a*a };
  integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) testfunc, &e, mpfrzero, mpfrone, 14, &integral, &ctx);
  analytic[0] = (a*b*c * integral - 1.0)/2.0;
  ctx.botVar = b*b;
  integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) testfunc, &e, mpfrzero, mpfrone, 14, &integral, &ctx);
  analytic[1] = (a*b*c * integral - 1.0)/2.0;
  ctx.botVar = c*c;
  integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) testfunc, &e, mpfrzero, mpfrone, 14, &integral, &ctx);
  analytic[2] = (a*b*c * integral - 1.0)/2.0;


  double *ev = (double*) malloc(sizeof(double)*(3));
  for(int p=0; p < 3; ++p)
    calcSurfaceOperatorEigenvalues(&e, n, p, l, 1, 1, ev+p); //ev[p] = calcSurfaceOperatorEigenvalues(&e, n, p, l, 1, 1);
  for(int p=0; p < 3; ++p) {
    if(fabs(ev[p] - analytic[p]) > tol) {
      printf("analytic (%8.8e) doesn't match numerical (%8.8e), error = (%8.8e) for integral %d\n", analytic[p], ev[p], fabs(ev[p] - analytic[p]), p+1);
      return 0;
    }
  }
  free(ev);
  return 1;

}

int compareIntegration()
{

  //int n, p given in calcNormalization
  int n = 1;
  int p = 1;

  int digits = 32;

  //set up ellipsoidal system
  EllipsoidalSystem e;
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  initEllipsoidalSystem(&e, a, b, c);
  
  mpfr_t bound_a, bound_b;  
  mpfr_inits2(4*digits, bound_a, bound_b, NULL);
  
  FuncInfo2 ctx1 = { .e = &e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };

  double integralMPFR, integralMidpoint;

  integrateMPFR((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, &e, e.hp_h, e.hp_k, 14, &integralMPFR, &ctx1);
  integrateMidpoint((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e.hp_h, e.hp_k, 10, &integralMidpoint, &ctx1);

  printf("ze MPFR integral is: %15.15f\n", integralMPFR);
  printf("ze Midpoint integral is: %15.15f\n", integralMidpoint);
  printf("oh wow\n");


  return 0;
}

int main(int argc, char **argv)
{
  mpfr_set_default_prec(4*64);
  /*
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  printf("h: %15.15f\nh2: %15.15f\nk: %15.15f\nk2: %15.15f\n", e.h, e.h2, e.k, e.k2);
  int n = 1;
  int p = 1;
  double l = 1.0;
  
  //printf("calcLame(1, 1, 1.0) = %15.15f\n", calcLame(&e, 1, 1, 1.0));
  
  printf("calcNormalization(%d, %d, %15.15f) = %15.15f\n", n, p, l, calcNormalization(&e, n, p));
  */
  
  if(testCoordinateTransform())
    printf("testCoordinateTransform         passed\n");
  else
    printf("testCoordinateTransform         failed\n");
  if(testLameType())
    printf("testLameType                    passed\n");
  else
    printf("testLameType                    failed\n");
  if(testI())
    printf("testI                           passed\n");
  else
    printf("testI                           failed\n");    
  if(testNormalization())
    printf("testNormalization               passed\n");
  else
    printf("testNormalization               failed\n");
  if(testSurfaceOperatorEigenvalues())
    printf("testSurfaceOperatorEigenvalues  passed\n");
  else
    printf("testSurfaceOperatorEigenvalues  failed\n");

  //compareIntegration();
  
  /* 
  double h = .1;
  mpfr_t x, val;
  EllipsoidalSystem e;
  initEllipsoidalSystem(&e, 3.0, 2.0, 1.0);
  FuncInfo2 ctx1 = { .e = &e, .n = 1, .p = 2, .numeratorType = 0, .denomSign = -1 };
  mpfr_inits2(4*16, x, val, NULL);
  for(int k=0; k < 10; ++k) {
    mpfr_set_d(x, k*h, MPFR_RNDN);
    normFunction1(&x, &val, &ctx1); 
    printf("normfunction1(%8.8f) = %15.15f\n", k*h, mpfr_get_d(val, MPFR_RNDN));
  }
  */
  return 0;
}
