/*
* Copyright 2021. TIBCO Software Inc.
* This file is subject to the license terms contained
* in the license file that is distributed with this file.
*/

#include <S.h>

/* Functions for computing maximum likelihood tilting weights using
   stratified data.

   Functions and subroutines in this file:

   Entry points from S:
     S_tiltFindLambdas    Find lambdas for all groups.

   Callable from other S routines:
     S_findLambda         Find lambdas for a single group.

   Internal routines (some of these are duplicated in tilt.h):
     square             square a double
     vectorMax          maximum of a vector
     printVectorDouble
     printVectorLong
     ml                 basic function of lambda: sum of weights = 1
     dml                derivative of ml with respect to lambda


   In this file, if weights is a null pointer, then take weights to
   be equal to 1/n, where n is the group size.

   -- Other quantities --
   xtol = on the scale of tau or lambda
   ytol = on the scale of tilted mean, sum of probs, or bootstrap prob
*/

/* Global variable (within this file)*/
static S_LONG debug;



static double square(double x)
{
  return(x*x);
}


static double vectorMax(S_LONG N, double *x){
  S_LONG i;
  double max;

  max = x[N-1];
  for(i = 0; i < N; i++)
    if(x[i] > max) max = x[i];
  return(max);
}

static void printVectorDouble(char *name, S_LONG N, double *x){
  S_LONG i;

  printf("%s", name);
  for(i=0; i<N; i++) printf(" %g", x[i]);
  printf("\n");
  return;
}


static void printVectorLong(char *name, S_LONG N, S_LONG *x){
  S_LONG i;

  printf("%s", name);
#ifdef WIN64 //S7_change_LONG
  for(i=0; i<N; i++) printf(" %lld", x[i]);
#else
  for(i=0; i<N; i++) printf(" %ld", x[i]);
#endif //WIN64
  printf("\n");
  return;
}


static double ml(double lambda,    /* scalar*/
		 double Ltau[],    /* vector[n], values of L*tau*/
		 double weights[], /* vector[n], should sum to 1*/
		 S_LONG   n)
{
  /* Sum of ML tilted weights:
       sum_i(weights[i]/(lambda - Ltau[i]))
  */
  S_LONG i = n;
  double result = 0.0;

  while(i--) result += weights[i]/(lambda - Ltau[i]);
  return(result);
}


static double dml(double lambda,
		  double Ltau[],
		  double weights[],
		  S_LONG   n)
{
  /* The derivative of ml (w/respect to lambda). */
  S_LONG i = n;
  double result = 0.0;

  while(i--) result -= weights[i]/square(lambda - Ltau[i]);
  return(result);
}



void S_findLambda(double Ltau[],     /* vector[n], values of L*tau*/
		  double weights[],  /* vector[n], should sum to 1 by group*/
		  S_LONG   n,          /* number of observations*/
		  double ytol,       /* tol for sum of weights*/
		  S_LONG   maxiter,
		  double *Lambda)    /* output, scalar*/
{
  /* Find lambda such that 
       sum_i weights[i]/(lambda - Ltau[i]) = 1

     Use Newton's method.
     Stop if sum of weights is within ytol of 1.
  */
  double lambdaMin, lambdaMax, lambdaNew, lambda;
  double f;
  S_LONG k;
  
  /* The acceptable range for lambda is an interval of length 1 with 
     the following left endpoint.  Anything smaller and weights are negative */
  lambdaMin = vectorMax(n, Ltau);
  lambdaMax = lambdaMin + 1.0;
  if(debug) printf("In S_findLambda, lambdaMin = %10.7f\n", lambdaMin);
  if(debug > 4) {
#ifdef WIN64 //S7_Change_LONG
    printf("n = %lld, ytol = %g, maxiter = %lld\n", n, ytol, maxiter);
#else
    printf("n = %ld, ytol = %g, maxiter = %ld\n", n, ytol, maxiter);
#endif //WIN64
    printVectorDouble("  Ltau", n, Ltau);
    printVectorDouble("  weights", n, weights);
  }

  /* With L centered, if tau=0 then lambda=1, at the maximum of the interval.
     As tau increases, lambda increases, initially at O(tau^2),
     but asympotically is just above the minimum.
     Initial guess:  */
  lambda    = sqrt(1.0 + square(lambdaMin));

  /* main loop */
  for(k = 1; k <= maxiter; k++){
    /* Solve f = ml(lambda, ...) = 1. */
    f = ml(lambda, Ltau, weights, n);
    if(debug > 1) printf("  lambda = %10.7f, f = %10.7f\n", lambda, f);

    if(fabs(f - 1.0) <= ytol) break;

    /* Newton update, based on  g = 1/f - 1 and its derivatives */
    lambdaNew = lambda - (1.0/f - 1.0)/
      (-dml(lambda, Ltau, weights, n)/square(f));

    if(debug && lambdaNew <= lambdaMin)
      printf("\ttoo low, lambdaNew = %g\n", lambdaNew);
    if(debug && lambdaNew > lambdaMax)
      printf("\ttoo high, lambdaNew = %g\n", lambdaNew);

    /* If result is out of bounds, make new guess halfway to the boundary */
    if(lambdaNew <= lambdaMin)
      lambdaNew = .5*(lambdaMin + lambda);
    /* g = 1/f - 1 is increasing and concave down, so Newton never
       overshoots to the right
       else if(lambdaNew > lambdaMax) lambdaNew = .5*(lambdaMax + lambda);
    */

    /* Set up for next iteration */
    lambda = lambdaNew;
  }
#ifdef WIN64 //S7_Change_LONG
  if(debug) printf("S_findLambda, %lld iterations\n", k);
#else
  if(debug) printf("S_findLambda, %ld iterations\n", k);
#endif //WIN64

  if(k > maxiter){
    PROBLEM
      "maximum iterations reached before convergence when finding lambda"
      WARNING(NULL_ENTRY);
  }

  *Lambda = lambda;

  return;
}
			 

void S_tiltFindLambdas(double Ltau[],
		       double weights[], 
		       S_LONG *pnGroups,
		       S_LONG groupSizes[],
		       double *pytol,       /* tol for sum of weights*/
		       S_LONG *pmaxiter,
		       double lambda[])
{
  /* Main entry point.
     Find lambda[g], g = 1..nGroups such that 
       sum_i weights[gi]/(lambda[g] - Ltau[gi]) = 1
     where sum_i weights[gi] = 1 for each g.
  */
  S_LONG nGroups = *pnGroups;
  S_LONG g, nG;

  for(g = 0; g < nGroups; g++){
    nG = groupSizes[g];
    S_findLambda(Ltau, weights, nG, *pytol, *pmaxiter, 
		 lambda + g);
    Ltau += nG;
    weights += nG;
  }
}
