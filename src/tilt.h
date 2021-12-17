/*
* Copyright 2021. TIBCO Software Inc.
* This file is subject to the license terms contained
* in the license file that is distributed with this file.
*/

/* This file contains routines used in other files:
     tiltBoot.c
     tiltMean.c

   Routines in this file:
     vectorMin          Minimum of a vector
     vectorMax          Maximum of a vector
     printVectorDouble  For debugging, print name and vector
     printVectorLong    For debugging, print name and vector
     weightedMean       weighted mean of a vector
     sumWeightedMeans   sum of weighted group means
     checkLambda
     mlg_setup


   Calling tree:
     mlg_setup
       checkLambda
       S_findLambda (in tilt.c)


*/


/* Declaration for routine defined in tilt.c */
void S_findLambda(double Ltau[],     /* vector[n], values of L*tau*/
		  double weights[],  /* vector[n], should sum to 1 by group*/
		  S_LONG   n,          /* number of observations*/
		  double ytol,       /* tol for sum of weights*/
		  S_LONG   maxiter,
		  double *Lambda)    /* output, scalar*/
     ;



static double vectorMin(S_LONG N, double *x){
  S_LONG i;
  double min;

  min = x[N-1];
  for(i = 0; i < N; i++)
    if(x[i] < min) min = x[i];
  return(min);
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
#ifdef WIN64 //S7_Change_LONG
  for(i=0; i<N; i++) printf(" %lld", x[i]);
#else
  for(i=0; i<N; i++) printf(" %ld", x[i]);
#endif //WIN64
  printf("\n");
  return;
}


static double weightedMean(S_LONG N,
			   double *L,       /* vector[n]*/
			   double *weights) /* vector[n], or NULL pointer*/
{
  /* Return the weighted mean of L
     If weights is a NULL pointer, then return the simple mean.
  */
  S_LONG i = N;
  double sumW = 0.0, sumWL = 0.0;

  if(weights){
    while(i--){
      sumW += weights[i];
      sumWL += weights[i] * L[i];
    }
    return(sumWL / sumW);
  }
  else {
    while(i--)
      sumWL += L[i];
    return(sumWL / (double)N);
  }
}



static double sumWeightedMeans(S_LONG N, double *L, double *weights,
			S_LONG nGroups, S_LONG *groupSizes)
{
  /* Return the sum of weighted grouped means of L.
     If weights is a NULL pointer, then return the simple mean.
  */
  S_LONG g, i, Ng;
  double sumW, sumWL, result=0.0;

  for(g=0; g<nGroups; g++){

    i = Ng = groupSizes[g];
    sumW = sumWL = 0.0;

    if(weights){
      while(i--){
	sumW += weights[i];
	sumWL += weights[i] * L[i];
      }
      result += sumWL / sumW;
      weights += Ng;    /* Set up for next group */
    }
    else {
      while(i--) sumWL += L[i];
      result += sumWL / (double) Ng;
    }
    L += Ng;     /* Set up for next group */
  }
  return(result);
}




static void checkLambda(double Ltau[],
			double weights[],
			double lambda[],
			S_LONG   nGroups,
			S_LONG   groupSizes[],
			double ytol)
{
  /* Check that given lambdas make ML weights sum to 1, in each group.
     weights should be a null pointer, or sum to 1 in each group. */
  S_LONG g, i, nG;
  double sumw, Lambda;

  for(g = 0; g < nGroups; g++){
    i = nG = groupSizes[g];
    sumw = 0.0;
    Lambda = lambda[g];
    if(weights)
      while(i--) sumw += weights[i] / (Lambda - Ltau[i]);
    else {
      while(i--) sumw += 1.0 / (Lambda - Ltau[i]);
      sumw /= nG;
    }
    if(fabs(sumw - 1.) > ytol)
      PROBLEM
	"given lambda does not produce weights that sum to 1"
	WARNING(NULL_ENTRY);
    Ltau += nG, weights += nG;  /* set up for next group */
  }
}



static void mlg_setup(double tau,
		      double L[],
		      double weights[],
		      S_LONG N,
		      S_LONG nGroups,
		      S_LONG groupSizes[],
		      S_LONG haveLambda,
		      double ytol,
		      S_LONG maxIter,
		      double lambda[],
		      double Ltau[])
{
  /* Calculate L * tau / h_i and store in vector Ltau,
     calculate lambda[g], for g = 1..nGroups.

     weights should sum to 1 in each group,

     L usually has weighted group means zero.
  */
  S_LONG g, i, Ng;
  double *Ltau2, taug;

  if(debug) printf("Entering mlg_setup, tau=%g\n", tau);
  if(debug > 3) printVectorDouble("  L", N, L);
  if(debug > 3) printVectorDouble("  weights in mlg_setup(1)", N, weights);
  if(debug > 6) printVectorLong("  groupSizes", nGroups, groupSizes);

  Ltau2 = Ltau;  /* save the pointer */

  for(g=0; g < nGroups; g++){
    Ng = groupSizes[g];

    taug = tau * N / (double) Ng;
    i = Ng;
    while(i--) Ltau[i] = L[i] * taug;

    L       += Ng;        /* Set up for next group */
    Ltau    += Ng;
  }

  /* Compute lambdas, if necessary.  Otherwise check that the given
     lambdas give weights that sum to 1. */
  if(haveLambda)
    checkLambda(Ltau2, weights, lambda, nGroups, groupSizes, ytol);
  else {
    Ltau = Ltau2;
    for(g = 0; g < nGroups; g++){
      Ng = groupSizes[g];
      S_findLambda(Ltau, weights, Ng, ytol, maxIter, lambda + g);
      Ltau += Ng;
      weights += Ng;
    }
  }
  if(debug > 3) printVectorDouble("  weights in mlg_setup(2)", N, weights-N);
  if(debug) printf("Exiting mlg_setup\n");
  return;
}
