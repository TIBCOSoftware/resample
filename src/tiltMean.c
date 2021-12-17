/*
* Copyright 2021. TIBCO Software Inc.
* This file is subject to the license terms contained
* in the license file that is distributed with this file.
*/

#include <S.h>
#include <machine.h>
static S_LONG debug;
#include "tilt.h"


/* Functions and subroutines in this file:

   Entry points:
     S_tiltMean_solve   tiltMeanSolve                 - given mean, find tau
     S_tiltMean_exp     tiltMean, exponential tilting - given tau, find mean
     S_tiltMean_ml      tiltMean, ML tilting          - given tau, find mean 

   Internal routines:
     solve_mlg_setup       setup for tiltMeanSolve, ML tilting with groups
     copy_dbl              copy one double array to another
     square                square a double
     Qlimits               sum of the max (and min) value of each group of L
     tiltMean_exp          tiltMean, exponential tilting
     tiltMean_ml           tiltMean, ML tilting, no groups
     tiltMean_mlg          tiltMean, ML tilting, groups
     bisect_tiltMean_exp   root-finding for tiltMeanSolve, exp tilting
     bisect_tiltMean_ml    root-finding for tiltMeanSolve, ML tilting
     uniroot_tiltMean_exp  root-finding for tiltMeanSolve, exp tilting
     uniroot_tiltMean_ml   root-finding for tiltMeanSolve, ML tilting
     update_lambda_end_pts find minimum possible lambda for ML tilting, groups
     reset_P               zero out global variables
     update_P              update P and its derivatives
     update_lambda         modified Newton step for updating lambda
     update_tau            Newton step for updating tau
     tiltMean_solve_ml     tiltMeanSolve for ML tilting, groups

   Routines defined in tilt.h
     vectorMin
     vectorMax
     weightedMean
     sumWeightedMeans
     checkLambda           check given lambas give weights that sum to 1
     mlg_setup             setup for tiltMean, ML tilting with groups

   Calling tree (most functions that don't call other functions are omitted)

     solve_mlg_setup
       S_tiltFindLambdas (in tilt.c)
     tiltMean_exp
       sumWeightedMeans
       Qlimits
       weightedMean
     tiltMean_ml
       weightedMean
     tiltMean_mlg
       mlg_setup
     bisect_tiltMean_exp
       tiltMean_exp
     bisect_tiltMean_ml
       tiltMean_ml
     uniroot_tiltMean_exp
       tiltMean_exp
       bisect_tiltMean_exp
     uniroot_tiltMean_ml
       tiltMean_ml
       bisect_tiltMean_ml
     update_lambda_end_pts
     update_lambda
       update_P
     tiltMean_solve_ml
       solve_mlg_setup
       update_lambda_end_pts
       update_P
       update_lambda

     S_tiltMean_solve
       Qlimits
       Mean
       uniroot_tiltMean_exp
       tiltMean_solve_ml
       uniroot_tiltMean_ml
     S_tiltMean_exp
       tiltMean_exp
     S_tiltMean_ml
       tiltMean_mlg
       tiltMean_ml



   ---------------------------------------------------------------------
   The basic relationship in this file is:
      Q = sum( v[i] * L[i] )
   where
      v[i] = c weights[i] * exp(tau * L[i])  (exponential tilting)
      v[i] = c weights[i] / (1 - tau * L[i]) (ML tilting, no groups)
   where
      c    normalizes the output to sum to 1.

   For groups,
      Q = sum_g (sum_i v[gi] L[gi]) = sum of weighted means
   where L[gi] is the i'th observation in group g, g=1 to G

   where, for exponential tilting,
      v[gi] = c[g] * weights[gi] * exp(tau / h[g] * L[gi]))
   and, for ML tilting,
      v[gi] = weights[gi] / (lambda[g] - tau / h[g] * L[gi])

   and where
      h[g] = N_g/N,
      c[g], lambda[g] are normalizing constants so sum(v[g*]) = 1.

   Q is the tilted mean, a function of L, tau, weights, group info
   tiltMean finds Q, given tau
   tiltMeanSolve finds tau given Q


   Note that:
      weights[] corresponds to the prior  weights $v$ in the tech report.
      v[]       corresponds to the output weights $w$ in the tech report.

   Group means are usually subtracted from L before calling this code.

   -- Other quantities --
   xtol = on the scale of tau or lambda
   ytol = on the scale of tilted mean, sum of probs, or bootstrap prob
*/


/* Defined in tilt.c */
void S_tiltFindLambdas(double Ltau[], double weights[],
		       S_LONG *pnGroups, S_LONG groupSizes[],
		       double *pytol, S_LONG *pmaxiter,
		       double lambda[]);


/* Global variables */
static char message[200];

/* Following are for tiltMeanSolve_ml */
static double  P0;
static double *Pi;
static double  Pnorm;
static double  dP0_dtau;
static double *dP0_dlambda;
static double *dPi_dtau;
static double *dPi_dlambda;






/*--------------------------------------------------*/
static void solve_mlg_setup(double tau,
			    double L[],
			    double weights[],
			    S_LONG N,
			    S_LONG nGroups,
			    S_LONG groupSizes[],
			    double ytol,   /* tol for sum of weights*/
			    S_LONG maxIter,
			    double Lmin[],
			    double Lmax[],
			    double lambda[])
{
  /* Return min/max of L for each group, and
     get initial guesses for lambda.

     weights are assumed to exist and to sum to 1 in each group.
  */
  S_LONG g, i, Ng;
  double *Ltau, *Ltau2, taug;

  if(debug>0) printf("Entering solve_mlg_setup\n");
  Ltau = Salloc(N, double);
  Ltau2 = Ltau;  /* save the pointer */

  for(g=0; g < nGroups; g++){
    Ng = groupSizes[g];
    Lmin[g] = vectorMin(Ng, L);
    Lmax[g] = vectorMax(Ng, L);

    i = Ng;
    taug = tau / Ng * (double) N;
    while(i--){
      Ltau[i] = L[i] * taug;
    }

    L       += Ng;        /* Set up for next group */
    Ltau    += Ng;
  }

  /* Get initial guess for lambdas */
  S_tiltFindLambdas(Ltau2, weights, &nGroups, groupSizes,
		    &ytol, &maxIter, lambda);
  return;
}


/*--------------------------------------------------*/
static void copy_dbl(double from[], double to[], S_LONG size)
{
  /* Copy one double array to another. */
  S_LONG i;

  i = size;
  while(i--) to[i] = from[i];
  return;
}


/*--------------------------------------------------*/
static double square(double x)
{
  return(x * x);
}


/*--------------------------------------------------*/
static void Qlimits(double L[],
		    S_LONG nGroups,
		    S_LONG groupSizes[],
		    double *Qmin,
		    double *Qmax)
{
  /* Return the sum of the maximum (and minimum) values of each group of L */
  S_LONG g, Ng;

  if(debug>0) printf("Entering Qlimits\n");
  *Qmin = 0;
  *Qmax = 0;
  for(g=0; g < nGroups; g++){
    Ng = groupSizes[g];
    *Qmin += vectorMin(Ng, L);
    *Qmax += vectorMax(Ng, L);
    L += Ng;
  }
}


/*--------------------------------------------------*/
static double tiltMean_exp(double tau, double *L, S_LONG N, double *weights,
			   S_LONG nGroups, S_LONG *groupSizes )
{
  /* Return exponential tilted mean, or sum of group means.
     Tilted weights are not explicitly computed; instead the
     tilted mean is computed on the fly, using
     numerically-stable calculations, guard against large tau.

     L and weights are vectors of length N, tau a scalar.
     Let v[i] = weights[i] * exp(tau * L[i])
     Return Q = sum(v[i] * L[i]) / sum(v[i])

     For stratified data (nGroups>1),
       Q = sum(Q[g]),
       Q[g] = sum_i( v[gi] * L[gi] ) / sum_g(v[gi]),
       g = 1 to G (number of groups)
       v[gi] = weights[gi] * exp(tau * L[gi] / h[g]),
     and where
       h[g] = N[g]/N, N[g] = size of group i

     If weights==NULL then no weights were supplied.
  */
  S_LONG g, i, Ng;
  S_LONG *pNg;
  double sumW, sumWL, taug, v, Lmax, result = 0.0;
  double Qmin, Qmax;

  if(debug>1) printf("Entering tiltMean_exp\n");
  if(tau == 0)
    return(sumWeightedMeans(N, L, weights, nGroups, groupSizes));

  if(is_inf(&tau, S_MODE_DOUBLE)){
    Qlimits(L, nGroups, groupSizes, &Qmin, &Qmax);
    return( (tau < 0) ? Qmin : Qmax );
  }

  /* In the following calculations, subtracting Lmax adds
     numerical stability.   Algebraically it cancels out.*/
  for(g=0; g<nGroups; g++){
    pNg = groupSizes + g; /* pointer to group size for this group */
    Ng = *pNg;
    taug = tau * (double) N / (double) Ng;

    sumW  = 0.0;
    sumWL = 0.0;

    if(tau < 0)
      Lmax = vectorMin(Ng, L);
    else
      Lmax = vectorMax(Ng, L);

    if(weights){
      for(i = 0; i < Ng; i++) {
	v = exp(taug * (L[i]-Lmax)) * weights[i];
	sumW  += v;
	sumWL += v * L[i];
      }
      weights += Ng; /* Set up for next group. */
    }
    else {                 /* weights is a NULL pointer, no weights */
      for(i = 0; i < Ng; i++) {
	v = exp(taug * (L[i]-Lmax));
	sumW  += v;
	sumWL += v * L[i];
      }
    }
    result += sumWL / sumW ;
    L += Ng;     /* Set up for next group */
  }
  return(result);
}


/*--------------------------------------------------*/
static double tiltMean_ml(double tau,
			  double L[],
			  S_LONG N,
			  double weights[])
{
  /* Compute tilted mean, ML tilting, no groups.
     Use a numerically-stable method in case tau is close to
     a value that makes the denominator blow up.

     L and weights are vectors of length N.
     Let v[i] = 1 / (1 - tau * L[i])
     Return Q = sum(weights[i] * v[i] * L[i]) / sum(weights[i] * v[i])

     If weights==NULL then no weights were supplied.
  */
  S_LONG i;
  double sumW, sumWL, tmax, ratio, denominator, minimumDenominator;

  if(debug>1) printf("Entering tiltMean_ml\n");
  if(tau == 0)
    return(weightedMean(N, L, weights));

  if(tau < 0)
    tmax = tau * vectorMin(N, L);
  else
    tmax = tau * vectorMax(N, L);

  if(tmax >= 1.0){
    na_set3(&denominator, S_MODE_DOUBLE, Is_NaN);
    PROBLEM
      "tau out of range to make weights positive; setting results to NA"
      WARNING(NULL_ENTRY);
    return(denominator);
  }

  /* minimumDenominator is for numerical stability; will cancel out. */
  minimumDenominator = 1.0 - tmax;
  sumW  = 0.0;
  sumWL = 0.0;

  if(weights){
    for(i = 0; i < N; i++) {
      denominator = 1.0 - tau * L[i];
      ratio = minimumDenominator/denominator * weights[i];
      sumW  += ratio;
      sumWL += ratio * L[i];
    }
  }
  else {                 /* weights is a NULL pointer, no weights */
    for(i = 0; i < N; i++) {
      denominator = 1.0 - tau * L[i];
      ratio = minimumDenominator/denominator;
      sumW  += ratio;
      sumWL += ratio * L[i];
    }
  }
  return(sumWL / sumW);
}


/*--------------------------------------------------*/
static double tiltMean_mlg(double tau,        /* scalar*/
			   double L[],        /* vector[N]*/
			   S_LONG N,            /* number of obs, all groups*/
			   double weights[],  /* vector[N]*/
			   double Ltau[],     /* vector[N], work space*/
			   S_LONG nGroups,      /* number of groups*/
			   S_LONG groupSizes[], /* number in each group*/
			   S_LONG haveLambda,   /* TRUE if lambda supplied*/
			   double ytol,       /* tol for sum of weights*/
			   S_LONG maxiter,
			   double lambda[])   /* vector[nGroups]*/
{
  /* Computed tilted mean, ML tilting, groups.

     Let v[gi] = weights[gi] / (lambda[g] - tau / h[g] * L[gi])),
     Return Q = sum_g (sum(v[g*] * L[g*])) = sum of weighted group means.

     The lambda[g] are computed and returned.

     It is assumed that the data is organized contiguously in groups, so that
     the first groupSizes[0] elements of L and weights belong to a group, the
     next groupSizes[1] belong to another, and so on.

     weights are assumed to exist and to sum to 1 within each group.
  */
  double result, Lambda;
  S_LONG g, i, Ng;

  if(debug>1) printf("Entering tiltMean_mlg\n");

  /* Compute Ltau = L * tau/h_i, and lambdas (if necessary) */
  mlg_setup(tau, L, weights, N, nGroups, groupSizes, haveLambda,
	    ytol, maxiter, lambda, Ltau);

  result = 0.;
  for(g = 0; g < nGroups; g++){
    Lambda = lambda[g];
    i = Ng = groupSizes[g];
    while(i--) result += L[i] * weights[i] / (Lambda - Ltau[i]);
    L += Ng, Ltau += Ng, weights += Ng;  /* set up for next group */
  }
  return(result);
}


/*--------------------------------------------------*/
static double bisect_tiltMean_exp(double a,
				  double b,
				  double Fa,
				  double Fb,
				  double Q,
				  double xtol,
				  double ytol,
				  S_LONG *iter,
				  double *L,
				  S_LONG lenL,
				  double *weights,
				  S_LONG nGroups,
				  S_LONG groupSizes[])
{
  /* Solve tiltMean(tau) = Q.

     a and b are endpoints of interval
     Fa = tiltMean(a) - Q, Fb = tiltMean(b) - Q
     Fa and Fb must have opposite signs

     Note that this is actually a bracketed secant method,
     modified to move the next point slightly toward the midpoint.
  */
  double c, Fc, lambda;

  if(debug>0) printf("Entering bisect_tiltMean_exp\n");
  if(Fa*(Fb/fabs(Fb)) > 0) {
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
	    "internal error, roots don't bracket answer in bisect"
	    WARNING(NULL_ENTRY);
  }
  xtol = xtol/2.;

  do {
    /* c = (a + b)/2.0; */  /* bisect */
    lambda = Fb / (Fb-Fa);
    lambda = .5 + .99 * (lambda - .5);
    c = lambda * a + (1.0 - lambda) * b;

    Fc = tiltMean_exp( c, L, lenL, weights, nGroups, groupSizes) - Q;
    (*iter)++;

    if(fabs(Fc) < ytol)
      return c;

    if(Fa*(Fc/fabs(Fc)) <  0.0) {
      b   = c;
      Fb = Fc;
    }
    else {
      a  = c;
      Fa = Fc;
    }
  }
  while(fabs(b-a) > xtol*(1.0 + fabs(c)));

  return c;
}


/*--------------------------------------------------*/
static double bisect_tiltMean_ml(double a, double b, double Fa, double Fb,
				 double Q, double xtol, double ytol, S_LONG *iter,
				 double *L, S_LONG lenL, double *weights)
{
  /* Solve tiltMean(tau) = Q.

     a and b are endpoints of interval
     Fa = tiltMean(a) - Q, Fb = tiltMean(b) - Q
     Fa and Fb must have opposite signs

     Note that this is actually a bracketed secant method,
     modified to move the next point slightly toward the midpoint.
  */
  double c, Fc, lambda;

  if(debug>0) printf("Entering bisect_tiltMean_ml\n");
  if(Fa*(Fb/fabs(Fb)) > 0) {
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
	    "internal error, roots don't bracket answer in bisect"
	    WARNING(NULL_ENTRY);
  }
  xtol = xtol/2.;

  do {
    lambda = Fb / (Fb-Fa);
    lambda = .5 + .99 * (lambda - .5);
    c = lambda * a + (1.0 - lambda) * b;

    Fc = tiltMean_ml(c, L, lenL, weights) - Q;
    (*iter)++;

    if(fabs(Fc) < ytol)
      return c;

    if(Fa*(Fc/fabs(Fc)) <  0.0) {
      b   = c;
      Fb = Fc;
    }
    else {
      a  = c;
      Fa = Fc;
    }
  }
  while(fabs(b-a) > xtol*(1.0 + fabs(c)));

  return c;
}


/*--------------------------------------------------*/
static double uniroot_tiltMean_exp(double a, double F0, double Q,
				   double xtol, double ytol, S_LONG *iter,
				   double *L, S_LONG lenL, double *weights,
				   S_LONG nGroups, S_LONG groupSizes[])
{
  /* Solve tiltMean(tau) = Q. */
  S_LONG    MaxIter;
  double Fa, b, Fb, d, Fd, S0, x;

  if(debug>0) printf("Entering uniroot_tiltMean_exp\n");
  MaxIter = *iter;
  *iter   = 0;

  F0 = F0 - Q;
  if(fabs(F0) < ytol) return(0.0);
  S0 = F0/fabs(F0);  /* This is positive when answer should be negative. */

  if(!a)
    a = .001;
  /* If a is on the wrong side of 0, reverse it. */
  if(S0 * a > 0)
    a = -a;

  Fa = tiltMean_exp(a, L, lenL, weights, nGroups, groupSizes) - Q;
  *iter = 1;
  if(fabs(Fa) < ytol) return(a);

  if(S0*Fa < 0) {
    /* If 0 and a bracket the root, use bisection. */
    x  = bisect_tiltMean_exp(a, 0.0, Fa, F0, Q, xtol, ytol, iter, L, lenL,
			     weights, nGroups, groupSizes);
    return x;
  }

  if(fabs(Fa) > fabs(F0)) {
    /* printf("Internal error in uniroot_tiltMean_exp\n"); */
    PROBLEM
	    "internal error in uniroot_tiltMean_exp"
	    WARNING(NULL_ENTRY);
  }

  /* At this point a & 0 don't bracket, and a is closer. */
  b = a;
  Fb = Fa;
  a = 0.0;
  Fa = F0;

  while(*iter <= MaxIter) {
    /* Use secant method to get new estimates;
       if estimates bracket the root then use bisection. */
    d  = b - a;
    Fd = Fb - Fa;
    if((Fb-Fa)/(b-a) <= 0){
      sprintf(message,
	      "nonmonotonicity (a, b, mean(a), mean(b)) = (%f, %f, %f, %f)",
	      a, b, Fa+Q, Fb+Q);
      PROBLEM
	message
	WARNING(NULL_ENTRY);
      na_set3(&b, S_MODE_DOUBLE, Is_NaN);
      return(b);
    }
    x  = b - 1.01 * Fb*d/Fd;
    a  = b;
    Fa = Fb;
    b  = x;
    Fb = tiltMean_exp(b, L, lenL, weights, nGroups, groupSizes) - Q;
    ++(*iter);
    if(fabs(Fb) < ytol) return(b);
    if(Fb * Fa < 0){ /* bracketed */
      x = bisect_tiltMean_exp(a, b, Fa, Fb, Q, xtol, ytol, iter, L, lenL,
			      weights, nGroups, groupSizes);
      return x;
    }
  }

  if(fabs(d) <= xtol*(1.+fabs(b))) return(b);

  PROBLEM
      "iteration limit reached in uniroot_tiltmean_exp"
      WARNING(NULL_ENTRY);
  return(b);
}


/*--------------------------------------------------*/
static double uniroot_tiltMean_ml(double a, double F0,
				  double Q, double xtol, double ytol,
				  double tauneg, double taupos, S_LONG *iter,
				  double *L, S_LONG lenL, double *weights)
{
  /* Solve tiltMean(tau) = Q. */
  S_LONG    MaxIter;
  double Fa, b, Fb, d, Fd, S0, x;

  if(debug>0) printf("Entering uniroot_tiltMean_ml\n");
  MaxIter = *iter;
  *iter   = 0;

  F0 = F0 - Q;
  if(fabs(F0) < ytol) return(0.0);
  S0 = F0/fabs(F0);  /* This is positive when answer should be negative. */

  if(!a)
    a = .001;
  if(S0 * a > 0)
    a = -a;
  if(a <= tauneg) a = tauneg/2.0;
  if(a >= taupos) a = taupos/2.0;

  Fa = tiltMean_ml(a, L, lenL, weights) - Q;
  *iter  = 1;
  if(fabs(Fa) < ytol) return(a);

  if(S0*Fa < 0) {
    /* If 0 and a bracket the root, use bisection. */
    x  = bisect_tiltMean_ml(a, 0.0, Fa, F0, Q, xtol, ytol, iter, L, lenL,
			    weights);
    return x;
  }

  if(fabs(Fa) > fabs(F0)) {
    PROBLEM
      "internal error in uniroot_tiltMean_ml"
      WARNING(NULL_ENTRY);
  }

  /* At this point a & 0 don't bracket, and a is closer. */
  b = a;
  Fb = Fa;
  a = 0.0;
  Fa = F0;

  while(*iter <= MaxIter) {
    /* Use secant method to get new estimates;
       if estimates bracket the root then use bisection. */
    d  = b - a;
    Fd = Fb - Fa;
    if((Fb-Fa)/(b-a) <= 0){
      sprintf(message,
	      "nonmonotonicity (a, b, mean(a), mean(b)) = (%f, %f, %f, %f)",
	      a, b, Fa+Q, Fb+Q);
      PROBLEM
	message
	WARNING(NULL_ENTRY);
      na_set3(&b, S_MODE_DOUBLE, Is_NaN);
      return(b);
    }
    x  = b - 1.00 * Fb*d/Fd;
    a  = b;
    Fa = Fb;
    b  = x;
    if(b <= tauneg) b = (a+tauneg)/2.0;
    if(b >= taupos) b = (a+taupos)/2.0;

    Fb = tiltMean_ml(b, L, lenL, weights) - Q;
    ++(*iter);
    if(fabs(Fb) < ytol) return(b);
    if(Fb * Fa < 0){ /* bracketed */
      x = bisect_tiltMean_ml(a, b, Fa, Fb, Q, xtol, ytol, iter, L, lenL,
			     weights);
      return x;
    }
  }

  if(fabs(d) <= xtol*(1.+fabs(b))) return(b);

  PROBLEM
    "iteration limit reached in uniroot_tiltmean_ml"
    WARNING(NULL_ENTRY);
  return(b);
}


/*--------------------------------------------------*/
static void update_lambda_end_pts(double tau, S_LONG N, S_LONG nGroups,
			       S_LONG groupSizes[], double Lmin[], double Lmax[],
			       double endPts[])
{
  /* Calculate the minimum value possible for lambda, given tau.
     The acceptable interval is of width 1.
  */
  S_LONG i;
  double *Lextr;

  Lextr = (tau < 0) ? Lmin : Lmax;
  i = nGroups;
  while(i--)
    endPts[i] = tau * (double) N / groupSizes[i] * Lextr[i];
  return;
}


/*--------------------------------------------------*/
static void reset_P(S_LONG nGroups)
{
  /* Zero out all global variables set by update_P */
  S_LONG i;

  P0 = 0;
  Pnorm = 0;
  dP0_dtau = 0;
  i = nGroups;
  while(i--){
    Pi[i] = 0;
    dP0_dlambda[i] = 0;
    dPi_dlambda[i] = 0;
    dPi_dtau[i] = 0;
  }
  return;
}


/*--------------------------------------------------*/
static void update_P(double tau,
		    double L[],
		    double lambda[],
		    double weights[],
		    double Q,
		    S_LONG N,
		    S_LONG nGroups,
		    S_LONG groupSizes[])
{
  /* Update P and its derivatives. */
  S_LONG gi, g, i;
  double h, denom, Lval;
  double w, wL, W, WL, WLL, sumWLL;

  reset_P(nGroups); /* zero out all the quantities */
  gi = 0;           /* This points to i'th observation in group g. */
  for(g = 0; g < nGroups; g++){
    h = groupSizes[g] / (double) N;
    sumWLL = 0.;
    for(i = 0; i < groupSizes[g]; i++){
      Lval = L[gi];
      denom = lambda[g] - tau * Lval / h;
      w = weights[gi] / denom;
      wL  = w  * Lval;
      W   = w  / denom;
      WL  = W  * Lval;
      WLL = WL * Lval;
      P0 += wL;
      Pi[g] += w;
      sumWLL += WLL;
      dPi_dtau[g] += WL;
      dPi_dlambda[g] -= W;
      gi++;
    }
    Pi[g] -= 1.;
    dP0_dtau += (sumWLL / h);
    dP0_dlambda[g] = -dPi_dtau[g];
    dPi_dtau[g] /= h;
    Pnorm += square(Pi[g]);
  }
  P0 -= Q;
  Pnorm = sqrt(Pnorm + P0 * P0);
}


/*--------------------------------------------------*/
static void update_lambda(double tau, double lambda[], double L[],
			 double weights[], double Q,
			 S_LONG N, S_LONG nGroups, S_LONG groupSizes[],
			 double endPts[], double endPtsNew[],
			 double lambdaNew[])
{
  /* modfied Newton step for updating lambda */
  S_LONG g;
  double *lambdaPN;
  lambdaPN = Salloc(nGroups, double);

  /* Locate lambda the same place relative to the new endpoints
     as the old lambda was relative to the old endpoints */
  for(g = 0; g < nGroups; g++)
    lambdaPN[g] = lambda[g] - endPts[g] + endPtsNew[g];

  /* Newton step */
  update_P(tau, L, lambdaPN, weights, Q, N, nGroups, groupSizes);
  for(g = 0; g < nGroups; g++){
    lambdaNew[g] = lambdaPN[g] - Pi[g] / dPi_dlambda[g];
    /* Adjust any out-of-bounds lambda by bringing it in halfway to
       the boundary from where it was previously. */
    if(lambdaNew[g] <= endPtsNew[g])
      lambdaNew[g] = .5*(endPtsNew[g] + lambdaPN[g]);
    else if(lambdaNew[g] >= endPtsNew[g] + 1.)
      lambdaNew[g] = .5*(endPtsNew[g] + lambdaPN[g] + 1.);
  }
}


/*--------------------------------------------------*/
static double update_tau(double tau,
			double lambda[],
			double Q,
			S_LONG N,
			double weights[],
			S_LONG nGroups,
			S_LONG groupSizes[])
{
  /* Newton step for the first component of solving P = 0 uses the
     first row of the inverse of the Jacobian of P.
  */
  double result, Jinv11;
  S_LONG g;

  /* The [1,1] entry of the inverse; it gets multiplied by everything */
  Jinv11 = dP0_dtau;
  for(g = 0; g < nGroups; g++)
    Jinv11 -= dP0_dlambda[g] / dPi_dlambda[g] * dPi_dtau[g];

  Jinv11 = 1 / Jinv11;
  /* Here begins the multiplication of the first row of the Jacobian
     inverse by P. */
  result = P0;
  for(g = 0; g < nGroups; g++)
    result -= dP0_dlambda[g] * dPi_dlambda[g] * Pi[g];

  /* Newton update: new = old - Jinv*P */
  return(tau - result * Jinv11);
}


/*--------------------------------------------------*/
static void tiltMean_solve_mlg(double tau,
			       double Q,
			       double xtol,   /* tol for difference in tau*/
			       double ytol,   /* tol for sum of weights, and ?*/
			       S_LONG *iter,
			       double L[],
			       S_LONG N,
			       double weights[],
			       S_LONG nGroups,
			       S_LONG groupSizes[],
			       double lambda[],
			       double *tau_out)
{
  /* Solve tilted mean = Q, for  maximum likelihood tilting and groups.
   Results are written to lambda and tau_out.

   This is a modified Newton's method for finding the root of the
   system of nGroups + 1 equations

     P(tau, lambda[i]) =

     {P0(tau, lambda[i]),
      Pi(tau, lambda[i])} =

     {sum_gi(v[gi](tau, lambda[g])*L[gi]) - Q  (weighted mean)
      sum_i(v[gi](tau, lambda[g]) - 1},        (weights sum to 1 in each group)

   for g = 1 to nGroups. The algorithm does a Newton update for tau
   only (using the first row of the inverse of the Jacobian for the
   full system), and given that tau, a Newton update for each of the
   nGroups equations for lambda_i (these are univariate equations, if
   tau is considered to be fixed, so Newton updates can be performed
   indpendently).

   An exception to the previous algorithm occurs if after the lambda
   updates the results are worse than the previous iteration.  In that
   case, tau is adjusted to be the average of the two previous updates
   and a Newton update for lambda is performed using the new tau.

   The two workhorse routines are update_tau and update_lambda.
  */

  double *lambdaNew, tauNew;
  double *Lmin, *Lmax, *endPts, *endPtsNew, curPnorm, tauDiff;
  S_LONG maxIter = *iter, i;
  int moreLambda;

  if(debug>0) printf("Entering tiltMean_solve_ml\n");

  /* Allocate global variables for evaluating P and its derivatives */
  Pi          = Salloc(nGroups, double);
  dP0_dlambda = Salloc(nGroups, double);
  dPi_dtau    = Salloc(nGroups, double);
  dPi_dlambda = Salloc(nGroups, double);

  /* Additional local storage */
  lambdaNew = Salloc(nGroups, double);
  Lmin      = Salloc(nGroups, double);
  Lmax      = Salloc(nGroups, double);
  endPts    = Salloc(nGroups, double);
  endPtsNew = Salloc(nGroups, double);

  /* Subtract off the group means from L, get max/min of each group,
     and get initial guesses for lambda. */
  solve_mlg_setup(tau, L, weights, N, nGroups, groupSizes, ytol,
		  maxIter, Lmin, Lmax, lambda);
  /* The acceptable interval for lambda is of width 1 with these left
     end points. */
  update_lambda_end_pts(tau, N, nGroups, groupSizes, Lmin, Lmax, endPts);
  /* Evaluate P and its derivatives for tau, lambda.  This also
     computes global variable Pnorm, the norm of P, which we want to
     be zero. */
  update_P(tau, L, lambda, weights, Q, N, nGroups, groupSizes);
  curPnorm = Pnorm;
  tauDiff = xtol + 1;

  /* Main loop */
  for(i = 1; (curPnorm > ytol || tauDiff > xtol) && i <= maxIter; i++){
    tauNew = update_tau(tau, lambda, Q, N, weights, nGroups, groupSizes);
    tauDiff = fabs((tau - tauNew)/tau);
    moreLambda = 1;
    while(moreLambda){
      update_lambda_end_pts(tauNew, N, nGroups, groupSizes, Lmin, Lmax,
			    endPtsNew);
      update_lambda(tauNew, lambda, L, weights, Q, N, nGroups, groupSizes,
		   endPts, endPtsNew, lambdaNew);
      update_P(tauNew, L, lambdaNew, weights, Q, N, nGroups, groupSizes);

      if(Pnorm < curPnorm || tauDiff < xtol){
	/* Ready for next iteration of tau */
	/* In latter case, we're probably stuck in a no-progress loop */
	moreLambda = 0;
	tau = tauNew;
	curPnorm = Pnorm;
 	copy_dbl(lambdaNew, lambda, nGroups);
	copy_dbl(endPtsNew, endPts, nGroups);
      }
      else{
	/* adjust tau, then go back and use old lambdas to find new ones */
	tauNew = .5 * (tau + tauNew);
	tauDiff = fabs((tau - tauNew)/tau);
      }
    }
  }

  if(i > maxIter){
    PROBLEM
      "maximum iterations reached before convergence"
      WARNING(NULL_ENTRY);
  }
  *iter = i;
  *tau_out = tau;
  return;
}


/*--------------------------------------------------*/
void S_tiltMean_solve(S_LONG *Nobs,
		      S_LONG *Ntau,
		      S_LONG *Nw,
		      double L[],
		      double weights[],
		      S_LONG *pnGroups,
		      S_LONG groupSizes[],
		      double *xTol,
		      double *yTol,
		      S_LONG *MaxEval,
		      S_LONG *wantML,
		      S_LONG *pDebug,
		      double expTau[],
		      double mlTau[],
		      double lambda[])
/* Find tau values (and possibly lambda) do that tilted mean is Q.

   L is a vector of length *Nobs
   expTau and mlTau are vectors of length *Ntau
   lambda is a vector of length (*Ntau) * (*pnGroups)

   On input, initial values for tau are in expTau and Q is in
   mlTau.
*/
{
  S_LONG N = *Nobs;
  S_LONG K = *Ntau;
  S_LONG nGroups = *pnGroups;
  double xtol  = *xTol;
  double ytol = *yTol;

  S_LONG k, iter, imax;
  double mean, tauneg, taupos, Qmin, Qmax;
  double tau, Q, alpha;


  if(*Nw == 0) weights = NULL;
  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltMean_solve\n");

  /* Find the maximum and minimum possible tilted means */
  Qlimits(L, nGroups, groupSizes, &Qmin, &Qmax);

  if(*wantML && nGroups <= 1){
    /* find domain for ML tilting, parameterization for no groups */
    tauneg = (Qmin < 0.) ? 1./Qmin : -log((double)SINGLE_XMAX)/2.;
    taupos = (Qmax > 0.) ? 1./Qmax :  log((double)SINGLE_XMAX)/2.;

    if((taupos-tauneg) < 2.*xtol) {
      PROBLEM
	"ML range < 2*xtol"
	WARNING(NULL_ENTRY);
    }
  }

  imax = 0;
  mean = sumWeightedMeans(N, L, weights, nGroups, groupSizes);
  for(k = 0; k < K; k++) {
    Q = mlTau[k];
    if(Q < Qmin - SINGLE_EPS || Q > Qmax + SINGLE_EPS){
      PROBLEM
	"q out of range; setting results to NA"
	WARNING(NULL_ENTRY);
      na_set3(mlTau+k, S_MODE_DOUBLE, Is_NaN);
      na_set3(expTau+k, S_MODE_DOUBLE, Is_NaN);
      continue;
    }
    else if(Q < Qmin + SINGLE_EPS){
      PROBLEM
	"q is the minimum possible value; setting tau accordingly"
	WARNING(NULL_ENTRY);
      inf_set(expTau+k, S_MODE_DOUBLE, -1);
      if(nGroups > 1)
	inf_set(mlTau+k, S_MODE_DOUBLE, -1);
      else
	mlTau[k] = tauneg;
      continue;
    }
    else if(Q > Qmax - SINGLE_EPS){
      PROBLEM
	"q is the maximum possible value; setting tau accordingly"
	WARNING(NULL_ENTRY);
      inf_set(expTau+k, S_MODE_DOUBLE, 1);
      if(nGroups > 1)
	inf_set(mlTau+k, S_MODE_DOUBLE, 1);
      else
	mlTau[k] = taupos;
      continue;
    }
    tau  = expTau[k];
    iter = *MaxEval;
    tau  = uniroot_tiltMean_exp(tau, mean, Q, xtol, ytol, &iter, L, N, weights,
				nGroups, groupSizes);
    expTau[k] = tau;
    if(iter > imax) imax = iter;

    if(!(*wantML))
      continue;
    iter  = *MaxEval;
    if(nGroups > 1){
      tiltMean_solve_mlg(tau, Q, xtol, ytol, &iter, L, N, weights,
			 nGroups, groupSizes, lambda, mlTau + k);
      lambda += nGroups;  /* set up for the next tau */
    }
    else{
      tau  = tau/(1. + tau*Q);
      alpha = (Q - Qmin) / (Qmax - Qmin);
      if(tau > alpha * Qmax)
	tau = alpha * Qmax;
      if(tau < (1.0-alpha) * Qmin)
	tau = (1.0-alpha) * Qmin;
      mlTau[k] = uniroot_tiltMean_ml(tau, mean, Q, xtol, ytol, tauneg, taupos,
				     &iter, L, N, weights);
    }
    if(iter > imax) imax = iter;
  }
  *MaxEval = imax;

  return;
}


/*--------------------------------------------------*/
void S_tiltMean_exp(S_LONG *Nobs, S_LONG *Ntau, S_LONG *Nw,
		    double *L, double *weights,
		    S_LONG *pnGroups, S_LONG groupSizes[],
		    S_LONG *pDebug,
		    double *Tau)
{
  /* Calculated tilted means.

     values in Tau on input; probabilities in Tau on output
  */
  S_LONG N, K, nGroups;
  S_LONG k;

  N = *Nobs;
  K = *Ntau;
  nGroups = *pnGroups;
  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltMean_exp\n");

  if(*Nw == 0) weights = NULL;

  for(k = 0; k < K; k++) Tau[k] = tiltMean_exp(Tau[k], L, N, weights,
					       nGroups, groupSizes);

  return;
}


/*--------------------------------------------------*/
void S_tiltMean_ml(S_LONG *Nobs,
		   S_LONG *Ntau,
		   S_LONG *Nw,
		   double L[],        /* length N */
		   double weights[],  /* length N, or NULL */
		   S_LONG *pnGroups,
		   S_LONG groupSizes[],
		   S_LONG *phaveLambda,
		   double *pytol,
		   S_LONG *pmaxiter,
		   S_LONG *pDebug,
		   double Tau[],
		   double lambda[])
{
  /* Calculated tilted means.

     values in Tau on input; probabilities in Tau on output
  */
  S_LONG N = *Nobs;
  S_LONG K = *Ntau;
  S_LONG nGroups = *pnGroups;
  S_LONG haveLambda = *phaveLambda;
  double ytol = *pytol;
  S_LONG maxiter = *pmaxiter;
  double *Ltau;
  S_LONG k;

  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltMean_ml\n");
  if(*Nw == 0) weights = NULL;
  Ltau = Salloc(N, double); /* Work space for tiltMean_mlg*/

  if(nGroups > 1){
    for(k = 0; k < K; k++){
      Tau[k] = 
	tiltMean_mlg(Tau[k], L, N, weights, Ltau, nGroups, groupSizes,
		     haveLambda, ytol, maxiter, lambda);
      lambda += nGroups;  /* set up for next group */
    }
  }
  else
    for(k = 0; k < K; k++)
      Tau[k] = tiltMean_ml(  Tau[k], L, N, weights);    

  return;
}
