/*
* Copyright 2021. TIBCO Software Inc.
* This file is subject to the license terms contained
* in the license file that is distributed with this file.
*/

#include <S.h>
#include <machine.h>

/* Global variables */
static S_LONG debug;
static double xTol, fTol;
static double s, epsilon, Y0, slope, quadratic, *sumWeights;
static double K1_s, K2_s, K3_s;
static char message[200];

/* Functions and subroutines in this file,
   in the order they appear, and what they call

   Calling trees:

   S_revSaddlepointPSolve_exp   Solve the exponential tilting parameter tau,
                                from the reverse saddlepoint cdf estimate
				probability {alpha}, with given the sum of
     (call:)                    total tilting premeters: {s}; called from S.
     initializeWeights
     initializeCenter
     uniroot_saddle         Solve Saddle2((s - tau), ...) = S_qnorm(1-alpha)

   initializeCenter         Assign global static variables values for quadratic
     (call:)                interpolation
     findKs                 Find K'(s), K''(s), K''(s)
     Saddle2                Saddle for the 2nd tilting premeter {t},
                            with given {tau} or {s},  different from the
                            code "Saddle" of "saddlepointP.c".

   Saddle2                  Returns qnorm(p(t)), where p(t) is the probability
                            of which exceeding the 2nd tilted mean, when given
                            sampling from the 1st tilted distribution.
     findK                  find K(s-t)

   uniroot_saddle           Solve a :  Saddle2(a, ...) = S_qnorm(1-alpha)
                            with given s, where	a is an initial guess of
     (call:)                (s - tau).
     Saddle2
     bisect_saddle          An bracketed secant method to solve t within (a, b)
                            in case there "Fa * Fb < 0"; where
			    Fa = Saddle2(a, ...) - S_qnorm(1-alpha),
			    Fb = Saddle2(b, ...) - S_qnorm(1-alpha).

   bisect_saddle
     (call:)
     Saddle2

   Rename original "Mean" and "Saddle" to be "preMean" and "Saddle2" to avoid
   that functions in this file may have names conflicted with corresponding
   functions in "saddle.c", "tiltBoot.c" and "tiltMean.c".

   ---------------------------------------------------------------------
   static variables and where they are set:
       s, epsilon, Y0, slope, quadratic, sumWeights, K1_s, K2_s
   All but {s} are set in "initializeCenter", which is called by
   entry point: "S_revSaddlepointPSolve_exp"; {s}, the sum of the 1st and 2nd
   exponential tilting premeters, will be called the earliest directly
   by "S_revSaddlepointPSolve_exp".

   ---------------------------------------------------------------------
   The basic relationship in this file is:

      K(tau) = n * log( sum(weights[i] * exp(tau * L[i])) / sum(weights))
      { K(t) | given tau } = K(t + tau) - K(tau)
      Q = { K'(t) | given tau }
        = n * sum( w[i] * L[i]) / sum( w[i] )
   where
      w[i] = weights[i] * exp((t + tau) * L[i])
      If (tau + t) == 0, then
      Q = the original weighted mean, and { K(t) | given tau } = - K(tau)
   also
      p(tau | given t) = 1 - pnorm(v + log(u/v)/v)                   if t != 0
                       = 1 - pnorm(K'''(t=0|tau)/K''(t=0|tau)^1.5)   if t == 0
   where
      u = t * sqrt({ K''(t) | given tau })
        = t * sqrt(K''( s ))
      v = sign(t) *
          sqrt(2 * (t * { K'(t) | given tau } - { K(t) | given tau }))
	= sign(t) * sqrt(2 * (t * K'(s) - K(s) + K(tau)))
   for
      s = t + tau

   If t is near 0, use quadratic interpolation
   Note that:
      { K'(t) | tau }   = n * (tilted mean with tilting s )
      { K''(t) | tau }  = n * (tilted variance with tilting s )
      { K'''(t) | tau } = n * (tilted centered third moment with tilting s )

   L should normally be centered by removing its (tau tilting weighted) mean
   prior to calling routines in this file.

   If there are multiple groups, then in place of tau use tau*N/N[g];
   however this is implemented by dividing the L's by (N[g]/N)
   prior to calling this function.  Hence this function effectively
   does calculations for sums of group sums, with the result that
   the S code does calculations for sums of group means.
*/



/******* Prototypes for functions in this file */
static void initializeWeights(double L[], double weights[],
			     S_LONG nGroups, S_LONG groupSizes[]);
static void initializeCenter(double L[], double weights[],
			     S_LONG nGroups, S_LONG groupSizes[]);
static double uniroot_saddle(double a, double Zalpha,
			     double L[], double weights[],
			     S_LONG nGroups, S_LONG groupSizes[], 
			     double xTol, double fTol, S_LONG *iter);
static double Saddle2(double t, double L[], double weights[],
		      S_LONG nGroups, S_LONG groupSizes[]);
static void findK(double tau, double L[], double weights[], 
		  S_LONG nGroups, S_LONG groupSizes[], double *K);
static void findKs(double L[], double weights[], 
		   S_LONG nGroups, S_LONG groupSizes[]);
static double bisect_saddle(double a, double b,
			    double Fa, double Fb,
			    double Zalpha,
			    double L[], double weights[],
			    S_LONG nGroups, S_LONG groupSizes[], 
			    double xTol, double fTol, S_LONG *iter);



/******* Begin function definitions *******/


void S_revSaddlepointPSolve_exp(S_LONG *pNtau, S_LONG *pNweights,
				double *L, double *weights, double *alpha, 
				double *pS, S_LONG *pnGroups, S_LONG groupSizes[], 
				double *Tau0, double *pXtol, double *pFtol,
				S_LONG *pMaxIter, S_LONG *pDebug,
				S_LONG *pAccIter, double *Tau)
{
  S_LONG M, i, iter, imax, MaxIter, nGroups;

  /* count = 0; */
  M = *pNtau;
  nGroups = *pnGroups;
  MaxIter = *pMaxIter;
  s = *pS;

  if(*pNweights == 0) weights = NULL;

  xTol   = *pXtol;
  fTol   = *pFtol;
  debug = *pDebug;

  imax  = 0;

  if(debug > 0) printf("Entering S_revSaddlepointPSolve_exp\n");

  initializeWeights(L, weights, nGroups, groupSizes);
  initializeCenter (L, weights, nGroups, groupSizes);

  for (i = 0; i < M; i++) {
    if(alpha[i] <= 0.0 || alpha[i] >= 1.0) {
      na_set3(Tau + i, S_MODE_DOUBLE, Is_NaN);
      continue;
    }

    iter = MaxIter;
    Tau[i] = s - uniroot_saddle(( s - Tau0[i] ), 
				S_qnorm( 1 - alpha[i] ), L, weights, 
				nGroups, groupSizes, xTol, fTol, &iter);
    if(iter > imax) imax = iter;
  }

  *pAccIter = imax;
  if(debug > 0) printf("Exiting S_revSaddlepointPSolve_exp\n");
}


static double Saddle2(double t, double *L, double *weights,
		      S_LONG nGroups, S_LONG groupSizes[])
{
  double u, v, sum0, K, tau;

  if(debug > 2) printf("Saddle2n: t = %g, epsilon = %g\n", t, epsilon);
  if(fabs(t) < epsilon) /* quadratic interpolation */
    return( Y0 + t * (slope + t * quadratic)); /* new */
/*      return( Y0 + t * (slope - t * quadratic)); */

  sum0 = 0.0;

  tau = s - t;
  findK(tau, L, weights, nGroups, groupSizes, &K);

  if(debug > 2)
    printf("Saddle2n: K = %g, K1_s = %g, K2_s = %g\n", K, K1_s, K2_s);
  u = t * sqrt(K2_s);
  v = t/fabs(t) * sqrt(2 * (t * K1_s - K));
  if(debug > 2) printf("Saddle2n: u = %g, v = %g, return = %g\n",
		       u, v, v + log(u/v)/v);
  return(v + log(u/v)/v);
}

/***********************************************************************
  findK

  Evaluate K(tau).
***********************************************************************/
static void findK(double tau, double L[], double weights[], 
		  S_LONG nGroups, S_LONG groupSizes[], double *K)
{
  S_LONG i, j, N;
  double w, sum0, sum1, mean;

  K[0] = 0.0;

  /* For grouped data, the value is just the sum for each group */
  for(j = 0; j < nGroups; j++){

    sum0 = sum1 = 0.0;
    N = groupSizes[j];

    if(weights){
      for(i=0; i<N; i++){
	w = weights[i] * exp(tau * L[i]);
	sum0 += w;
	sum1 += w * L[i];
      }
      weights += N; /* set up for next group */
    }
    else {
      for(i=0; i<N; i++){
	w = exp(tau * L[i]);
	sum0 += w;
	sum1 += w * L[i];
      }
    }
    mean = sum1 / sum0;
    K[0] += N * log(sumWeights[j] / sum0);
    L += N; /* set up for next group */
  }
  return;
}



static void findKs(double L[], double weights[], 
		   S_LONG nGroups, S_LONG groupSizes[])
{
  S_LONG i, j, N;
  double sum0, sum1, sum2, sum3, temp, temp2, mean;

  K1_s = K2_s = K3_s = 0.0;

  /* For grouped data, the values are just the sum of the values for each
     group. */
  if(s == 0.0){
    for(j=0; j<nGroups; j++){
      sum0 = sum1 = sum2 = sum3 = 0.0;
      N = groupSizes[j];
      if(weights){
	for(i=0; i<N; i++){
	  sum0 += weights[i];
	  sum1 += weights[i] * L[i];
	}
      }
      else {
	sum0 = N;
	for(i=0; i<N; i++)
	  sum1 += L[i];
      }
      mean = sum1 / sum0;
      K1_s += N * mean;
      if(weights){
	for(i = 0; i < N; i++){
	  temp = (L[i] - mean);
	  temp2 = weights[i] * temp * temp;
	  sum2 += temp2;
	  sum3 += temp2 * temp;
	}
	weights += N; /* set up for next group */
      }
      else{
	for(i = 0; i < N; i++){
	  temp = (L[i] - mean);
	  temp2 = temp * temp;
	  sum2 += temp2;
	}
      }
      K2_s += N * sum2 / sum0;
      K3_s += N * sum3 / sum0;

      L += N; /* set up for next group */
    }
  }
  else{
    for(j=0; j<nGroups; j++){
      sum0 = sum1 = sum2 = sum3 = 0.0;
      N = groupSizes[j];
      if(weights){
	for(i=0; i<N; i++){
	  sum0 += weights[i] * exp(s * L[i]);
	  sum1 += weights[i] * exp(s * L[i]) * L[i];
	}
      }
      else {
	sum0 = N;
	for(i=0; i<N; i++){
	  sum0 += exp(s * L[i]);
	  sum1 += exp(s * L[i]) * L[i];
	  sum1 += L[i];
	}
      }
      mean = sum1 / sum0;
      K1_s += N * mean;
      if(weights){
	for(i = 0; i < N; i++){
	  temp = L[i] - mean;
	  temp2 = weights[i] * temp * temp * exp(s * L[i]);
	  sum2 += temp2;
	  sum3 += temp2 * temp;
	}
	weights += N; /* set up for next group */
      }
      else{
	for(i = 0; i < N; i++){
	  temp = L[i] - mean;
	  temp2 = temp * temp * exp(s * L[i]);
	  sum2 += temp2;
	}
      }
      K2_s += N * sum2 / sum0;
      K3_s += N * sum3 / sum0;

      L += N; /* set up for next group */
    }
  }

  return;
}

/***********************************************************************
  initializeWeights

  Compute global variable sumWeights: sum of weights for each group
***********************************************************************/
static void initializeWeights(double L[], double weights[], 
			      S_LONG nGroups, S_LONG groupSizes[])
{
  S_LONG i, j;

  sumWeights = Salloc(nGroups, double);  /* global variable set here */

#ifdef WIN64 //S7_Change_LONG
  if(debug > 1) printf("Entering initializeWeights: s = %g, nGroups = %lld\n", s,nGroups);
#else
  if(debug > 1) printf("Entering initializeWeights: s = %g, nGroups = %ld\n", s,nGroups);
#endif //WIN64
  if(s == 0.0){
    if(weights){
      for(i=0; i < nGroups; i++){
	for(j=0; j < groupSizes[i]; j++){
	  sumWeights[i] += weights[j];
	}
	weights += groupSizes[i]; /* set up for next group */
      }
    }
    else{
      for(i=0; i < nGroups; i++){
	sumWeights[i] = (double)(groupSizes[i]);
      }
    }
  }
  else{
    if(weights){
      for(i=0; i < nGroups; i++){
#ifdef WIN64 //S7_Change_LONG
	if(debug > 3) printf("groupSizes[i] = %lld\n", groupSizes[i]);
#else
    if(debug > 3) printf("groupSizes[i] = %ld\n", groupSizes[i]);
#endif //WIN64
	for(j=0; j < groupSizes[i]; j++){
	  sumWeights[i] += weights[j] * exp(s * L[j]);;
	}
	if(debug > 3) printf("sumWeights[i] = %g\n", sumWeights[i]);
	weights += groupSizes[i]; /* set up for next group */
	L += groupSizes[i]; 
      }
    }
    else{
      for(i=0; i < nGroups; i++){
	for(j=0; j < groupSizes[i]; j++){
	  sumWeights[i] += exp(s * L[j]);;
	}
	L += groupSizes[i]; /* set up for next group */
      }
    }
  }
  if(debug > 1) printf("Exiting initializeWeights\n");
  return;
}

/***********************************************************************
  initializeCenter

  Define quantities used for quadratic interpolation for "t = s-tau" near 0:
     epsilon, Y0, slope, quadratic, K1_s, K2_s.
***********************************************************************/
static void initializeCenter(double *L, double *weights,
			     S_LONG nGroups, S_LONG groupSizes[])
{
  double sigma, temp, y1, y2;

  findKs(L, weights, nGroups, groupSizes);
  sigma = sqrt(K2_s);  /* sigma of a sample sum */
  Y0 = K3_s / (6 * K2_s * sigma);
  /* Saddle interpolates if fabs(tau) <= epsilon; prevent that now */
  epsilon = 0.0;
  temp = 0.003 / sigma;
  y1 = Saddle2( -temp, L, weights, nGroups, groupSizes) / 2.0;
  y2 = Saddle2(  temp, L, weights, nGroups, groupSizes) / 2.0;
  /* divide by 2 in previous to simplify calculations below */
  /* Define domain, slope, and quadratic term for interpolation */
  epsilon = temp;
  slope = (y2 - y1) / ( epsilon);
  quadratic = (y1 + y2 - Y0) / (epsilon * epsilon);
  return;
}


static double uniroot_saddle(double a, double Zalpha,
			     double *L, double *weights,
			     S_LONG nGroups, S_LONG groupSizes[], 
			     double xTol, double fTol, S_LONG *iter)

     /* Solve a, such that Saddle2(a, ...) = Zalpha with given s, where
	a is an initial guess of (s - tau).
     */
{
  S_LONG    MaxIter;
  double F0, Fa, b, Fb, d, Fd, S0, x, eps;

  if(debug > 2) printf("Entering uniroot_saddle\n");
  F0 = Y0 - Zalpha;  /* Zalpha is treated as Ya, if a is the solution. */
  if(fabs(F0) < fTol) return(0.0); /* i.e. a = 0 => tau = s */

  S0 = (F0 > 0.0) ? 1.0 : -1.0;  /* sign of F0 */

  if(debug > 3)
    printf("Y0 = %g, Zalpha = %g, S0 = %g, F0 = %g\n", Y0, Zalpha, S0, F0);
  if(fabs(a) < epsilon) {
  /* If a near 0, change it to one a point used for quadratic interp */
    a = (S0 > 0) ? -epsilon : epsilon;
  }
  else if(S0 * a > 0){
  /* i.e. a is on the wrong side of 0, reverse it. */
    a = -a;
  }

  MaxIter = *iter;

  Fa = Saddle2(a, L, weights, nGroups, groupSizes) - Zalpha;
  *iter = 1;

  if(fabs(Fa) < fTol)
    return(a);

  if(S0*Fa < 0) {
  /* If 0 and a bracket the root, use bisection. */
    x = bisect_saddle(a, 0.0, Fa, F0, Zalpha, L, weights, 
		      nGroups, groupSizes, xTol, fTol, iter);
    return (x);
  }

  /* At this point a & 0 don't bracket, and a is closer. */
  b = a;
  Fb = Fa;
  a = 0.0;
  Fa = F0;
  if(debug > 3)
    printf("before loop: a = %g, b = %g, Fa = %g, Fb = %g\n", a, b, Fa, Fb);
  
  xTol = DOUBLE_EPS + xTol;
  eps = xTol*(1. + fabs(b))/2.;

  while(*iter <= MaxIter) {
    /* Use secant method to get new estimates;
       if estimates bracket the root then use bisection. */
    if(debug > 4) printf("secant %f %f %f %f\n", a, b, Fa, Fb);
    d  = b - a;
    Fd = Fb - Fa;
    if((Fb-Fa)/(b-a) <= 0){
      /* printf("WARNING:  nonmonotonicity (a, b, qnorm(p(a)), qnorm(p(b)))
	 %f %f %f %f\n", a, b, Fa+Zalpha, Fb+Zalpha);
      */
      sprintf(message,
	      "nonmonotonicity (a, b, qnorm(p(a)), qnorm(p(b))) = (%f, %f, %f, %f)",
		a, b, Fa+Zalpha, Fb+Zalpha);
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
    if(debug > 4) printf("uniroot, before Saddle2: a = %g, b = %g\n", a, b);
    Fb = Saddle2(b, L, weights, nGroups, groupSizes) - Zalpha;
    if(debug > 4) printf("Fa = %g, Fb = %g\n", Fa, Fb);
    *iter += 1;
    if(fabs(Fb) < fTol)
      return(b);
    if(Fb * Fa < 0){ /* bracketed */
      x = bisect_saddle(a, b, Fa, Fb, Zalpha, L, weights,
			nGroups, groupSizes, xTol, fTol, iter);
      return(x);
    }
  }

  if(fabs(d) <= eps || fabs(Fb) <= fTol) 
    return(b);

  /* printf("WARNING: iteration limit %ld reached in uniroot_saddle\n",
     MaxIter);
  */

#ifdef WIN64 //S7_Change_LONG
  sprintf(message,
	  "Maxmum iteration limit %lld reached in uniroot_saddle", MaxIter);
#else
  sprintf(message,
	  "Maxmum iteration limit %ld reached in uniroot_saddle", MaxIter);
#endif //WIN64
  PROBLEM
    message
    WARNING(NULL_ENTRY);
  return (b);
}


static double bisect_saddle(double a, double b,
			    double Fa, double Fb,
			    double Zalpha,
			    double *L, double *weights,
			    S_LONG nGroups, S_LONG groupSizes[], 
			    double xTol, double fTol, S_LONG *iter)

     /* Solve Saddle2(t, ..., tau=s-t) = Zalpha.
	a and b are endpoints of interval
	Fa = Saddle2(a, ...) - Zalpha, Fb = Saddle2(b, ...) - Zalpha
	Fa and Fb must have opposite signs
	Note that this is actually a bracketed secant method,
	modified to move the next point slightly toward the midpoint.
     */

{
  double d, c, Fc, lambda;

  if(Fa*(Fb/fabs(Fb)) > 0){
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
	"internal error, roots don't bracket answer in bisect"
	WARNING(NULL_ENTRY);
  }
  xTol = xTol/2.;

  while (1) {
    if(debug > 4) printf("bisect %f %f %f %f\n", a, b, Fa, Fb);
    if(fabs(Fa) < fabs(Fb)) {
      c  =  a;
      a  =  b;
      b  =  c;
      c  = Fa;
      Fa = Fb;
      Fb = c;
    }

    if(fabs(Fb) < fTol) 
      return(b);

    d = b - a;

    if(fabs(d) < xTol * (1. + fabs(b))) 
      return(b);

    /*    c = (a + b)/2.;*/  /* bisect */
    lambda = Fb / (Fb-Fa);
    lambda = .5 + .99 * (lambda - .5);
    c = lambda * a + (1.0 - lambda) * b;

    Fc = Saddle2( c, L, weights, nGroups, groupSizes) - Zalpha;
    ++(*iter);

    if(fabs(Fc) < fTol) 
      return(c);

    if(Fa*(Fc/fabs(Fc)) <  0.) {
      b  = c;
      Fb = Fc;
    }
    else {
      a  = c;
      Fa = Fc;
    }
    if((Fb-Fa)/(b-a) < 0){
      /* printf("WARNING:  nonmonotonicity (a, b, qnorm(p(a)), qnorm(p(b)))
	 %g %g %g %g\n", a, b, Fa+Zalpha, Fb+Zalpha);
      */
      sprintf(message,
	"nonmonotonicity (a, b, qnorm(p(a)), qnorm(p(b))) = (%g, %g, %g, %g)",
	      a, b, Fa+Zalpha, Fb+Zalpha);
      PROBLEM
	message
	WARNING(NULL_ENTRY);
    }
  }
}





