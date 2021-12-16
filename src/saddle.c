#include <S.h>
#include <machine.h>

/* Saddlepoint calculations for a sum, mean,
   sum of group sums, or sum of group means. */


/***********************************************************************

   The entry points from S-PLUS are:
     S_saddlepointP
     S_saddlepointD
     S_saddlepointPSolve


   Notation:  p(tau) is the saddlepoint cdf estimate
              K is the cgf (cumulant generating function)


   Functions and subroutines in this file, in the order they appear,
     (what they call)    (and what they do):
   
   SaddleP               Returns qnorm(p(tau))
     findK
   SaddleD               Returns d(tau), the saddlepoint density estimate
     findK
   findK                 Evaluate K and its first 3 derivatives at tau
   findK0                Evaluate K'' and K''' at 0
   S_saddlepointP        Evaluate p(tau) for vector tau, called from S
     initializeWeights
     rescaleL
     initializeCenter
     SaddleP
   S_saddlepointD        Evaluate d(tau) for vector tau, called from S
     initializeWeights
     rescaleL
     findK0
     convolutionVariance
     SaddleD
   S_saddlepointPSolve   Solve p(tau) = alpha, vector alpha, called from S
     initializeWeights
     rescaleL
     initializeCenter
     uniroot_saddle
   initializeWeights     assign global vector, sum of the weights by group
   rescaleL              rescale L (in multiple groups case, for sample mean)
   initializeCenter      assign global variables for quadratic interp. near 0
     convolutionVariance
     SaddleP
   convolutionVariance   Return the var of the normal obs for convolution
   uniroot_saddle
     SaddleP
     bracket_saddle
   bracket_saddle
     SaddleP

***********************************************************************/

static char message[200];

/***********************************************************************
   Prototypes:
***********************************************************************/
static double SaddleP(double tau, double L[], double weights[], S_LONG size,
		      S_LONG nGroups, S_LONG groupSizes[], double Work[],
		      double Knvar);
static double SaddleD(double tau, double L[], double weights[], S_LONG size,
		      S_LONG nGroups, S_LONG groupSizes[], double Work[],
		      double Knvar);
static void findK(double tau, double L[], double weights[], S_LONG size,
		  S_LONG nGroups, S_LONG groupSizes[], double Work[],
		  double Knvar, double *K);
static void findK0(double L[], double weights[], S_LONG size,
		   S_LONG nGroups, S_LONG groupSizes[], double Work[],
		   double *K0);
static void initializeWeights(S_LONG N, double weights[], S_LONG nGroups,
			      S_LONG groupSizes[]);
static void rescaleL(S_LONG N, double L[], S_LONG nGroups, S_LONG groupSizes[]);
static void initializeCenter(S_LONG N, double L[], double weights[], S_LONG size,
			     S_LONG nGroups, S_LONG groupSizes[],
			     double conv_factor, double Work[], double *Knvar);
static double convolutionVariance(double *K0, double conv_factor,
				  S_LONG size, S_LONG N);
static double uniroot_saddle(double a, double Zalpha, S_LONG n, double L[],
			     double weights[], S_LONG size, S_LONG nGroups,
			     S_LONG groupSizes[], double tol, double ytol,
			     S_LONG *iter, double Work[], double Knvar);
static double bracket_saddle(double a, double b, double Fa, double Fb,
			     double Zalpha, S_LONG n, double L[],
			     double weights[], S_LONG size, S_LONG nGroups,
			     S_LONG groupSizes[], double tol, double ytol,
			     S_LONG *iter, double Work[], double Knvar);



/***********************************************************************
   Global variables:
   are set by one of two functions:
     initializeWeights:
       sets: sumWeights  (allocates space, initializes vector)
       called by:  all entry points
                   (S_saddlepointPSolve, S_saddlepointP, S_saddlepointD)
     initializeCenter:
       sets: epsilon, Y0, slope, quadratic
       called by:  (S_saddlepointPSolve, S_saddlepointP)

 ***********************************************************************/
static double epsilon, Y0, slope, quadratic, *sumWeights;


/***********************************************************************

   The basic relationship in this file is:
      K(tau) = size * K1(tau)    where K1 is the cgf for a single observation
             = size * log( sum(weights[i] * exp(tau * L[i])) / sum(weights))
      Q = K'(tau)
        = size * sum( w[i] * L[i]) / sum( w[i] )
   where
      w[i] = weights[i] * exp(tau * L[i])

   K is the cumulant generating function for the sum of `size'
   observations drawn from L.  

   We are interested in the density or cdf of the mean (if mean=T) or sum.
   We work with the cgf for the sum, which is simpler.
   This makes it easier to relate to the original observations;
   this tau corresponds to tilting individual obs by exp(tau x).

   The density estimate is

      d(tau) = size * exp(K(tau) - tau * K'(tau)) / sqrt(2pi * K''(tau))

   The p-coordinate is given by Barndorff-Nielsen 1986 Biometrika
   approximation, which is similar to the Lugannani-Rice estimate.
   (Kolassa, p. 91 calls this r*.  BN 1986 has r* as something else.)

      p(tau) = pnorm(v + log(u/v)/v)        if tau != 0
             = pnorm(K'''(0)/K''(0)^1.5)    if tau == 0
   where

      u = tau sqrt(K''(tau))
      v = sign(tau) * sqrt(2 * (tau * K'(tau) - K(tau)))

   If tau is near 0, use quadratic interpolation
   Note that:

      K'(tau)   = n * (tilted mean) = tilted mean for a sample sum
      K''(tau)  = n * (tilted variance) = tilted variance for sum
      K'''(tau) = n * (tilted centered third moment)



   For grouped data, sum of group sums ("calculateMean" is false):
   K is defined as above for each group g (with n replaced by
   the size of the group and the sums covering only the members of the group)
   to give K_g, and then K = sum_g K_g.  You still get the same
   relationships between derivatives of K and tilted mean, variance, etc.

   For grouped data, sum of group means ("calculateMean" is true):
   the input tau is mapped nominally mapped to tau / (n[g]/n),
   on the scale of each individual sample.  In practice we implement
   by dividing each group by (n[g]/n).  Then proceed as for sum of sums.

   It is assumed that the data is organized contiguously in groups, so that
   the first groupSizes[0] elements of L and weights belong to a group, the
   next groupSizes[1] belong to another, and so on.

   L should normally be centered by removing its (weighted) mean
   prior to calling routines in this file.  This has been tested
   with uncentered L and gives good results.

   The convolution factor `conv_factor' refers to a single, normally
   distributed observation with mean 0 and variance equal to
   Knvar = conv_factor * var(L). This normal distribution has cgf equal to
   Kn = .5 * (*Knvar) * tau * tau, which is added to K to get the final cgf.
   conv_factor = 0 corresponds to no convolution.

 ***********************************************************************/



/***********************************************************************
  SaddleP
***********************************************************************/
static double SaddleP(double tau, double L[], double weights[], S_LONG size,
		      S_LONG nGroups, S_LONG groupSizes[], double Work[],
		      double Knvar)
{
  double u, v, z, du, dv, dz, K[4];

  if(fabs(tau) < epsilon) /* quadratic interpolation */
    return( Y0 + tau * (slope + tau * quadratic));

  /* Get the values of K and its derivatives */
  findK(tau, L, weights, size, nGroups, groupSizes, Work, Knvar, K);

  u = tau * sqrt(K[2]);
  v = (tau > 0) ? sqrt(2 * (tau * K[1] - K[0])) :
                 -sqrt(2 * (tau * K[1] - K[0]));
  z = v + log(u / v) / v;

  /* Compute the derivative of z to check for monotonicity */
  du = (2 * K[2] + tau * K[3])/ 2.0 / sqrt(K[2]);
  dv = tau * K[2] / v;
  dz = dv * (2 - 1 / v / v) + (du / u - dv * z) / v;
  if(dz < 0){
    /* printf("WARNING: approximation decreasing at tau = %f (z'(tau) = %f)\n",
       tau, dz); */
    /* printf("Consider increasing the convolution factor.\n"); */
    sprintf(message,
	    "approximation decreasing at tau = %f (z'(tau) = %f). Consider increasing the convolution factor.", tau, dz);
    PROBLEM
      message
      WARNING(NULL_ENTRY);
  }

  return z;
}

/***********************************************************************
  SaddleD
***********************************************************************/
static double SaddleD(double tau, double L[], double weights[], S_LONG size,
		      S_LONG nGroups, S_LONG groupSizes[], double Work[],
		      double Knvar)
{
  double K[4];

  /* Get the values of K and its derivatives */
  findK(tau, L, weights, size, nGroups, groupSizes, Work, Knvar, K);

  /* Return density for the sample sum */
  return exp(K[0] - tau * K[1])/sqrt(2 * PI * K[2]);
}

/***********************************************************************
  findK

  Evaluate the values of K and its first three derivatives at tau.
  Return in the vector K.
***********************************************************************/
static void findK(double tau, double L[], double weights[], S_LONG size,
		  S_LONG nGroups, S_LONG groupSizes[], double Work[], double Knvar,
		  double *K)
{
  S_LONG i, j, N, sampleSize;
  double w, sum0, sum1, sum2, sum3, temp1, temp2, mean;
  double Kn, Kn1, Kn2;

  K[0] = K[1] = K[2] = K[3] = 0.0;

  /* cgf for a single, normally distributed observation, to be added
     into the observed cgf.
  */
  Kn = .5 * (Knvar) * tau * tau;
  Kn1 = Knvar * tau;
  Kn2 = Knvar;

  /* For grouped data, the values are just the sum of the values for each
     group */

  for(j = 0; j < nGroups; j++){

    sum0 = sum1 = sum2 = sum3 = 0.0;
    N = groupSizes[j];
    /* For stratified data, sample sizes are the group sizes.  */
    sampleSize = (nGroups > 1) ? N : size;

    if(weights){
      for(i=0; i<N; i++){
	w = weights[i] * exp(tau * L[i]);
	Work[i] = w;
	sum0 += w;
	sum1 += w * L[i];
      }
      K[0] += sampleSize * log(sum0 / sumWeights[j]);
      weights += N; /* set up for next group */
    }
    else {
      for(i=0; i<N; i++){
	w = exp(tau * L[i]);
	Work[i] = w;
	sum0 += w;
	sum1 += w * L[i];
      }
      K[0] += sampleSize * log(sum0 / N);
    }
    mean = sum1 / sum0;
    K[1] += sampleSize * mean;
    for(i=0; i<N; i++){
      temp1 = L[i] - mean;
      temp2 = Work[i] * temp1 * temp1;
      sum2 += temp2;
      sum3 += temp2 * temp1;
    }
    K[2] += sampleSize * sum2 / sum0;
    K[3] += sampleSize * sum3 / sum0;

    L += N; /* set up for next group */
  }
  K[0] += Kn;
  K[1] += Kn1;
  K[2] += Kn2;

  return;
}


/***********************************************************************
  findK0

  Find the values of K and its first three derivatives at tau=0.  Return in
  vector K0.

  In fact, we only need K''(0), and K'''(0), so only K0[2] and K0[3] are
  modified.
***********************************************************************/
static void findK0(double L[], double weights[], S_LONG size,
			 S_LONG nGroups, S_LONG groupSizes[], double Work[],
			 double *K0)
{
  S_LONG i, j, N, sampleSize;
  double w, sum0, sum1, sum2, sum3, temp, temp2, mean;

  K0[2] = K0[3] = 0.0;

  /* For grouped data, the values are just the sum of the values for each
     group. */
  for(j=0; j<nGroups; j++){

    sum0 = sum1 = sum2 = sum3 = 0.0;
    N = groupSizes[j];
    /* For stratified data, sample sizes are the group sizes.  */
    sampleSize = (nGroups > 1) ? N : size;

    if(weights){
      for(i=0; i<N; i++){
	w = weights[i];
	Work[i] = w;
	sum0 += w;
	sum1 += w * L[i];
      }
      weights += N; /* set up for next group */
    }
    else {
      sum0 = N;
      for(i=0; i<N; i++){
	Work[i] = 1;
	sum1 += L[i];
      }
    }
    mean = sum1 / sum0;
    for(i = 0; i < N; i++){
      temp = (L[i] - mean);
      temp2 = Work[i] * temp * temp;
      sum2 += temp2;
      sum3 += temp2 * temp;
    }
    K0[2] += sampleSize * sum2 / sum0;
    K0[3] += sampleSize * sum3 / sum0;

    L += N; /* set up for next group */
  }

  return;
}

/***********************************************************************
  S_saddlepointP
***********************************************************************/
void S_saddlepointP(double Tau[], S_LONG *Ntau, double L[], S_LONG *Nobs,
		    S_LONG *psize, double weights[], S_LONG *Nweights,
		    S_LONG *pnGroups, S_LONG groupSizes[],
		    S_LONG *calculateMean,
		    double *pconv_factor, double alpha[])
{
  S_LONG i, M, N, size, nGroups;
  double *Work, Knvar, conv_factor;

  N           = *Nobs;
  M           = *Ntau;
  size        = *psize;
  nGroups     = *pnGroups;
  conv_factor = *pconv_factor;
  Work        = Salloc(N, double);
  if(*Nweights == 0) weights = NULL;

  initializeWeights(N, weights, nGroups, groupSizes);
  if(*calculateMean) rescaleL(N, L, nGroups, groupSizes);
  initializeCenter(N, L, weights, size, nGroups, groupSizes,
		   conv_factor, Work, &Knvar);

  for (i = 0; i < M; i++){
    alpha[i] = S_pnorm(SaddleP(Tau[i], L, weights, size,
			       nGroups, groupSizes, Work, Knvar));
  }
  return;
}


/***********************************************************************
  S_saddlepointD
***********************************************************************/
void S_saddlepointD(double Tau[], S_LONG *Ntau, double L[], S_LONG *Nobs,
		    S_LONG *psize, double weights[], S_LONG *Nweights,
		    S_LONG *pnGroups, S_LONG groupSizes[],
		    S_LONG *calculateMean,
		    double *pconv_factor, double alpha[])
{
  S_LONG i, N, M, size, nGroups;
  double *Work, K0[4], Knvar, conv_factor;

  N           = *Nobs;
  M           = *Ntau;
  size        = *psize;
  nGroups     = *pnGroups;
  conv_factor = *pconv_factor;
  Work        = Salloc(N, double);
  if(*Nweights == 0) weights = NULL;

  initializeWeights(N, weights, nGroups, groupSizes);
  if(*calculateMean) rescaleL(N, L, nGroups, groupSizes);
  findK0(L, weights, size, nGroups, groupSizes, Work, K0);
  Knvar = convolutionVariance(K0, conv_factor, size, N);

  for (i = 0; i < M; i++)
    alpha[i] = SaddleD(Tau[i], L, weights, size, nGroups, groupSizes, Work,
		       Knvar);
  if(*calculateMean){
    for (i = 0; i < M; i++)
      alpha[i] *= size;
    /* Use same factor in both group and single-sample case, because
       of how the data were rescaled in the group case. */
  }
  return;
}


/***********************************************************************
  S_saddlepointPSolve
***********************************************************************/
void S_saddlepointPSolve(double alpha[], S_LONG *Ntau, double L[], S_LONG *Nobs,
			 S_LONG *psize, double weights[], S_LONG *Nweights,
			 S_LONG *pnGroups, S_LONG groupSizes[],
			 S_LONG *calculateMean,
			 double *pconv_factor, double initial[],
			 double *ptol, double *pytol, S_LONG *MaxEval,
			 double Tau[])
{
  S_LONG N, M, size, nGroups;
  S_LONG i, iter, imax;

  double tol, ytol;
  double *Work, Knvar, conv_factor;

  N       = *Nobs;
  M       = *Ntau;
  size    = *psize;
  nGroups = *pnGroups;
  conv_factor = *pconv_factor;
  tol     = *ptol;
  ytol    = *pytol;
  Work    = Salloc(N, double);/* global variable for weights computations*/
  if(*Nweights == 0) weights = NULL;

  imax  = 0;

  initializeWeights(N, weights, nGroups, groupSizes);
  if(*calculateMean) rescaleL(N, L, nGroups, groupSizes);
  initializeCenter(N, L, weights, size, nGroups, groupSizes, conv_factor,
		   Work, &Knvar);
  for (i = 0; i < M; i++) {
    if(alpha[i] <= 0.0 || alpha[i] >= 1.0) {
      na_set3(Tau + i, S_MODE_DOUBLE, Is_NaN);
      continue;
    }
    iter = *MaxEval;
    Tau[i] = uniroot_saddle(initial[i],
			    S_qnorm( alpha[i] ),
			    N, L, weights, size, nGroups, groupSizes,
			    tol, ytol, &iter, Work, Knvar);
    if(iter > imax) imax = iter;
  }
  *MaxEval = imax;
}


/***********************************************************************
  initializeWeights

  Compute global variable sumWeights: sum of weights for each group
***********************************************************************/
static void initializeWeights(S_LONG N, double weights[], S_LONG nGroups,
			      S_LONG groupSizes[])
{
  S_LONG i, gi, j;

  if(weights){/* calculate the sum of the weights for each group */
    sumWeights = Salloc( N, double);  /* global variable set here */
    gi = 0;
    for(i=0; i < nGroups; i++){
      for(j=0; j < groupSizes[i]; j++){
	sumWeights[i] += weights[gi+j];
      }
      gi += groupSizes[i];
    }
  }
  return;
}

/***********************************************************************
  rescaleL

  Rescale L, for computations in sample mean case with multiple groups.
***********************************************************************/
static void rescaleL(S_LONG N, double L[],
		     S_LONG nGroups, S_LONG groupSizes[])
{
  S_LONG i, gi, j;
  double factor;

  if(nGroups == 1) return;  /* rescaling would do nothing */
  gi = 0;
  for(i=0; i < nGroups; i++){
    factor = N / (double) groupSizes[i];
    for(j=0; j < groupSizes[i]; j++){
      L[gi + j] *= factor;
    }
    gi += groupSizes[i];
  }
  return;
}


/***********************************************************************
  initializeCenter

  Define quantities used for quadratic interpolation for tau near 0:
     epsilon, Y0, slope, quadratic
***********************************************************************/
static void initializeCenter(S_LONG N, double L[], double weights[], S_LONG size,
			     S_LONG nGroups, S_LONG groupSizes[],
			     double conv_factor, double Work[], double *Knvar)
{
  double K0[4];
  double sigma, temp, y1, y2;

  /* interpolation based on values of K''(0) and K'''(0). */
  findK0(L, weights, size, nGroups, groupSizes, Work, K0);
  sigma = sqrt(K0[2]);  /* sigma of a sample sum */
  /* var of a single observation */
  *Knvar = convolutionVariance(K0, conv_factor, size, N);
  Y0 = K0[3] / (6 * K0[2] * sigma);

  /* SaddleP interpolates if fabs(tau) <= epsilon; prevent that now */
  epsilon = 0.0;
  temp = 0.003 / sigma;
  y1 = SaddleP(-temp, L, weights, size, nGroups, groupSizes, Work, *Knvar) /
    2.0;
  y2 = SaddleP( temp, L, weights, size, nGroups, groupSizes, Work, *Knvar) /
    2.0;
  /* divide by 2 in previous to simplify calculations below */
  /* Define domain, slope, and quadratic term for interpolation */
  epsilon = temp;
  slope = (y2 - y1) / epsilon;
  quadratic = (y1 + y2 - Y0) / (epsilon * epsilon);
  return;
}

/***********************************************************************

  convolutionVariance

  compute the variance to be used for convoution = weighted variance of
  of original distribution (K''(0)/size) times convultion factor input.

***********************************************************************/
static double convolutionVariance(double *K0, double conv_factor,
				  S_LONG size, S_LONG N)
{
  return K0[2] * conv_factor / size;
}

/***********************************************************************
  uniroot_saddle

  Solve SaddleP(tau) = Zalpha.  a is an initial guess
***********************************************************************/
static double uniroot_saddle(double a, double Zalpha,
			     S_LONG n, double L[], double weights[], S_LONG size,
			     S_LONG nGroups, S_LONG groupSizes[],
			     double tol, double ytol, S_LONG *iter,
			     double Work[], double Knvar)
{
  S_LONG   MaxIter;
  double F0, Fa, b, Fb, d, Fd, S0, x, eps;

  F0 = Y0 - Zalpha;
  if(fabs(F0) < ytol) return(0.0);
  S0 = F0 / fabs(F0);  /* This is positive when answer should be negative */

  if(fabs(a) < epsilon)
    /* If a near 0, change it to one a point used for quadratic interp */
    a = (S0 > 0) ? -epsilon : epsilon;
  else if(S0 * a > 0)
    /* If a is on the wrong side of 0, reverse it. */
    a = -a;

  MaxIter = *iter;
  Fa = SaddleP(a, L, weights, size, nGroups, groupSizes, Work, Knvar) - Zalpha;
  *iter = 1;
  if(fabs(Fa) < ytol) return(a);

  if(S0 * Fa < 0) {
    /* If 0 and a bracket the root, use bracketed secant. */
    x = bracket_saddle(a, 0.0, Fa, F0, Zalpha, n, L, weights, size,
		       nGroups, groupSizes, tol, ytol, iter, Work, Knvar);
    return x;
  }

  /* At this point a & 0 don't bracket, and a is closer. */
  b  = a;
  Fb = Fa;
  a  = 0.0;
  Fa = F0;

  tol = DOUBLE_EPS + tol;
  eps = tol * (1. + fabs(b)) / 2.;

  while(*iter <= MaxIter) {
    /* Use secant method to get new estimates;
       if estimates bracket the root then use bracketed secant. */
    d  = b - a;
    Fd = Fb - Fa;
    if(Fd / d <= 0){
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
      return b;
    }
    x  = b - 1.01 * Fb * d / Fd;
    a  = b;
    Fa = Fb;
    b  = x;
    Fb = SaddleP(b, L, weights, size, nGroups, groupSizes, Work, Knvar) -
      Zalpha;
    *iter += 1;
    if(fabs(Fb) < ytol) return(b);
    if(Fb * Fa < 0){ /* bracketed */
      x = bracket_saddle(a, b, Fa, Fb, Zalpha, n, L, weights, size,
			 nGroups, groupSizes, tol, ytol, iter, Work, Knvar);
      return x;
    }
  }

  if(fabs(d) <= eps || fabs(Fb) <= ytol) {
    return b;
  }

  /* printf("WARNING: iteration limit reached in uniroot_saddle\n"); */
  PROBLEM
	"iteration limit reached in uniroot_saddle"
	WARNING(NULL_ENTRY);
  return b;
}

/***********************************************************************
  bracket_saddle

  Solve SaddleP(tau) = Zalpha.
  a and b are endpoints of interval
  Fa = SaddleP(a) - Zalpha, Fb = SaddleP(b) - Zalpha
  Fa and Fb must have opposite signs

  This is a bracketed secant method, modified to move the next point
  slightly toward the midpoint.
***********************************************************************/
static double bracket_saddle(double a, double b, double Fa, double Fb,
			     double Zalpha,
			     S_LONG n, double L[], double weights[], S_LONG size,
			     S_LONG nGroups, S_LONG groupSizes[],
			     double tol, double ytol, S_LONG *iter, double Work[],
			     double Knvar)
{
  double d, c, Fc, lambda;

  if(Fa * (Fb / fabs(Fb)) > 0) {
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
	"internal error, roots don't bracket answer in bisect"
	WARNING(NULL_ENTRY);
  }
  tol = tol / 2.;

  while (1) {
    if(fabs(Fa) < fabs(Fb)) {
      c  =  a;
      a  =  b;
      b  =  c;
      c  = Fa;
      Fa = Fb;
      Fb = c;
    }

    if(fabs(Fb) < ytol) {
      return b;
    }

    d = b - a;

    if(fabs(d) < tol * (1. + fabs(b))) {
      return b;
    }

    lambda = Fb / (Fb - Fa);
    lambda = .5 + .99 * (lambda - .5);
    c = lambda * a + (1.0 - lambda) * b;

    Fc = SaddleP(c, L, weights, size, nGroups, groupSizes, Work, Knvar) -
      Zalpha;
    ++(*iter);

    if(fabs(Fc) < ytol) {
      return c;
    }

    if(Fa * (Fc / fabs(Fc)) <  0.) {
      b  = c;
      Fb = Fc;
    }
    else {
      a  = c;
      Fa = Fc;
    }
    if((Fb - Fa) / (b - a) < 0){
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
