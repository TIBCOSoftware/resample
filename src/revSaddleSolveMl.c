#include <S.h>
#include <machine.h>
#define NR_END 1
#define FREE_ARG void*

static S_LONG debug;
static char message[200];

/*
   Functions and subroutines in this file,
   in the order they appear, and what they call


   vector                     Memory allocate vector

   free_vector                Free memory allocated vector

   S_revSaddlepointPSolve_ml  Solve the maximum likelihood tilting premeter
                              {tau}, from the reverse saddlepoint cdf estimate
			      probability {alpha}, called from S.
     find_roots

   find_roots                 Solve {tau} by iteration with Newton method.
     upDate
     findJac

   findJac                    Estimate Jacobi matrix by difference method.
     upDate

   upDate                     Update values of zero functions, by new value
                              of input X = {tau, tau2}.
     vector
     free_vector

   ---------------------------------------------------------------------

   static variables and where they are set:

   `Nobs'     : the size of input sample `L', i.e. the `length' of `L'.
   `Nweights' : the `length' of input `Weights' associated with `L',
                either the same as `Nobs' if `Weights' supplied,
                or zero for `NULL' `Weights'.
   `NiterMax' : the maximum number allowed for iteration.
   `Q', the supplied numerical value, treated as the mean of 2nd tilting.
   `Alpha', the supplied `revSaddlepoint properbility'.
   `xTol', the supplied tolrence for solution `{tau, tau2}', in meaning of
    absolute norm, default is 1e-6.
   `fTol', the supplied tolrence for zero functions, in meaning of absolute
    norm, default is 1e-8.
   `sumW', sum of `Weights', it will be defined as `(double)Nobs',
    if no `Weights' supplied.
   `sumL', sum of `L', or weighted sum of `L', if `Weights' supplied.
   `meanL' sample mean of `L', or weighted sample mean, if `Weights' supplied.
   `upB', `lowB' are the boundary for `tau', since in ML tilting case,
    it is required that `(1/min(L-meanL)) < tau < (1/max(L-meanL))', and
    upbound = 1/max(L-meanL), lowbound = 1/min(L-meanL).

   ---------------------------------------------------------------------

   The Newton method in this file is based on following relationship :

   W = tiltWeights(Tau, L, tilt = "ml", weights = Weights)
   Q = tiltMean(Tau2, L, weights = W)
   Alpha = 1 - saddlepointP(t, L, size = Nobs, weights = W)

   Assign:    X = {Tau, Tau2}
              Fvec(X) = { (Q - S_tiltMean_exp(Tau2, L, W(Tau))),
                          (S_saddlepointP(Tau2, L, W(Tau)) - (1 - Alpha)) }

   Target:    Find X to make Fvec = 0.

   Notice:    Inorder to assure positive W, the initial of Tau is required
              to be within a region specified by input L and Weights.
	      Otherwise, as if the given Alpha out of region of (0,1), a
	      value of {NA, NA} for {Tau, Tau2} will be returned.

*/

static S_LONG Nobs, Nweights, NiterMax;
static double Q, *L, *Weights, upB, lowB;
static double sumW, sumL, meanL, xTol, fTol;

double *vector(S_LONG nl, S_LONG nh)
{
  double *v;

  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if(!v){
    sprintf(message, "\n allocation failure in vector()\n");
    PROBLEM message WARNING(NULL_ENTRY);
  }
  return v-nl+NR_END;
}

void free_vector(double *v, S_LONG nl, S_LONG nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}


void S_revSaddlepointPSolve_ml(S_LONG *pNobs, S_LONG *pNtau, S_LONG *pNweights,
			       double *pL, double *pWeights,
			       double *Alpha, double *pQ,
			       double *Tau0, double *Tau20, double *pXtol,
			       double *pFtol, S_LONG *pMaxIter, S_LONG *pDebug,
			       S_LONG *pAccIter,
			       double *Tau, double *Tau2)
{
  /*
     Define following auxiliary variables:
     `Ntau' : the size of `Tau', it was copied from the length of `Alpha'
     `i' : index for `for' loop;
     `Niter' : actual iteration number for each `Alpha';
     `imax' : maximum iteration number for any of `Alpha'.
  */

  S_LONG Ntau;
  S_LONG i, j, Niter, imax;
  double alpha, Ldiff;

  void find_roots(double, double*, double*, double, double, S_LONG*);

  /* Define values for global variables: */

  Nobs = *pNobs;
  Ntau = *pNtau;
  Nweights = *pNweights;
  NiterMax = *pMaxIter;
  Q = *pQ;
  xTol = *pXtol;
  fTol = *pFtol;
  L = pL;
  Weights = pWeights;
  debug = *pDebug;
  if(debug > 0) printf("Entering S_revSaddlepointPSolve_ml\n");

  /* sumW, sumL and meanL will have to change for groups */
  sumW = 0;
  sumL = 0;
  if(Nweights == 0) {
    Weights = NULL;
    sumW = (double)Nobs;
  }
  else
    for(i=0; i < Nobs; i++)
      sumW += Weights[i];
  for(i=0; i < Nobs; i++) {
    if(Nweights == 0)
      sumL += L[i];
    else
      sumL += Weights[i] * L[i];
  }
  meanL = sumL / sumW;

  /* Find the boundary of `tau' in order to assure positive tilt weights : */
  /* This will have to change for groups:  loop over groups */
  upB = 0;
  lowB = 0;
  for(i=0; i < Nobs; i++){
    Ldiff = (L[i] - meanL);
    if(upB < Ldiff)
      upB = Ldiff;
    if(lowB > Ldiff)
      lowB = Ldiff;
  }
  upB = 1/upB;
  lowB = 1/lowB;

  /* Find solution for each alpha (revSaddlepoint probability) : */

  imax  = 0;
  alpha = 0;
  Niter = 0;

  for (j = 0; j < Ntau; j++) {
    alpha = Alpha[j];
    if(alpha <= 0.0 || alpha >= 1.0) {
      sprintf(message, "The probs = %f is out of (0, 1). \n", alpha);
      PROBLEM
	message
	WARNING(NULL_ENTRY);
      na_set3(Tau + j, S_MODE_DOUBLE, Is_NaN);
      na_set3(Tau2 + j, S_MODE_DOUBLE, Is_NaN);
    }
    else {
      Niter = NiterMax;
      find_roots(alpha, &Tau[j], &Tau2[j], Tau0[j], Tau20[j], &Niter);
      if(Niter > imax)
	imax = Niter;
    }
  }

  *pAccIter = imax;
  if(debug > 0) printf("Exiting S_revSaddlepointPSolve_ml\n");
  return;
}


void find_roots(double alpha, double *ptau, double *ptau2,
		double tau0, double tau20, S_LONG *iter)
{
  S_LONG its, MAXITS, i;
  double test, temp, lumda, delta, pTemp;

  double x[2], xold[2], fvec[2], fvecold[2], p[2], fjac[2][2];

  void upDate(double*, double, double*);   /* update `fvec' by values of `x' */
  void findJac(double*, double*, double, double fjac[2][2]);/* estimate fjac */

  MAXITS = *iter;
  x[0] = tau0;
  x[1] = tau20;

  /* test if input of `tau' is out of a required region, return NA if so: */

  if(x[0] > upB || x[0] < lowB) {
    sprintf(message, "For probs = %f, the initial: `tau' = %f  was out of (%f, %f), and caused negative tilt weight.", alpha, x[0], lowB, upB);
    PROBLEM
      message
      WARNING(NULL_ENTRY);
    na_set3(x, S_MODE_DOUBLE, Is_NaN);
    na_set3(x+1, S_MODE_DOUBLE, Is_NaN);

    *iter = 0;
    *ptau = x[0];
    *ptau2 = x[1];

    return;
  }

  upDate(x, alpha, fvec);

  /* test if the initial is a solution, in the meaning of within tolrence : */

  test = fabs(fvec[0]);
  if(fabs(fvec[1]) > test)
    test = fabs(fvec[1]);
  if (test < fTol) {
    /*
    PROBLEM
      "Initial pick is a root, no more iteration."
      WARNING(NULL_ENTRY);
    */
    *iter = 0;
    *ptau = x[0];
    *ptau2 = x[1];

    return;
  }

  /* Iteration : */

  for(its=0; its<MAXITS; its++){
    for(i=0; i<2; i++){
      xold[i] = x[i];
      p[i] = -fvec[i];
    }

    findJac(x, fvec, alpha, fjac);      /* estimate Jacobi */

    /* solve p as new shift for x, so it may change to the desired root :  */

    delta = fjac[0][0]*fjac[1][1] - fjac[0][1]*fjac[1][0];
    if(delta == 0.0)
      PROBLEM
	"The Jacobi matrix is a Singular matrix."
	RECOVER(NULL_ENTRY);
    pTemp = (fjac[0][0]*p[1] - p[0]*fjac[1][0]) / delta;
    p[0]  = (p[0]*fjac[1][1] - fjac[0][1]*p[1]) / delta;
    p[1]  = pTemp;

    /* modified to assure iteration stay with the boundary of `tau' : */

    temp = xold[0] + p[0];
    lumda = 1;

    if(temp >= upB) {
      lumda = 0.9 * (upB - xold[0]) / p[0];
      p[0] *= lumda;
      p[1] *= lumda;
      x[0] = xold[0] + p[0];
      x[1] = xold[1] + p[1];
    }
    else if(temp <= lowB) {
      lumda = 0.9 * (lowB - xold[0]) / p[0];
      p[0] *= lumda;
      p[1] *= lumda;
      x[0] = xold[0] + p[0];
      x[1] = xold[1] + p[1];
    }
    else {
      x[0] = xold[0] + p[0];
      x[1] = xold[1] + p[1];
    }

    /* test if the X converges : */

    test = fabs(p[0]);
    if(fabs(p[1]) > test)
      test = fabs(p[1]);
    if(test < xTol) {
      *iter = its;
      /*
      PROBLEM
	"X converges"
	WARNING(NULL_ENTRY);
      */
      *ptau = x[0];
      *ptau2 = x[1];

      return;
    }

    for(i = 0; i < 2; i++)
      fvecold[i] = fvec[i];

    upDate(x, alpha, fvec);

    /* test if the new x is a solution, in the meaning of within tolrence : */

    test = fabs(fvec[0]);
    if(fabs(fvec[1]) > test)
      test = fabs(fvec[1]);
    if(test < fTol) {
      *iter = its;
      *ptau = x[0];
      *ptau2 = x[1];

      return;
    }

    /* test if the Fvec converges : */

    test = fabs(fvec[0] - fvecold[0]);
    if(fabs(fvec[1] - fvecold[1]) > test)
      test = fabs(fvec[1] - fvecold[1]);
    if(test < fTol) {
      *iter = its;
      /*
      PROBLEM
	"Function converges"
	WARNING(NULL_ENTRY);
      */
      *ptau = x[0];
      *ptau2 = x[1];

      return;
    }
  }

  *iter = its;

  PROBLEM
    "Maximum iteration is exceeded"
    WARNING(NULL_ENTRY);

  *ptau = x[0];
  *ptau2 = x[1];

  return;
}


void findJac(double x[2], double fvec[2], double alpha, double fjac[2][2])
{
  /* forward-difference approximation to Jacobian */

  /* Parameter epsilon */
  static double epsilon = 1E-6;

  S_LONG j, k;
  double h, temp, fvecTemp[2];

  void upDate(double*, double, double*);

  for(j = 0; j < 2; j++) {
    temp = x[j];
    h = epsilon * fabs(temp);
    if(h == 0.0)
      h = epsilon;
    x[j] = temp + h;

    upDate(x, alpha, fvecTemp);

    x[j] = temp;
    for(k=0; k < 2; k++)
      fjac[k][j] = (fvecTemp[k]/h) - (fvec[k]/h);
  }

  return;
}


void upDate(double x[2], double alpha, double fvec[2])
{
  S_LONG i;
  double tau, tau2, *W;

  S_LONG Nt, pnG;
  double Tau2Q, pConv, saddleP;
  S_LONG zero = 0;

  void S_tiltMean_exp(S_LONG *Nobs, S_LONG *Ntau, S_LONG *Nw,
		      double *L, double *weights,
		      S_LONG *pnGroups, S_LONG groupSizes[],
		      S_LONG *pDebug,
		      double *Tau);
  void S_saddlepointP(double Tau[], S_LONG *Ntau, double L[], S_LONG *Nobs,
		      S_LONG *psize, double weights[], S_LONG *Nweights,
		      S_LONG *pnGroups, S_LONG groupSizes[],
		      S_LONG *calculateMean,
		      double *pconv_factor, double alpha[]);

  W = vector(0, Nobs-1);

  /* update W = ml tilting weights: */

  tau = x[0];

  /* The following must change for groups */
  for(i=0; i<Nobs; i++){

    W[i] = 1 / (1 - tau * (L[i] - meanL));

    if(W[i] < 0) {
#ifdef WIN64 //S7_Change_LONG
      sprintf(message, "The %lld th tilting weight is negative, need to check. The guessed \"tau\" is %f. ", i, tau);
#else
      sprintf(message, "The %ld th tilting weight is negative, need to check. The guessed \"tau\" is %f. ", i, tau);
#endif //WIN64
      PROBLEM
	message
	RECOVER(NULL_ENTRY);
    }

    if(Nweights)
      W[i] *= Weights[i];
  }

  /* update the tilting back mean: */

  tau2 = x[1];
  Nt = (S_LONG)1;
  pnG = (S_LONG)1;

  Tau2Q = tau2;

  /* Following will have to change for groups */
  /* adopt `S_tiltMean_exp', match case of only one group. */
  S_tiltMean_exp(&Nobs, &Nt, &Nweights, L, W, &pnG,
		 &Nobs, &debug, &Tau2Q);

  /* Note: the `&Tau2Q' input tau2, and return tiltMean. */

  fvec[0] = Tau2Q - Q;

  /* update the tilting back saddle: */

  saddleP = 0;
  pConv = 0;

  /* Following will have to change for groups */
  /* adopt `S_saddlepointP', match case of only one group. */
  S_saddlepointP(&tau2, &Nt, L, &Nobs, &Nobs, W, &Nweights,
		 &pnG, &Nobs, &pnG, &pConv, &saddleP);

  fvec[1] = saddleP - (1 - alpha);
  free_vector(W, 0, Nobs-1);

  return;
}
