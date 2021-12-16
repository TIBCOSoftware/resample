#include <S.h>
#include <machine.h>
static S_LONG debug;  /* global within this file*/
#include "tilt.h"

/* Functions and subroutines in this file:
   (entry points are marked with *)

   findDomain         Find range of L, and domain for ML tau
   makeSubsets        Two subsets of 1:B, with theta <= observed, or >=

   tiltBoot_exp       given tau, find prob
   tiltBoot_ml        given tau, find prob

   my_qnorm           ifelse(x < .5, qnorm(x), linear approx)

 * S_tiltBoot_exp     given taus, find probs; calls tiltBoot_exp
 * S_tiltBoot_ml      given taus, find probs; calls tiltBoot_ml
 * S_tiltBootSolve    given probs, find taus

   solve_tiltBoot_exp
   solve_tiltBoot_ml
   bisect_tiltBoot_exp
   bisect_tiltBoot_ml

   --------------------------------------------------
   Functions in tilt.h -- must be included

   vectorMin
   vectorMax
   weightedMean       weighted mean of a vector
   sumWeightedMeans   sum of weighted group means
   checkLambda
   mlg_setup          compute Ltau, find lambdas


   --------------------------------------------------
   Functions in tiltMean.c -- must be linked:
   S_tiltMean_exp       tilted mean

   --------------------------------------------------
   Functions in basic.c -- must be linked:
   S_bootstrapSums      bootstrap sums:     answer[i,j] = sum(x[indices[,i],j])
   S_bootstrapMeans     bootstrap means:    answer[i,j] = mean(x[indices[,i],j])
   S_bootstrapProducts  bootstrap products: answer[i,j] = prod(x[indices[,i],j])
   S_bootstrapProducts_log  bootstrap products, but use logs


   --------------------------------------------------
   Calling trees (not including calls to various utility functions)
   S_tiltBoot_*
     tiltBoot_*

   S_tiltBootSolve
     solve_tiltBoot_*
       tiltBoot_*
       bisect_tiltBoot_*
         tiltBoot_*

   --------------------------------------------------

   Variables used throughout this file:

   -- Dimension variables --
   N = *Nobs = sample size
   B = *Nboot = number of bootstrap replications
   P = *Nstat = number of statistics
   K = *Ntau = *Nprob = number of probabilities or values of tau
   B2 = number of bootstrap replications in a subset, determined by either
        (replicates >= observed) or (replicates <= observed); roughly B/2
   nGroups = *pnGroups = number of groups
   groupSizes[g] = number of observations in group g

   -- Input vectors and matrices --
   L[N, P]          (in S_tiltBootSolve)
   L[N]             (elsewhere)
   Replicates[B, P] (in S_tiltBootSolve)
   Replicates[B]    (elsewhere)
   Observed[P]
   expTau[P, K]
   mlTau[P, K]
   Indices[N, B]    (in the range 1..N -- S convention, not C convention)
   weights[B]       (weights on bootstrap samples, not observations)

   -- Index variables --
   i = 1 to N       (0 to (N-1))
   j = 1 to P       (0 to (P-1))
   k = 1 to K       (0 to (K-1))
   b = 1 to B or B2 (0 to (B-1))       (b is an interval endpoint in solve*)
   g = 1 to nGroups (0 to (nGroups-1))
   i = 1 to groupSizes[g]

   -- Other quantities --
   tau = scalar, tilting parameter
   Tau[] = vector
   xtol = on the scale of tau or lambda
   ytol = on the scale of tilted mean, sum of probs, or bootstrap prob
*/




/* Global variables (also `debug', above) */
static char message[200];
static S_LONG nGroups;
static S_LONG *groupSizes;   /* [nGroups] */
static double *Lambda;     /* [nGroups], normalizing constants for ML tilting */
static double *V0;         /* [N],       V0[gi] = 1/groupSizes[g]; pass to mlg_setup */


/* Prototypes for functions external to this file. */
void S_bootstrapSums(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		     double *x, S_LONG *indices,
		     double *answer);
void S_tiltMean_exp(S_LONG *Nobs, S_LONG *Ntau, S_LONG *Nw,
		    double *L, double *weights,
		    S_LONG *pnGroups, S_LONG groupSizes[],
		    S_LONG *pDebug,
		    double *Tau);


/* Prototypes for functions in this file. */
static void findDomain(S_LONG, double *, double *, double *);
static void makeSubsets(S_LONG, double *, double,	S_LONG *, S_LONG *, S_LONG *);
static double tiltBoot_exp(double tau,
			   S_LONG N, double *L,
			   S_LONG B, double *Sj, double *weights,
			   S_LONG B2, S_LONG *subset);
static double tiltBoot_ml(double tau,
			  S_LONG N, double *L, double *Ltau,
			  S_LONG B, S_LONG *Indices, double *bweights,
			  S_LONG B2, S_LONG *subset,
			  double *W,
			  double ytol, S_LONG maxiter);
static double my_qnorm(double);
static double solve_tiltBoot_exp(double a, double F0,
				 double v,
				 double xtol, double ytol, S_LONG *iter,
				 S_LONG sgn,
				 S_LONG N, double *L,
				 S_LONG B, double *Sj, double *weights,
				 S_LONG B2, S_LONG *subset);
static double solve_tiltBoot_ml(double a, double F0,
				double v,
				double xtol, double ytol, S_LONG *iter,
				double tauneg, double taupos,
				S_LONG sgn,
				S_LONG N, double *L, double *Ltau,
				S_LONG B, S_LONG *Indices, double *weights,
				S_LONG B2, S_LONG *subset,
				double *W);
static double bisect_tiltBoot_exp(double a, double b, double Fa, double Fb,
				  double v,
				  double xtol, double ytol, S_LONG *iter,
				  S_LONG N, double *L,
				  S_LONG B, double *Sj, double *weights,
				  S_LONG B2, S_LONG *subset);
static double bisect_tiltBoot_ml(double a, double b, double Fa, double Fb,
				 double v,
				 double xtol, double ytol, S_LONG maxiter,
				 S_LONG *iter,
				 S_LONG N, double *L, double *Ltau,
				 S_LONG B, S_LONG *Indices, double *weights,
				 S_LONG B2, S_LONG *subset,
				 double *W);


/*--------------------------------------------------*/



static void findDomain(S_LONG N, double *L,
		       double *range, double *domain)
{
  /* Find range of L and domain = (1/min(L), 1/max(L)). */

  /* range and domain may be NULL; if so don't assign them. */
  double Lmin, Lmax;

  if(debug > 0) printf("Entering findDomain{\n");
  Lmin = vectorMin(N, L);
  Lmax = vectorMax(N, L);
  if(range){
    range[0] = Lmin;
    range[1] = Lmax;
  }
  if(domain){
    domain[0] = (Lmin < 0.0) ? 1.0/Lmin : -HUGE_VAL;
    domain[1] = (Lmax > 0.0) ? 1.0/Lmax : HUGE_VAL; /* Infinity */
  }
  if(debug > 0) printf("}Exiting findDomain\n");
  /* Use this code if need to make domain finite (shouldn't need to).
  domain[0] = (Lmin < 0.0) ? 1.0/Lmin : -log(HUGE)/2.0;
  domain[1] = (Lmax > 0.0) ? 1.0/Lmax : log(HUGE)/2.0;
  */  /* HUGE=3.4e+38 */
}



static void makeSubsets(S_LONG B, double *Replicates, double observed,
			S_LONG *subset, S_LONG *B2_lo, S_LONG *B2_hi)
{
  /* Determine two subsets of (1:B), one for which
     (Replicates <= Observed) and one for which
     (Replicates >= Observed).
     subset is of length B.
     Put indices for the first subset in the first B2_lo positions of
     subset, and indices for the second subset in the last B2_hi positions,
     so tie cases are in the middle.
  */
  S_LONG lo, hi, eq, b;
  double difference;

  if(debug > 0) printf("Entering makeSubsets{\n");
  lo = 0;
  hi = B;
  eq = 0;
  for(b = 0; b < B; b++) {
    difference = Replicates[b] - observed;
    if(difference < 0.0)
      subset[lo++] = b;
    else if(difference > 0.0)
      subset[--hi] = b;
    else
      ++eq;  /* Place these later. */
  }

  *B2_lo = lo + eq;
  *B2_hi = B - lo;

  /* If there were equality cases, handle them now. */
  if(eq)
    for(b = 0; b < B; b++)
      if(Replicates[b] == observed)
	subset[lo++] = b;
  /* printf("B2_lo = %ld, B2_hi = %ld, indices are:\n", *B2_lo, *B2_hi); */
  /* for(b=0; b<B; b++) printf("%ld ", subset[b]); printf("\n"); */
  if(debug > 0) printf("}Exiting makeSubsets\n");
}




static double tiltBoot_exp(double tau,
			   S_LONG N, double *L,
			   S_LONG B, double *Sj, double *weights,
			   S_LONG B2, S_LONG *subset)
{
  /* tau     [1]  = tilting parameter (scalar)
     L       [N]  = values to be tilted
     Sj      [B]  = bootstrap sums of L
     weights [B]  = weights - don't normalize, mean about 1
     subset  [B2] = indices, for subscripting Sj & weights
         the subset is for replicates either >= or <= observed
     Some global variables are also used.

     Return bootstrap tilting probability,
         (1/B) sum_{subset of b's} weights_b prod_{i=1}^N w^*_{i,b}
     where w_i = c exp(tau * L_i), and c normalizes to mean 1.

     Note that
         prod w^*_{i,b} = exp(N log(c) + sum(L^*_{i,b}));
     the latter sum is in Sj.  Hence in this function we don't need
     the indices themselves.
  */

  S_LONG g, gi, i, b;
  double sum, tmax, normalizingTerm;

  if(debug>1) printf("Entering tiltBoot_exp{\n");
#ifdef WIN64 //S7_Change_LONG
  if(debug>3) printf("  tau=%g, N=%lld, B=%lld, B2=%lld\n", tau, N, B, B2);
#else
  if(debug>3) printf("  tau=%g, N=%ld, B=%ld, B2=%ld\n", tau, N, B, B2);
#endif //WIN64
  if(debug>4) printVectorDouble("  L", N, L);
  if(debug>4) {
    if(weights) printVectorDouble("  weights", N, weights);
    else printf("  No weights\n");
  }
  if(debug>7) printVectorDouble("  Sj", B, Sj);
  if(debug>7) printVectorLong("  subset", B2, subset);
      
  if(tau == 0){
    if(!weights)
      return((double)B2/(double)B);
    sum = 0.0;
    for(b = 0; b < B2; b++) sum += weights[subset[b]];
    return(sum/(double)B);
  }

  /* Use a numerically stable method, in case tau is large.  Start by
     finding max(tau*L) (over all groups)
  */
  if(tau < 0)
    tmax = tau * vectorMin(N, L);
  else
    tmax = tau * vectorMax(N, L);


  if(nGroups == 1){
    /* sum = \sum_{i=1}^N \exp(tau * L_i - tmax) */
    sum = 0.0;
    for(i = 0; i < N; i++)
      sum += exp(tau * L[i] - tmax);
    /* Let c = N/(sum * exp(tmax)), c makes weights have mean 1, we want
       prod (c * exp(L * tau))
         = c^N exp(tau * sum(L))
         = exp(tau * sum(L) + normalizingTerm)
       normalizingTerm = log(c^N) = N log(c)
    */
    normalizingTerm = (double)N * (log((double)N / sum) - tmax);
  }
  else {
    /* Let 
         sum_g = sum_{group g} exp(tau * L_i - tmax)
       then the constant that makes group g have mean 1 is
         N_g / (sum_g * exp(tmax))
       and we need to compute
         prod_g [N_g / (sum_g * exp(tmax))]^{N_g} exp(tau * sum(L))
	 = prod_g (N_g/sum_g)^{N_g}   exp(-N tmax) exp(tau * sum(L))
       normalizingTerm = sum_g N_g log(N_g/sum_g) - N tmax
    */
    normalizingTerm = 0.0;
    gi = 0;
    for(g = 0; g < nGroups; g++){
      sum = 0.0;
      for(i = 0; i < groupSizes[g]; i++)
	sum += exp(tau * L[gi++] - tmax);
      normalizingTerm += groupSizes[g] * log((double) groupSizes[g] / sum);
    }
    normalizingTerm -= (double) N * tmax;
  }
  if(debug > 2) printf("  normalizingTerm=%g\n", normalizingTerm);

  sum = 0.0;
  if(weights)
    for(b = 0; b < B2; b++) sum += exp(tau * Sj[subset[b]] + normalizingTerm)
			      * weights[subset[b]];
  else
    for(b = 0; b < B2; b++) sum += exp(tau * Sj[subset[b]] + normalizingTerm);

  if(debug>1) printf("}Exiting tiltBoot_exp\n");
  return(sum/(double)B);
}



static double tiltBoot_ml(double tau,
			  S_LONG N, double *L, double *Ltau,
			  S_LONG B, S_LONG *Indices, double *bweights,
			  S_LONG B2, S_LONG *subset,
			  double *W,
			  double ytol, S_LONG maxiter)
{
  /* tau      [1]   = tilting parameter (scalar)
     L        [N]   = values to be tilted
     Ltau     [N]   = L*tau[g] (for group ML tilting only, else NULL pointer)
     Indices  [N,B] = bootstrap indices
     bweights [B]   = bootstrap weights - don't normalize, mean about 1
     subset   [B2]  = indices, for subscripting Indices & bweights
         the subset is for replicates either >= or <= observed
     W        [N]   = work space
     ytol, maxiter = govern numerical search for Lambda

     Return bootstrap tilting probability,
         (1/B) sum_{subset of b's} bweights_b prod_{i=1}^N (w^*_{i,b})
     where 
         w_i = c /(1 - tau * L_i), and c normalizes to mean 1.
     In the multiple group case
         w_i = 1 /(Lambda[g] - tau[g] * L_i)
         tau[g] = tau / (groupSizes[g]/N)
     and Lambda[g] normalizes to mean 1.

     Note:  w_i have sum=1 or mean=1 in each group, depending
     on whether numerator is 1/groupSizes[g] or 1.
  */

  S_LONG i, b, g, gi, Ng;
  S_LONG *Indices_b;
  double tmax, temp, minimumDenominator, denominator, sum, lambdag, logc;
  double *W_offset;

#ifdef WIN64 //S7_Change_LONG
  if(debug>1) printf("Entering tiltBoot_ml{, tau=%g, N=%lld, B=%lld, B2=%lld\n",
		     tau, N, B, B2);
#else
  if(debug>1) printf("Entering tiltBoot_ml{, tau=%g, N=%ld, B=%ld, B2=%ld\n",
		     tau, N, B, B2);
#endif //WIN64
  if(debug>5) {
    printVectorDouble("  L", N, L);
    if(bweights) printVectorDouble("  bweights", B, bweights);
    printVectorLong("  subset", B2, subset);
  }
  if(debug>6)
    printVectorLong("  Indices", N*B, Indices);

  if(tau == 0){
    if(!bweights)
      return((double)B2/(double)B);
    sum = 0.0;
    for(b = 0; b < B2; b++) sum += bweights[subset[b]];
    return(sum/(double)B);
  }

  /* Compute log(w_i), store in W[i]. */
  if(nGroups == 1){   /*---------------single-group case-----------------*/
    /* Use a numerically stable method, in case tau is close to one
       end of its range.  Start by finding max(tau*L). */
    if(tau < 0)
      tmax = tau * vectorMin(N, L);
    else
      tmax = tau * vectorMax(N, L);
    if(tmax >= 1.0){
      na_set3(&denominator, S_MODE_DOUBLE, Is_NaN);
      PROBLEM
	"tau out of range to keep weights positive; setting results to NA"
	WARNING(NULL_ENTRY);
      return(denominator);
    }

    minimumDenominator = 1.0 - tmax;
    sum = 0.0;
    for(i = 0; i < N; i++) {
      denominator = 1.0 - tau * L[i];
      W[i] = denominator;         /* Temporarily W contains denominators */
      sum += minimumDenominator/denominator; /* This weight over largest. */
    }
    /* normalizing constant is c = N * minimumDenominator / sum
       The true weights are c / W[i].
       Find   prod((bootstrap sample of true weights))
    */

    /* Let W[i] = log(weight normalized to mean 1),
       so the product for boot sample b is exp(sum(W^*_{i,b})) */
    temp = N * minimumDenominator / sum;
    for(i=0; i<N; i++)
      W[i] = log(temp / W[i]);
  }
  else {  /*-------------------------multiple-group case---------------*/
    /* Compute Ltau = L * tau/h_i, and lambdas */
    if(debug > 4) printVectorDouble("  V0 in tiltBoot_ml(1)", N, V0);
    mlg_setup(tau, L, V0, N, nGroups, groupSizes, 0,
	      ytol, maxiter, Lambda, Ltau);
    if(debug > 4) printVectorDouble("  V0 in tiltBoot_ml(2)", N, V0);

    /* Let W[i] = log( c / (Lambda[g] - tau[g] * L_i))
       where c normalizes to mean 1.
       Algebraically c=1 after mlg_setup, but do not count on that. */
    for(g = 0, gi = 0; g < nGroups; g++){
      Ng = groupSizes[g];
      lambdag = Lambda[g];
      sum = 0.0;
      for(i = 0; i < Ng; i++)
	sum += 1/(lambdag - Ltau[gi++]);
      gi -= Ng;
      logc = log( (double)Ng / sum);
      for(i = 0; i < Ng; i++, gi++)
	W[gi] = logc - log(lambdag - Ltau[gi]);
    }
    if(debug > 4) printVectorDouble("  V0 in tiltBoot_ml(3)", N, V0);
  }

  /* Final calculations : */

  /* Since Indices are in 1:N, not 0:N-1, use an offset */
  W_offset = W - 1;

  sum = 0.0;
  for(b = 0; b < B2; b++) {
    temp = 0.0;
    Indices_b = Indices + N * subset[b]; /* Pointer to column of Indices */
    for(i = 0; i<N; i++)
      /* temp += W[ Indices_b[i] ]; */
      temp += W_offset[ Indices_b[i] ];
    if(bweights)
      sum += exp(temp) * bweights[subset[b]];
    else
      sum += exp(temp);
  }
  if(debug>1) printf("}Exiting tiltBoot_ml\n");
  return(sum/(double)B);
}




static double my_qnorm(double x)
{
  /* If x<=.5, return qnorm(x), else use linear approx
     If x>=1, print warning.
  */
  if(debug>2) printf("Entering my_qnorm\n");
  if(x <= 0.5)
    return(S_qnorm(x));
  if(x >= 1.0) {
    /*    printf("Warning, requested qnorm(p) with p = %g > 1\n",x); */
    sprintf(message, "requested qnorm(p) with p = %g > 1", x);
    PROBLEM
      message
      WARNING(NULL_ENTRY);
  }
  return(2.506628 * (x-0.5));
}


void S_tiltBoot_exp(S_LONG *Nobs,
		    S_LONG *Nboot,
		    S_LONG *Ntau,
		    double *Observed,
		    double *Replicates,
		    double *L,
		    S_LONG *Indices,
		    S_LONG *Nweights,
		    double *weights,
		    S_LONG *pnGroups,
		    S_LONG pGroupSizes[],
		    S_LONG *pDebug,
		    double *Tau)
     /* Calculate bootstrap tilting probabilities.
	See documentation for tiltBoot_exp

	L should have mean 0 (in each group).

	Tau contains tau on input and probabilities on output.
	L is also modified (when there are multiple groups).
     */
{
  S_LONG N, B, K;
  S_LONG one = 1;
  S_LONG g, gi, i, k;
  S_LONG B2_lo, B2_hi;
  double factor, tau;
  S_LONG   *subset, *subset_hi;
  double *Sj;

  N = *Nobs;
  B = *Nboot;
  K = *Ntau;
  if(*Nweights == 0) weights = NULL;
  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltBoot_exp{\n");
#ifdef WIN64 //S7_Change_LONG
  if(debug>3) printf("  N=%lld, B=%lld, K=%lld, NWeights=%lld, nGroups=%lld, Observed=%g\n",
		     N, B, K, *Nweights, *pnGroups, *Observed);
#else
  if(debug>3) printf("  N=%ld, B=%ld, K=%ld, NWeights=%ld, nGroups=%ld, Observed=%g\n",
		     N, B, K, *Nweights, *pnGroups, *Observed);
#endif //WIN64
  if(debug>4) printVectorDouble("  Replicates", B, Replicates);
  if(debug>3) printVectorDouble("  L", N, L);
  if(debug>3) printVectorDouble("  tau", K, Tau);

  subset = Salloc( B, S_LONG);
  Sj   = Salloc( B, double);

  nGroups = *pnGroups;      /* Set global variables */
  groupSizes = pGroupSizes;

  makeSubsets(B, Replicates, *Observed, subset, &B2_lo, &B2_hi);
  subset_hi = subset + (B - B2_hi);
#ifdef WIN64 //S7_Change_LONG
  if(debug > 2) printf("  B=%lld, B2_lo=%lld, B2_hi=%lld\n", B, B2_lo, B2_hi);
#else
  if(debug > 2) printf("  B=%ld, B2_lo=%ld, B2_hi=%ld\n", B, B2_lo, B2_hi);
#endif //WIN64
  if(debug > 4) printVectorLong("  subset", B, subset);

  if(nGroups > 1){
    /* Divide L values by groupSizes/N, instead of later using tau[g] */
    for(g = 0, gi = 0; g < nGroups; g++){
      factor = (double) N / (double) groupSizes[g];
      i = groupSizes[g];
      for(i=0; i<groupSizes[g]; i++)
	L[gi++] *= factor;
    }
    if(debug>3) printVectorDouble("  L/(Ng/N)", N, L);
  }
  S_bootstrapSums(&N, &one, &N, &B, L, Indices, Sj);
  if(debug > 4) printVectorDouble("  Sj", B, Sj);

  for(k = 0; k < K; k++) {
    tau = Tau[k];
    if(tau < 0.0)
      Tau[k] =       tiltBoot_exp(tau, N, L, B, Sj, weights, B2_hi, subset_hi);
    else
      Tau[k] = 1.0 - tiltBoot_exp(tau, N, L, B, Sj, weights, B2_lo, subset);
  }
  if(debug>0) printf("}Exiting S_tiltBoot_exp\n");
  return;
}



void S_tiltBoot_ml(S_LONG *Nobs, S_LONG *Nboot, S_LONG *Ntau,
		   double *Observed, double *Replicates,
		   double *L, S_LONG *Indices,
		   S_LONG *Nweights, double *weights,
		   S_LONG *pnGroups, S_LONG pGroupSizes[],
		   double *pytol, S_LONG *pmaxiter,
		   S_LONG *pDebug,
		   double *Tau)
     /* Calculate bootstrap tilting probabilities.
	See documentation for tiltBoot_ml.

	Tau contains tau on input and probabilities on output.
     */
{
  S_LONG N = *Nobs;
  S_LONG B = *Nboot;
  S_LONG K = *Ntau;
  double ytol  = *pytol;
  S_LONG maxiter = *pmaxiter;
  S_LONG g, i, gi, k;
  S_LONG B2_lo, B2_hi;
  double tau, domain[2];
  S_LONG   *subset, *subset_hi;
  double *W, *Ltau = (double *) NULL;

  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltBoot_ml{\n");
  if(*Nweights == 0) weights = NULL;

  nGroups = *pnGroups;      /* Set global variables */
  groupSizes = pGroupSizes;

  subset = Salloc( B, S_LONG);
  W    = Salloc(N, double);
  if(nGroups > 1){
    Ltau = Salloc(N, double);
    V0   = Salloc(N, double);
    Lambda = Salloc(nGroups, double);
    /* Initialize V0 to contain observation weights, sum to 1 in each group. */
    for(g = 0, gi = 0; g < nGroups; g++)
      for(i = 0; i < groupSizes[g]; i++)
	V0[ gi++ ] = 1.0 / (double) groupSizes[g];
  }

  makeSubsets(B, Replicates, *Observed, subset, &B2_lo, &B2_hi);
  subset_hi = subset + (B - B2_hi);

  findDomain(N, L, NULL, domain);


  for(k = 0; k < K; k++) {
    tau = Tau[k];
    if(tau <= domain[0] || tau >= domain[1])
      na_set3(Tau+k, S_MODE_DOUBLE, Is_NaN);
    else if(tau < 0.0)
      Tau[k] =       tiltBoot_ml(tau, N, L, Ltau, B, Indices, weights,
				 B2_hi, subset_hi, W,
				 ytol, maxiter);
    else
      Tau[k] = 1.0 - tiltBoot_ml(tau, N, L, Ltau, B, Indices, weights,
				 B2_lo, subset, W,
				 ytol, maxiter);
  }
  if(debug>0) printf("}Exiting S_tiltBoot_ml\n");
  return;
}



void S_tiltBootSolve(S_LONG *Nobs, S_LONG *Nboot, S_LONG *Nstat, S_LONG *Nprob,
		     double *Observed, double *Replicates, double *Probs,
		     double *L, S_LONG *Indices,
		     S_LONG *Nweights, double *weights,
		     S_LONG *pnGroups, S_LONG pGroupSizes[],
		     double *xTol, double *yTol,
		     S_LONG *pDebug, S_LONG *MaxIter,
		     double *expTau, double *mlTau)
{
  /* Solve for tau, given desired probabilities.
     See documentation for tiltBoot_exp and tiltBoot_ml

     To be called from S, supports multiple parameters.

     col means are assumed to be swept out of L (haven't tested otherwise)

     On input the initial guesses for tau are in expTau
     Output is tau values, in expTau and mlTau
  */
  S_LONG N = *Nobs;
  S_LONG B = *Nboot;
  S_LONG P = *Nstat;
  S_LONG K = *Nprob;
  double xtol = *xTol;
  double ytol = *yTol;
  S_LONG one = 1;
  S_LONG g, gi, i, j, k, jk, imax, iter;
  S_LONG B2_hi, B2_lo, sgn;

  double Lrange[2], domain[2];
  double alpha, tau, Q, fraction, factor;
  double F0lo, F0hi;

  S_LONG    *subset, *subset_hi;
  double *Sj, *W, *Ltau = NULL, *Lgroup;

  if(*Nweights == 0) weights = NULL;
  debug = *pDebug;
  if(debug>0) printf("Entering S_tiltBootSolve{\n");
#ifdef WIN64 //S7_Change_LONG
  if(debug>3) printf("  N=%lld, B=%lld, P=%lld, K=%lld, NWeights=%lld, nGroups=%lld\n",
		     N, B, P, K, *Nweights, *pnGroups);
#else
  if(debug>3) printf("  N=%ld, B=%ld, P=%ld, K=%ld, NWeights=%ld, nGroups=%ld\n",
		     N, B, P, K, *Nweights, *pnGroups);
#endif //WIN64
  if(debug>3) printVectorDouble("  Observed", P, Observed);
  if(debug>4) printVectorDouble("  Replicates", B*P, Replicates);
  if(debug>3) printVectorDouble("  L", N, L);


  nGroups = *pnGroups;      /* Set global variables */
  groupSizes = pGroupSizes;

  subset = Salloc( B, S_LONG);
  Sj   = Salloc( B, double);
  W    = Salloc( N, double);
  if(nGroups > 1){
    Ltau = Salloc(N, double);
    V0   = Salloc(N, double);
    Lambda = Salloc(nGroups, double);
    Lgroup = Salloc(N, double);
    /* Initialize V0 to contain observation weights, sum to 1 in each group. */
    for(g = 0, gi = 0; g < nGroups; g++)
      for(i = 0; i < groupSizes[g]; i++)
	V0[ gi++ ] = 1.0 / (double) groupSizes[g];
  }

  imax = 0;
  for(j = 0; j < P; j++) {  /* j'th statistic */
    /* L and Replicates are matrices, L[N, P] and Replicates[B, P]
       We work with only the j'th column in this loop.
       So think of these as vectors in this loop.
       At end of loop increment these so they point to the next columns. */
#ifdef WIN64 //S7_Change_LONG
    if(debug > 1) printf(" Starting statistic %lld\n", j+1);
#else
    if(debug > 1) printf(" Starting statistic %ld\n", j+1);
#endif //WIN64

    makeSubsets(B, Replicates, Observed[j], subset, &B2_lo, &B2_hi);
    subset_hi = subset + (B - B2_hi);
#ifdef WIN64 //S7_change_LONG
    if(debug > 2) printf("  B=%lld, B2_lo=%lld, B2_hi=%lld\n", B, B2_lo, B2_hi);
#else
	if(debug > 2) printf("  B=%ld, B2_lo=%ld, B2_hi=%ld\n", B, B2_lo, B2_hi);
#endif //WIN64
    if(debug > 4) printVectorLong("  subset", B, subset);

    F0lo = (double)B2_hi/(double)B;
    F0hi = (double)B2_lo/(double)B;

    if(nGroups > 1){
      /* While doing exponential tilting, let Lgroup contain values of
	 L / (groupSizes/N), so can multiply by common tau instead of tau[g].*/
      for(g = 0, gi = 0; g < nGroups; g++){
	factor = (double) N / (double) groupSizes[g];
	for(i=0; i<groupSizes[g]; i++, gi++)
	  Lgroup[gi] = factor * L[gi];
      }
      if(debug>3) printVectorDouble("  Lgroup = L/(Ng/N)", N, Lgroup);
    }
    else  /* Let Lgroup point to L, to avoid if/else below. */
      Lgroup = L;
    S_bootstrapSums(&N, &one, &N, &B, Lgroup, Indices, Sj);

    if(debug > 4) printVectorDouble("  Sj", B, Sj);

    /* Find range of L and domain for tau */
    if(nGroups == 1)
      findDomain(N, L, Lrange, domain);
    else {
      domain[0] = -HUGE_VAL;
      domain[1] =  HUGE_VAL;
    }

    for(k = 0; k < K; k++) {  /* k'th probability */
#ifdef WIN64 //S7_Change_LONG
      if(debug > 1) printf(" Starting prob number %lld, alpha = %g\n",
			   k+1, Probs[k]);
#else
	  if(debug > 1) printf(" Starting prob number %ld, alpha = %g\n",
			   k+1, Probs[k]);
#endif //WIN64
      if(debug > 4) printVectorDouble("  V0 in S_tiltBootSolve(1)", N, V0);
      jk = j + P*k;

      alpha = Probs[k];
      if(alpha < F0lo && alpha > (1-F0hi)){
	sprintf(message, "solution is not uniquely-defined for alpha = %g, due to replications which equal the observed value", alpha);
	PROBLEM
	  message
	  WARNING(NULL_ENTRY);
      }


      /* Find tau for exponential tilting */
      tau = expTau[jk];  /* Initial guess */
      iter = *MaxIter;
      if(alpha < F0lo) { /* low probability; tau negative */
        sgn  = -1;
        tau  = solve_tiltBoot_exp(tau, my_qnorm(F0lo),
				  my_qnorm(alpha),
				  xtol, ytol, &iter, sgn,
				  N, Lgroup, B, Sj, weights,
				  B2_hi, subset_hi);
      }
      else { /* high probability; tau positive */
        sgn  = 1;
        tau  = solve_tiltBoot_exp(tau, my_qnorm(F0hi),
				  my_qnorm(1.0 - alpha),
				  xtol, ytol, &iter, sgn,
				  N, Lgroup, B, Sj, weights, B2_lo, subset);
      }

      if(debug > 1)
	printf("exp tilt initial guess=%g, final=%g\n", expTau[jk], tau);
      expTau[jk] = tau;
      if(iter > imax) imax = iter;

      /* Find tau for ML tilting. */
      if(nGroups == 1){ /* Transform exp tau for starting value. */
	/* If multiple groups, then do not transform for starting value. */
	Q = tau;
	S_tiltMean_exp(Nobs, &one, Nweights, L, weights, 
		       pnGroups, groupSizes, pDebug,
		       &Q); /* Now Q = tilted mean. */
	tau  = tau/(1.0 + tau*Q);

	/* Make sure that tau is within the interior of the domain */
  	fraction = (Q - Lrange[0]) / (Lrange[1] - Lrange[0]);
  	if(tau > fraction * domain[1])
  	  tau = fraction * domain[1];
  	if(tau < (1.0-fraction) * domain[0])
  	  tau = (1.0-fraction) * domain[0];
      }
      if(debug > 4) printVectorDouble("  V0 in S_tiltBootSolve(2)", N, V0);

      iter  = *MaxIter;
      if(alpha <= F0lo) { /* low probability; tau negative */
	mlTau[jk] =
	  solve_tiltBoot_ml(tau, F0lo,
			    alpha, xtol, ytol, &iter,
			    domain[0], domain[1], sgn,
			    N, L, Ltau, B, Indices, weights,
			    B2_hi, subset_hi, W);
      }
      else { /* high probability; tau positive */
	mlTau[jk] =
	  solve_tiltBoot_ml(tau, F0hi,
			    (1.0 - alpha), xtol, ytol, &iter,
			    domain[0], domain[1], sgn,
			    N, L, Ltau, B, Indices, weights,
			    B2_lo, subset, W);
      }
      if(debug > 1)
	printf("ml tilt initial guess=%g, final=%g\n", tau, mlTau[jk]);
      if(debug > 4) printVectorDouble("  V0 in S_tiltBootSolve(3)", N, V0);
      if(iter > imax) imax = iter;

    }  /* End of k loop */

    /* Make L and Replicates point to the next columns of their data */
    L += N;
    Replicates += B;
  }    /* End of j loop */
  if(debug>0) printf("}Exiting S_tiltBootSolve\n");
}


static double solve_tiltBoot_exp(double a, double F0,
				 double v,
				 double xtol, double ytol, S_LONG *iter,
				 S_LONG sgn,
				 S_LONG N, double *L,
				 S_LONG B, double *Sj, double *weights,
				 S_LONG B2, S_LONG *subset)
{
  /* Solve f(tau)=v;  F(tau)=f(tau)-v,
     f is qnorm(bootstrap tilting probs).
     a = initial guess
     F0 = value at 0
     v = goal, either qnorm(alpha) or qnorm(1-alpha).
     sgn = 1 if answer should be positive, else -1

     This should only be called with one of two subsets of the
     bootstrap samples, those for which
     replicate >= (or <=) observed.  For either subset,
     the bootstrap tilting probs -> 0 as abs(tau) increases,
     though not necessarily monotonely.
  */

  S_LONG   MaxIter;
  double Fa, b, Fb, c, Fc, S0, x;

  MaxIter = *iter;
  *iter  = 0;

  F0 = F0 - v;
  if(fabs(F0) < ytol) return(0.0);
  S0 = F0/fabs(F0);  /* This should be positive */

  if(a * sgn < 0.0){
    PROBLEM
      "initial guess has wrong sign in solve_tiltBoot_exp"
      WARNING(NULL_ENTRY);
    a = -a;
  }
  if(!a) a = .001 * sgn;

  Fa = my_qnorm(tiltBoot_exp(a, N, L, B, Sj, weights, B2, subset)) - v;
  *iter += 1;
  if(fabs(Fa) < ytol) return(a);

  if(S0*Fa < 0) {
    x  = bisect_tiltBoot_exp(a, 0.0, Fa, F0, v, xtol, ytol, iter,
			     N, L, B, Sj, weights, B2, subset);
    return x;
  }

  /* a and 0 don't bracket, and a should be closer to answer.
     Function could be non-monotone, so don't use secant method.
     Instead take wider intervals to bracket answer.
     Keep track of best answer (c, Fc) in case fail to converge.
  */
  if(fabs(Fa) < fabs(F0)){
    c = a;
    Fc = Fa;
  }
  else {
    c = 0.0;
    Fc = F0;
  }
#ifdef WIN64 //S7_Change_LONG
  if(debug > 3) printf("In solve_tiltBoot_exp, goal=%g, maxiter=%lld\n",
		       S_pnorm(v), *iter);
#else
  if(debug > 3) printf("In solve_tiltBoot_exp, goal=%g, maxiter=%ld\n",
		       S_pnorm(v), *iter);
#endif //WIN64
  while(*iter <= MaxIter){
    b = a;
    Fb = Fa;
    a = a * 2.0;
    Fa = my_qnorm(tiltBoot_exp(a, N, L, B, Sj, weights, B2, subset)) - v;
    *iter += 1;
    if(debug > 4) printf("In solve_tiltBoot_exp, tau=%g, prob(tau)=%g\n",
			 a, S_pnorm(Fa+v));
    if(fabs(Fa) < ytol) return(a);
    if(Fb * (Fa/fabs(Fa)) < 0) {
      x  = bisect_tiltBoot_exp(a, b, Fa, Fb, v, xtol, ytol, iter,
			       N, L, B, Sj, weights, B2, subset);
      return x;
    }
    if(fabs(Fa) < fabs(Fc)){
      c = a;
      Fc = Fa;
    }
  }
  /* printf("Warning, failed to converge in solve_tiltBoot_exp\n"); */
  /* printf("  Best answer is tau=%f, prob(tau)=%f\n",c,S_pnorm(Fc+v)); */
  sprintf(message, "failed to converge in solve_tiltBoot_exp. \n Best answer is tau = %f, prob(tau) = %f", c, S_pnorm(Fc+v));
  PROBLEM
    message
    WARNING(NULL_ENTRY);
  return(c);
}


static double solve_tiltBoot_ml(double a, double F0,
				double v,
				double xtol, double ytol, S_LONG *iter,
				double tauneg, double taupos,
				S_LONG sgn,
				S_LONG N, double *L, double *Ltau,
				S_LONG B, S_LONG *Indices, double *weights,
				S_LONG B2, S_LONG *subset,
				double *W)
     /* Solve f(tau)=v;  F(tau)=f(tau)-v, f is bootstrap tilting probs.
	a = initial guess
	F0 = value at 0
	v = goal
	sgn = 1 if answer should be positive, else -1
	W = work space for tiltBoot_ml
	iter = max iterations allowed on input, actual iterations on output.

	This should only be called with one of two subsets of the
	bootstrap samples, those for which
	replicate >= (or <=) observed.  For either subset,
	the bootstrap tilting probs -> 0 as abs(tau) increases,
	though not necessarily monotonely.
     */
{
  S_LONG   MaxIter = *iter;
  S_LONG   maxiter = *iter;  /* Number of iterations for finding lambda */
  double Fa, b, Fb, c, Fc, S0, x;

  *iter = 0;

  if(debug > 1)
    printf("Entering solve_tiltBoot_ml{, a= %g, F0= %g, v=%g\n",
	   a, F0, v);

  F0 = F0 - v;
  if(fabs(F0) < ytol) return(0.0);
  S0 = F0/fabs(F0);  /* This should be positive */

  if(a * sgn < 0.0){
    PROBLEM
      "initial guess has wrong sign in solve_tiltBoot_ml"
      WARNING(NULL_ENTRY);
    a = -a;
  }
  if(!a) a = .001 * sgn;
  if(debug > 4) printVectorDouble("  V0 in solve_tiltBoot_ml(1)", N, V0);
  Fa = tiltBoot_ml(a, N, L, Ltau, B, Indices, weights,
		   B2, subset, W,
		   ytol, maxiter) - v;
  *iter += 1;
  if(debug > 4) printVectorDouble("  V0 in solve_tiltBoot_ml(2)", N, V0);
  if(fabs(Fa) < ytol) return(a);

  if(S0*Fa < 0) {
    x  = bisect_tiltBoot_ml(a, 0.0, Fa, F0, v, xtol, ytol, maxiter, iter,
			    N, L, Ltau, B, Indices, weights,
			    B2, subset, W);
    if(debug > 4) printVectorDouble("  V0 in solve_tiltBoot_ml(3)", N, V0);
    return x;
  }

  /* a and 0 don't bracket, and a should be closer to answer.
     Function could be non-monotone, so don't use secant method.
     Instead take wider intervals to bracket answer.
     Keep track of best answer (c, F(c)) in case fail to converge.
  */
  if(fabs(Fa) < fabs(F0)){
    c = a;
    Fc = Fa;
  }
  else {
    c = 0.0;
    Fc = F0;
  }
  while(*iter <= MaxIter){
    b = a;
    Fb = Fa;
    a = a * 2.0;
    if(a <= tauneg) a = (tauneg + b)/2.0;
    if(a >= taupos) a = (taupos + b)/2.0;
    if(debug > 4) printVectorDouble("  V0 in solve_tiltBoot_ml(4)", N, V0);
    Fa = tiltBoot_ml(a, N, L, Ltau,
		     B, Indices, weights, B2, subset, W,
		     ytol, maxiter) - v;
    *iter += 1;
    if(fabs(Fa) < ytol) return(a);
    if(Fb * (Fa/fabs(Fa)) < 0) {
      x  = bisect_tiltBoot_ml(a, b, Fa, Fb, v, xtol, ytol, maxiter, iter,
			      N, L, Ltau, B, Indices, weights,
			      B2, subset, W);
      if(debug > 4) printVectorDouble("  V0 in solve_tiltBoot_ml(5)", N, V0);
      return x;
    }
    if(fabs(Fa) < fabs(Fc)){
      c = a;
      Fc = Fa;
    }
  }
  /* printf("Warning, failed to converge in solve_tiltBoot_ml\n"); */
  /* printf("  Best answer is tau=%f, prob(tau)=%f\n",c,Fc+v); */
  sprintf(message, "failed to converge in solve_tiltBoot_ml. \n Best answer is tau = %f, prob(tau) = %f", c, Fc+v);
  PROBLEM
    message
    WARNING(NULL_ENTRY);
  if(debug > 1) printf("}Exiting solve_tiltBoot_ml\n");
  return(c);
}



static double bisect_tiltBoot_exp(double a, double b, double Fa, double Fb,
				  double v,
				  double xtol, double ytol, S_LONG *iter,
				  S_LONG N, double *L,
				  S_LONG B, double *Sj, double *weights,
				  S_LONG B2, S_LONG *subset)
     /* Solve f(tau) = v, where f(tau) = qnorm(tiltBoot(tau)).
	a and b are endpoints of interval
	Fa and Fb must have opposite signs
	Fa = f(a) - v, Fb = f(b) - v

	Note that this is actually a bracketed secant method,
	modified to move the next point slightly toward the midpoint.

	iter is incremented by the number of function evaluations
     */
{
  double c, Fc, fraction;

  if(Fa*(Fb/fabs(Fb)) > 0){
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
      "internal error, roots don't bracket answer in bisect"
      WARNING(NULL_ENTRY);
  }
  xtol = xtol/2.0;

  if(Fa*(Fb/fabs(Fb)) > 0) {
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
      "internal error, roots don't bracket answer in bisect"
      WARNING(NULL_ENTRY);
  }
  do {
    /* c = (a + b)/2.0; */  /* bisect */
    fraction = Fb / (Fb-Fa);
    fraction = .5 + .99 * (fraction - .5);
    c = fraction * a + (1.0 - fraction) * b;

    Fc = my_qnorm(tiltBoot_exp(c, N, L, B, Sj, weights, B2, subset)) - v;
    *iter += 1;

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

static double bisect_tiltBoot_ml(double a, double b, double Fa, double Fb,
				 double v,
				 double xtol, double ytol,
				 S_LONG maxiter, /* for finding Lambda */
				 S_LONG *iter,
				 S_LONG N, double *L, double *Ltau,
				 S_LONG B, S_LONG *Indices, double *weights,
				 S_LONG B2, S_LONG *subset,
				 double *W)
     /* Solve f(tau) = v, where f(tau) = tiltBoot(tau).
	a and b are endpoints of interval
	Fa and Fb must have opposite signs
	Fa = f(a) - v, Fb = f(b) - v

	Note that this is actually a bracketed secant method,
	modified to move the next point slightly toward the midpoint.

	iter is incremented by the number of function evaluations
     */
{
  double c, Fc, fraction;

  if(debug > 1) printf("Entering bisect_tiltBoot_ml{\n");

  if(Fa*(Fb/fabs(Fb)) > 0) {
    /* printf("WARNING:  internal error,
       roots don't bracket answer in bisect\n");
    */
    PROBLEM
      "internal error, roots don't bracket answer in bisect"
      WARNING(NULL_ENTRY);
  }
  xtol = xtol/2.0;

  do {
    /* c = (a + b)/2.0; */  /* bisect */
    fraction = Fb / (Fb-Fa);
    fraction = .5 + .99 * (fraction - .5);
    c = fraction * a + (1.0 - fraction) * b;

    if(debug > 4) printVectorDouble("  V0 in bisect_tiltBoot_ml loop", N, V0);
    Fc = tiltBoot_ml(c, N, L, Ltau,
		     B, Indices, weights, B2, subset, W,
		     ytol, maxiter) - v;
    *iter += 1;

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
  if(debug > 1) printf("}Exiting bisect_tiltBoot_ml\n");
  return c;
}
