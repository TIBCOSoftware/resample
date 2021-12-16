
#include <S.h>
#if defined(SPLUS_VERSION) && SPLUS_VERSION >= 4000
#  include <newredef.h>
#endif

/*
 * The main entry point is S_pTestSums
 *    It in turn calls permTestExchange
 */



/* computeSumsSmaller:

   Compute the sums of the `nsmaller' rows of `xcols' indexed by `inds'.  
   Put result in the `sumsrow'th row of `sumscols'.  
   
   Assumes `sumscols' is already zeroed out. `p' = number of columns of 
   `xcols' and `sumscols'. 
*/
static void
computeSumsSmaller(double **xcols, 
		   S_LONG p, 
		   S_LONG *inds, 
		   S_LONG nsmaller, 
		   double **sumscols,
		   S_LONG sumsrow)
{
  S_LONG i, j;

  for(i = 0; i < nsmaller; i++)
    for(j = 0; j < p; j++)
      sumscols[j][sumsrow] += xcols[j][inds[i]];

  return;
}

/* computeSumsLarger:

   Compute the sums of the n-`nsmaller' rows of `xcols' not indexed by `inds'
   by subtracting the sum of the `nsmaller' rows from the total sum.
*/
static void
computeSumsLarger(double **xcols, 
		  S_LONG p, 
		  S_LONG *inds, 
		  S_LONG nsmaller, 
		  double *sumx, 
		  double **sumscols,
		  S_LONG sumsrow)
{
  S_LONG j;

  computeSumsSmaller(xcols, p, inds, nsmaller, sumscols, sumsrow);
  for(j = 0; j < p; j++)
    sumscols[j][sumsrow] = sumx[j] - sumscols[j][sumsrow];
  return;
}



/* samplePermuteNoSeed

   Perform a partial permutation of x, in place.
   First k elements of result are drawn randomly from x.
   Elements originally in those locations are now elsewhere in x. 
*/
static void
samplePermuteNoSeed(S_LONG n,
		    S_LONG *x,
		    S_LONG k)
{
  S_LONG j, i, xtemp;

  if(n == 1)
    return;
  for(i = 0 ; i < k; i++){
    j = i + (S_LONG) floor((n - i) * unif_rand(S_evaluator));
    /* random integer in [i, n - 1]  */
    xtemp = x[i];
    x[i] = x[j];
    x[j] = xtemp;
  }
  return;
}

/*
  permTestExchange
*/
static void
permTestExchange(S_LONG B,
		 S_LONG n, 
		 S_LONG p, 
		 S_LONG nsmaller,
		 double **xcols,
		 double *sumx,
		 int larger,
		 S_LONG *inds,
		 double **sumscols)
{
  S_LONG i;

  /* Note: last argument to computeSumsxxx is row number offset for 
     sumscols.  First row of sumscols are observed values, so skip it.
  */
  if(larger){
    for(i = 1; i < B + 1; i++) {
      samplePermuteNoSeed(n, inds, nsmaller);
      computeSumsLarger(xcols, p, inds, nsmaller, sumx, sumscols, i);
    }
  }
  else{
    for(i = 1; i < B + 1; i++) {
      samplePermuteNoSeed(n, inds, nsmaller);
      computeSumsSmaller(xcols, p, inds, nsmaller, sumscols, i);
    }
  }
  return;
}


/*
  S_pTestSums

  Main entry point, called from S

 */
void 
S_pTestSums(double *x, 
	    S_LONG *pn, 
	    S_LONG *pp, 
	    S_LONG *set1,
	    S_LONG *pB, 
	    double *sums)
{
  S_LONG i, j;
  S_LONG n = *pn;  /* number of rows of x */
  S_LONG p = *pp;  /* number of cols of x */
  S_LONG B = *pB;  /* number of permutations = number of rows of sums - 1 */
  S_LONG n0, n1;   
  S_LONG nsmaller;
  int larger;
  double *sumx, **xcols, **sumscols;

  /* space for overall sum */
  sumx = Salloc(p, double);  /* this zeros out sumx */

  /* arrays of pointers to the columns of x, sums */
  xcols = Salloc(p, double *);  
  sumscols = Salloc(p, double *);  

  for(j = 0; j < p; j++){
    xcols[j] = x + j*n;              /* x has n rows */
    sumscols[j] = sums + j*(B+1);    /* sums has B+1 rows */
  }

  /* Compute overall sum and the number in and sum for set1 == T */
  n1 = 0;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++)
      sumx[j] += xcols[j][i];  
    if(set1[i]){
      n1++;
      for(j = 0; j < p; j++)
	sumscols[j][0] += xcols[j][i]; /* first row of sums contains observed */
    }
    set1[i] = i;  /* turn set1 into an array of indices to use later */
  }

  n0 = n - n1;
  larger = (n1 > n0);  /* Is set1 == T the larger group?  */
  nsmaller = (larger) ? n0 : n1;

  seed_in((S_LONG *)NULL, S_evaluator); /* needed for unif_rand */

  permTestExchange(B, n, p, nsmaller, xcols, sumx, larger, set1, sumscols);

  seed_out((S_LONG *)NULL, S_evaluator);

  return;
}



