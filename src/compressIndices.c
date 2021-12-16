#include <S.h>
#include <machine.h>

/* Entry points:
     S_compressIndices
     S_uncompressIndices

   There are no other functions or subroutines in this file.
*/


/* Input is an [m,B] matrix with maximum value n (column order, as in S).
   Store the number of times each of the n possibilities occur,
   in a shorter vector, of length K.
   Output is an [K,B] matrix.

   If an original vector is 5 4 5 2 2 5 with max=6,
   the frequences are 0 2 0 1 3 0
   and we'll store this as:    |++||+|+++||
   where | indicates borders and number of repetitions of + indicate frequency.
   This takes (m+n) bits; K = ceiling((m+n)/32).
*/
void S_compressIndices(S_LONG originalDim[], S_LONG x[], S_LONG *pn,
		       S_LONG *pK, S_ULONG y[])
{
  S_LONG m, B, n, K;           /* Input integers */
  S_LONG *xj, xij;             /* Pointer to parts of input data, & single val */
  S_LONG i, j, mn, ibit;       /* Working integers */
  S_LONG *frequency;           /* Working vectors */
  S_ULONG *yj, yij, bits[32]; /* ptr to part of output vector, work vec */

  m = originalDim[0];
  B = originalDim[1];
  n = *pn;
  K = *pK;
  bits[0] = 1;
  for(i=1; i<32; i++){
    bits[i] = bits[i-1] * 2;
  }

  mn = m+n;
  if(mn > 32*K)
    PROBLEM
      "Not enough space in output matrix; need k >= ceiling((m+n)/32)"
      RECOVER(NULL_ENTRY);

  frequency = Salloc(n+1, S_LONG);
  /* one extra spot, in case some input indices are zero.  Legal are 0...n */
  for(i=0; i <= n; i++) frequency[i] = 0;

  for(j=0; j<B; j++){  /* Compress the j'th vector. */
    xj = x + m * j;  /* pointer to start of this vector */
    yj = y + K * j;  /* pointer to start of this vector */

    /* Calculate frequencies */
    for(i=0; i < m; i++){
      xij = xj[i];
      if(xij > n)
	PROBLEM "Index value greater than maximum" RECOVER(NULL_ENTRY);
      if(xij < 0)
	PROBLEM "Index value less than zero" RECOVER(NULL_ENTRY);
      ++frequency[xij];
    }

    /* Convert the frequencies (of 1:m, not 0:m) to zeros and ones, 
       with frequencies encoded by lengths of strings of 1's between 0's.
       Store 32 such bits in each integer. */
    ibit = 0; /* This will run from 0 to 31, index to current bit. */
    yij = 0;  /* Temporary storage, for 32 bits */
    for(i=1; i <= n; i++){
      while(frequency[i] > 0){
	yij += bits[ibit++];
	--frequency[i];
	if(ibit == 32){  /* Store this number, start with bit 0 on the next. */
	  *yj++ = yij;
	  yij = 0;
	  ibit = 0;
	}
      }
      ibit++; /* leave a zero after a string of 1's */
      if(ibit == 32){
	*yj++ = yij;
	yij = 0;
	ibit = 0;
      }
    }
    if(ibit)
      *yj++ = yij;
  }
  return;
}


/* Convert from compressed form into ordinary indices.
   These will be sorted, instead of matching the original.

   Input y is a [K,B] matrix.
   Output x is an [m,B] matrix.
*/
void S_uncompressIndices(S_LONG originalDim[], S_ULONG y[], S_LONG *pn,
			 S_LONG *pK, S_LONG x[])
{
  S_LONG m, B, n, K;           /* Input integers */
  S_LONG *xj;                  /* Pointer to part of output vector */
  S_LONG i, i2, j;             /* Working integers */
  S_LONG *zeroone;             /* Working vector */
  S_LONG *pzo;                 /* Pointer to parts of the zeroone vector. */
  S_ULONG *yj, yij;    /* pointer to part of input, and single value */

  m = originalDim[0];
  B = originalDim[1];
  n = *pn;
  K = *pK;

  zeroone = Salloc(32*K, S_LONG);

  for(j=0; j<B; j++){  /* Uncompress the j'th vector. */
    xj = x + m * j;  /* pointer to start of this vector */
    yj = y + K * j;  /* pointer to start of this vector */

    pzo = zeroone;
    for(i=0; i<K; i++){
      yij = yj[i];
      /* Uncompress the bit patterns in yj into 32 zeros and ones. */
      for(i2=0; i2<32; i2++){
	/*	*pzo++ = yij & 01;*/
	/* yij >> 1; */
	*pzo++ = yij % 2;
	yij /= 2;
      }
    }
    /* Use the zeros and ones to fill in the indices matrix x. */
    pzo = zeroone;
    for(i=1; i<=n; i++){
      while(*pzo++)
	*xj++ = i;
    }
  }
  return;
}
