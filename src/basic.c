#include <S.h>
#include <machine.h>

/* Functions in this file:
   S_bootstrapSums      bootstrap sums:     answer[i,j] = sum(x[indices[,i],j])
   S_bootstrapMeans     bootstrap means:    answer[i,j] = mean(x[indices[,i],j])
   S_bootstrapProducts  bootstrap products: answer[i,j] = prod(x[indices[,i],j])
   S_bootstrapProducts_log  bootstrap products, but use logs
   S_bootstrapVars      bootstrap variances

   These functions have a common structure.
   Arguments are
       N = number of observations in x
       P = number of columns      in x (often 1)
       M = size of each bootstrap sample (often N)
       B = number of bootstrap samples
       x[N,P],
       indices[M,B]  (matrix data in column order, as in S)
   Return
       y[B,P]
   where y[b,j] = f( x[indices[,b], j])  and f is sum, mean, prod, etc.
*/

void S_bootstrapSums(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		   double *x, S_LONG *indices,
		   double *answer)
{
  /* x[N,P], indices[M,B]  (matrix data in column order, as in S)
     Return y[B,P], where y[i,j] = sum( x[indices[,i],j])
  */
  S_LONG i, j, k, N, P, M, B;
  S_LONG ind;
  double *ans, sum;

  ans = answer;
  N =  *NN;
  P =  *PP;
  M =  *MM;
  B =  *BB;

  x--;  /* because indices are in range 1:N, not 0:(N-1) */

  for(k=0; k<P; k++){   /* For each of the P trials in x, */
    for(j=0; j<B; j++){ /* the B bootstrap samples        */
      sum = 0.0;
      for(i = 0; i < M; i++){
	ind = indices[i];
	if(ind > 0 && ind <= N)
	  sum += x[ind];
      }
      *(ans++)=sum;
      indices += M;  /* Point to the next column of indices */
    }
    x+=N;
    indices -= (M * B); /* Point back to the original start */
  }
  return;
}

void S_bootstrapMeans(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		   double *x, S_LONG *indices,
		   double *answer)
{
  /* x[N,P], indices[M,B]  (matrix data in column order, as in S)
     Return y[B,P], where y[i,j] = mean( x[indices[,i],j])
  */
  S_LONG i, j, k, N, P, M, B, sumi;
  S_LONG ind;
  double *ans, sum;

  ans = answer;
  N =  *NN;
  P =  *PP;
  M =  *MM;
  B =  *BB;

  x--;  /* because indices are in range 1:N, not 0:(N-1) */

  for(k=0; k<P; k++){   /* For each of the P trials in x, */
    for(j=0; j<B; j++){ /* the B bootstrap samples        */
      sum = 0.0;
      sumi = 0;
      for(i = 0; i < M; i++){
	ind = indices[i];
	if(ind > 0 && ind <= N){
	  sum += x[ind];
	  ++sumi;
	}
      }
      *(ans++) = sum / (double) sumi;
      indices += M;  /* Point to the next column of indices */
    }
    x+=N;
    indices -= (M * B); /* Point back to the original start */
  }
  return;
}

void S_bootstrapProducts(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		   double *x, S_LONG *indices,
		   double *answer)
{
  /* x[N,P], indices[M,B]  (matrix data in column order, as in S)
     Return y[B,P], where y[i,j] = product( x[indices[,i],j])
  */
  S_LONG i, j, k, N, P, M, B;
  S_LONG ind;
  double *ans, product;

  ans = answer;
  N =  *NN;
  P =  *PP;
  M =  *MM;
  B =  *BB;

  x--;  /* because indices are in range 1:N, not 0:(N-1) */

  for(k=0; k<P; k++){   /* For each of the P trials in x, */
    for(j=0; j<B; j++){ /* the B bootstrap samples        */
      product = 1.0;
      for(i = 0; i < M; i++){
	ind = indices[i];
	if(ind > 0 && ind <= N)
	  product *= x[ind];
      }
      *(ans++) = product;
      indices += M;  /* Point to the next column of indices */
    }
    x+=N;
    indices -= (M * B); /* Point back to the original start */
  }
  return;
}


void S_bootstrapProducts_log(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		   double *x, S_LONG *indices,
		   double *answer)
{
  /* x[N,P], indices[N,B]  (matrix data in column order, as in S)
     Return y[B,P], where y[i,j] = product( x[indices[,i],j])
     Use logs; x must be positive.
  */

  S_LONG i, j, k, N, P, M, B;
  S_LONG ind;
  double *ans, sum;

  ans = answer;
  N =  *NN;
  P =  *PP;
  M =  *MM;
  B =  *BB;

  /* Use logs for calculations */

  k=N*P;
  for(i=0; i<k; i++)
    x[i] = (x[i] <= 0) ? -HUGE_VAL : log( x[i] );

  x--;  /* because indices are in range 1:N, not 0:(N-1) */

  for(k=0; k<P; k++){   /* For each of the P trials in x, */
    for(j=0; j<B; j++){ /* the B bootstrap samples        */
      sum = 0.0;
      for(i = 0; i < M; i++){
	ind = indices[i];
	if(ind > 0 && ind <= N)
	  sum += x[ind];
      }
      *(ans++) = exp(sum);
      indices += M;  /* Point to the next column of indices */
    }
    x+=N;
    indices -= (M * B); /* Point back to the original start */
  }
  return;
}


void S_bootstrapVars(S_LONG *NN, S_LONG *PP, S_LONG *MM, S_LONG *BB,
		   double *x, S_LONG *indices,
		   double *answer)
{
  /* x[N,P], indices[N,B]  (matrix data in column order, as in S)
     Return y[B,P], where y[i,j] = var( x[indices[,i],j])
     Use 1/(n-1) divisor.
  */
  S_LONG i, j, k, N, P, M, B, sumi;
  S_LONG ind;
  double *ans, sum, mean, temp;

  ans = answer;
  N =  *NN;
  P =  *PP;
  M =  *MM;
  B =  *BB;

  x--;  /* because indices are in range 1:N, not 0:(N-1) */

  for(k=0; k<P; k++){   /* For each of the P trials in x, */
    for(j=0; j<B; j++){ /* the B bootstrap samples        */
      sum = 0.0;
      sumi = 0;
      for(i = 0; i < M; i++){
	ind = indices[i];
	if(ind > 0 && ind <= N){
	  sum += x[ind];
	  ++sumi;
	}
      }
      mean = sum / (double) sumi;
      sum = 0.0;
      if(sumi == M) /* Don't need to check for indices out of range */
	for(i = 0; i < M; i++) {
	  temp = x[indices[i]] - mean;
	  sum += temp * temp;
	}
      else
	for(i = 0; i < M; i++){
	  ind = indices[i];
	  if(ind > 0 && ind <= N){
	    temp = x[indices[i]] - mean;
	    sum += temp * temp;
	  }
	}
      *(ans++)=sum/(sumi - 1.0);
      indices += M;  /* Point to the next column of indices */
    }
    x += N;
    indices -= (M * B); /* Point back to the original start */
  }
  return;
}
