/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s6castelja.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6DECASTELJAU

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    s6deCasteljau(double C[], double a, double b, double t, int k, double D[],
		     int* jstat)
#else
void s6deCasteljau(C,a,b,t,k,D,jstat)
     double C[];
     double a,b,t;
     int k;
     double D[];
     int* jstat;
#endif
/*
***************************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To subdivide a 1-dim Bezier curve 'f' at a point 't'
*	       using the deCasteljau lgorithm, and to calculate 
*              the value f(t).
*              
*
*
*
* INPUT      : C[0:k-1]	- Bezier coefficients of f relative to the 
*			  intervall [a,b].
*              a	- start of parameter intervall.
*              b	- end of parameter intervall.
*              t 	- parameter value at which to subdivide.
*                         Assumes: a < b , and a <= t <= b.
*              k	- polynomial order (=degree+1) of f.
*              D[0:2*k-1] - allocated space.
*
*
* OUTPUT     : D[] 	- Bezier coefficients for the subdivided repr of 'f'.
*			  D[r],  r=0,...,k-1	-	coeffs. on [a,t],
*			  D[k+r], r=0,...,k-1	-	coeffs. on [t,b],
*                         D[k-1] = D[k]		-	value of f at t.
*
* METHOD     : 		- The multiaffine blossom  F of f at argument
*                         bags consisting of multipla of a,c, and t are 
*                         calculated using the deCasteljau algorithm. These
*                         values are stored in an local array A[0:k*k-1]
*                         as follows:
*			  A[k*r+j] = F(a,..,a,t,..,t,b,..,b)
*                             where a is repeated k-1-j times,
*                                   t is repeated r times,
*                                   b is repeated j-r times.
*			  In particular
*			  A[k*r+r] r=0,...,k-1 are Bezier coeffs. on [a,t],
*			  A[k*(k-1-r)+k-1] r=0,...,k-1 : coeffs. on [t,b],
*			  A[k*(k-1)+k-1] : value of f at t.
*                        

*              jstat     - status messages 
*                          < 0 - error
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Kyrre Strom, SI, 93-01.
*
*
****************************************************************************
*/
{
  int r,j,kk=k*k,kr;
  double alpha;
  double Al[16];
  double* A = SISL_NULL;


  *jstat = 1;
  if (a > b || DEQUAL(a,b) ) goto err109;

  if (k > 4 )
    {
      A = newarray(kk,double);
      if (A == SISL_NULL) goto err101;
    }
  else
    A = Al;

  for (j=0; j<k; j++)
    A[j] = C[j];

  alpha = (b-t)/(b-a);
  for (r = 1; r < k; r++)
    for (j = r; j < k; j++)
      A[k*r+j] = alpha*A[k*(r-1)+j-1] + (1-alpha)*A[k*(r-1)+j];
  
  for (kk--,kr=r=0; r<k; r++,kr+=k)
    {
      D[r] = A[kr+r];
      D[k+r] = A[kk-kr];
    }

  goto out;

 err109: *jstat = -109;
  goto out;

 err101: *jstat = -101;
  goto out;

 out: 
  if (A != SISL_NULL && A != Al)
    freearray(A);
  return ;
}

