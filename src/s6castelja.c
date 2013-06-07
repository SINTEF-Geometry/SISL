/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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

