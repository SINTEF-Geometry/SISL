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
 * $Id: s6ratder.c,v 1.3 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6RATDER

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6ratder(double eder[],int idim,int ider,double gder[],int *jstat)
#else
void s6ratder(eder,idim,ider,gder,jstat)
     double eder[];
     int    idim;
     int    ider;
     double gder[];
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate the ider derivatives of a rational
*              point described in homogenous coordinates
*
* INPUT      : eder    - The derivatives in homogenous coordinates
*                        In sequence:
*                         Position (x,y,...h)
*                         1st der (x,y,...h)
*                         2nd der (x,y,...h)
*                         etc.
*              idim    - The dimension of the non homogenous space
*              ider    - The number of input derivatives
*
*
* OUTPUT     : jstat   - Status message
*                                        >0      : Warning
*                                        =0      : ok
*                                        <0      : Error
*              gder    - The derivatives in the nonhomogenous space
*
*
* METHOD     :  The curve P(u) can be written as the quotient
*               P(u) = T(u) / w(u) where T and w are ordinary splines.
*               The dimensions of T and w are idim and 1
*               respectively. The array eder contains position
*               and derivatives of the idim+1 dimensional curve
*               (T(u),w(u)).
*
*               Now, since wP = T, we find, by the Leibnitz formula,
*
*                 k
*                         k!     (k-i) (i)         (k)
*                sum   -------- w     P       =   T    .
*                      i!(k-i)!
*                i=0
*
*               Therefore
*
*
*                   --         k-1                      --
*             (k)   |   (k)             k!     (k-i) (i) |
*            P    = |  T    -  sum   -------- w     P    | / w .
*                   |                i!(k-i)!            |
*                   --         i=0                      --
*
*               This formula is applied recursively to evaluate P's derivatives.
*
*                                                          MF.
*
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-des-1988
* REVISED BY : Michael Floater, SI, 30/9/91 Removed division by t=1.
* REWRITTEN BY : Michael Floater, SI, 16/12/91. New algorithm.
* REWRITTEN BY : Michael Floater, SI, 25/8/92. Extend to arbitrary
*                   number of derivatives (by Leibnitz). Finally!
* REVISED BY : Paal Fugelli, SINTEF, 07/07-94. Added free'ing of binom and
*              initiation to SISL_NULL to avoid memory leakage.
*********************************************************************
*/
{
  int kpos=0;          /* Position of error.                     */
  double w0;           /* The denominator.                       */
  int ki;              /* Count through dimensions.              */
  int id;              /* Count through derivatives.             */
  int *binom = SISL_NULL;   /* Array for binomial coefficients.       */
  double sum;          /* Binomial (Leibnitz) expansion.         */
  int idimp1;          /* idim + 1.                              */
  int iw;              /* Pointer to a weight.                   */
  int igder;           /* Pointer to already calculated derivs.  */
  int i,j,k;           /* Counters.                              */
  int iwfix;           /* Initial value of iw in Leibnitz loop.  */

  if (ider<0) goto err178;
  if (idim<1) goto err102;

  idimp1 = idim + 1;

  /* Find denominator. */

  w0 = eder[idim];
  if (DEQUAL(w0,DZERO)) w0 = (double)1.0;

  /* Set up initial binomial coefficient (1). */

  binom = newarray(ider+1, INT);
  if(binom == SISL_NULL) goto err179;

  binom[0] = 1;


  /* Calculate position first. */

  for(ki=0; ki<idim; ki++)
  {
      gder[ki] = eder[ki] / w0;
  }



  /* Then derivatives if there are any. */

  for(id=1,j=idim,k=idimp1; id<=ider; id++,k++)
  {
      /* Calculate the new row of binomial coefficients. */

      binom[id] = 1;

      for(i=id-1; i>=1; i--)
      {
	  binom[i] += binom[i-1];
      }


      /* Run through the idim dimensions, calculating each
	 coefficient of the id'th derivative of
	 the rational curve (in gder). */

      iwfix = k + idim;

      for(ki=0; ki<idim; ki++,j++,k++)
      {
	  /* Calculate the Leibnitz sum (the binomial
	     coefficient in the first term is always 1). */

	  sum = eder[iwfix] * gder[ki];

          for(i=1,igder=idim+ki,iw=iwfix-idimp1;
	  i<id;
	  i++,igder+=idim,iw-=idimp1)
          {
	      sum += (double)binom[i] * eder[iw] * gder[igder];
          }

	  gder[j] = (eder[k] - sum) / w0;

      }

  }

  /* Done. */


  *jstat = 0;
  goto out;

/* idim less than 1. */
 err102: *jstat = -102;
         s6err("s6ratder",*jstat,kpos);
         goto out;

/* Derivative negative */
 err178: *jstat = -178;
         s6err("s6ratder",*jstat,kpos);
         goto out;


/* Not enough memory */
 err179: *jstat = -179;
         s6err("s6ratder",*jstat,kpos);
         goto out;


out:
  if (binom != SISL_NULL) freearray(binom);

  return;
}
