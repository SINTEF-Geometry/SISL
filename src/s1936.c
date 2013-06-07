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
 * $Id: s1936.c,v 1.4 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1936

#include "sislP.h"
#define MAX_SIZE  50



#if defined(SISLNEEDPROTOTYPES)
void
s1936 (SISLCurve * crv, double etd[], int ind, double *curvd, int *jstat)
#else
void
s1936 (crv, etd, ind, curvd, jstat)
     SISLCurve *crv;
     double etd[];
     int ind;
     double *curvd;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To express a B-spline curve using a refined basis.
*	   of the refinement etref of the array et.
*
*
* INPUT:   crv	- The original curve.
*	   etd	- The knots of the subdivided curve.
*	   ind	- The number of vertices of the subdivided curve.
*
*
* OUTPUT:  curvd - Array containing the vertices of the refined curve.
*          jstat - Output status:
*                   < 0: Error.
*                   = 0: Ok.
*                   > o: Warning.
*
* METHOD:  The vertices are calculated using the "Oslo"-algorithm
*	   developped by Cohen, Lyche and Riesenfeld.
*
*
* REFERENCES: Cohen, Lyche, Riesenfeld: Discrete B-splines and subdivision
*	      techniques in computer aided geometric design, computer
*             graphics and image processing, vol 14, no.2 (1980)
*
* CALLS: s1937, s6err.
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07.
* CORRECTED BY : Paal Fugelli, SINTEF, Oslo 1994-07. Added test for
*                equality using DEQUAL in knot verification loop.
*
*********************************************************************
*/
{
  int kpos = 0;			/* Error position indicator	*/
  int ki, kj, kr, low;		/* Loop control parameters 	*/
  int nu;			/* Pointer into knot vector:
				 * knt[nu-1]<=etd[kj]<knt[nu]	*/
  int ins;			/* Number of vertices of curve	*/
  int iordr;			/* Order of curve		*/
  int idim;			/* Dimension of space where the
				 * curve lies			*/
  double *knt = SISL_NULL;		/* Original knot-vector. */
  double *coef = SISL_NULL;		/* Original coefficient arrays.	*/
  double sum;			/* Used to compute vertices of
				 * new curve			*/
  double sarray[MAX_SIZE];
  int alloc_needed=FALSE;
  double *alfa = SISL_NULL;		/* Array needed in subroutine
				 * s1937 (Oslo-algorithm)	*/

  *jstat = 0;


  /* Initialization. */

  knt = crv->et;
  ins = crv->in;
  iordr = crv->ik;
  idim = crv->idim;
  coef = crv->ecoef;


  /* Test if legal input. */

  if (iordr < 1)
    goto err110;
  if (ins < iordr || ind < iordr)
    goto err111;
  if (idim < 1)
    goto err102;


  /* Allocate array for internal use only. */

  if (iordr > MAX_SIZE)
    {
       if ((alfa = newarray (iordr, DOUBLE)) == SISL_NULL)
	 goto err101;
       alloc_needed = TRUE;
    }
  else
    alfa = sarray;

  /* Find if etd is a refinement of the original knot vector knt. */

  kj = 0;

  for (ki = 0; kj < ins; ki++)
    {
      if (ki >= ind)
	goto err111;
      if ( DEQUAL(knt[kj], etd[ki]) )  /* PFU 25/07-1994 */
      {
	kj++;
	continue;
      }
      if (knt[kj] > etd[ki])
	continue;
      if (knt[kj] < etd[ki])
	goto err112;
    }

  /* etd is a refinement of original knot vector knt
   * Produce refined curve. */

  nu = 1;
  for (kj = 0; kj < ind; kj++)
    {

      /* We want to find  knt[nu-1] <= etd[kj] < knt[nu]
         The copying of knots guarantees the nu-value to be found.
         Since kj is increasing, the nu-values will be increasing
         due to copying of knots. */

       /* for (; (((knt[nu - 1] > etd[kj]) || (etd[kj] >= knt[nu])) && (nu != ins));
	   nu++) ; */
       for (; (((knt[nu - 1] > (double)0.5*(etd[kj]+etd[kj+1])) 
		 || ((double)0.5*(etd[kj]+etd[kj+1]) >= knt[nu])) && (nu != ins)); nu++) ; 

      /* Now we have  knt[nu-1] <= etd[kj] < knt[nu],
         so the discrete B-splines can be calculated. */

      s1937 (knt, iordr, kj + 1, nu, alfa, etd);


      /* Compute the coefficients of etd. */

      low = nu - iordr;
      for (ki = 0; ki < idim; ki++)
	{
	  sum = (double) 0.0;
	  for (kr = MAX (0, low); kr < nu; kr++)
	    sum += alfa[kr - low] * coef[kr * idim + ki];
	  curvd[kj * idim + ki] = sum;
	}
    }

  /* OK. */

  goto out;


  /* Memory error */

err101:
  *jstat = -101;
  s6err ("s1936", *jstat, kpos);
  goto out;

  /* Error in B-spline curve description:
     Dimension less than 1. */

err102:
  *jstat = -102;
  s6err ("s1936", *jstat, kpos);
  goto out;

  /* Order less than 1. */

err110:
  *jstat = -110;
  s6err ("s1936", *jstat, kpos);
  goto out;

  /* No. of vertices less than order. */

err111:
  *jstat = -111;
  s6err ("s1936", *jstat, kpos);
  goto out;

  /* Error in knot-vector. */

err112:
  *jstat = -112;
  s6err ("s1936", *jstat, kpos);

out:
  if (alloc_needed)
    freearray (alfa);

  return;
}
#undef MAX_SIZE
