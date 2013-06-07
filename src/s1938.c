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
 * $Id: s1938.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1938

#include "sislP.h"
#define MAX_SIZE 50


#if defined(SISLNEEDPROTOTYPES)
void
s1938 (SISLSurf * srf, double etr1[], int inr1, double etr2[], int inr2,
       double **surfr, int *jstat)
#else
void
s1938 (srf, etr1, inr1, etr2, inr2, surfr, jstat)
     SISLSurf *srf;
     double etr1[];
     int inr1;
     double etr2[];
     int inr2;
     double **surfr;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To express a B-spline surface curve using a refined basis.
*
*
* INPUT:   srf	- The original surface.
*	   etr1	- The refined knot vector in the first parameter
*		  direction.
*	   inr1	- The number of vertices of the subdivided surface
*		  in first parameter direction.
*	   etr2	- The refined knot vector in the second parameter
*		  direction.
*	   inr2	- The number of vertices of the subdivided surface
*		  in the second parameter direction.
*
* OUTPUT:  surfr - Array containing the vertices of the refined
*		   surface.
*          jstat - Status variable.
*                   < 0: Error.
*	            = 0: Ok.
*                   > 0: Warning.
*
* METHOD:  The vertices are calculated using the "Oslo"-algorithm
*	   developped by Cohen, Lyche and Riesenfeld.
*
*
* REFERENCES: Cohen, Lyche, Riesenfeld: Discrete B-splines and subdivision
*	      techniques in computer aided geometric design, computer
*	      graphics and image processing, vol 14, no.2 (1980)
*
* CALLS: s1937, s6err.
*
*
* WRITTEN BY:  Christophe R. Birkeland, SI, 1991-08
* REVISED BY:  Johannes Kaasa, SI, May 1992 (Introduced NURBS)
*
*********************************************************************
*/
{
  int kpos = 0;
  int ki, kj, kl, kr, low;	/* Loop control parameters 		*/
  int start;
  int rem;			/* Used to store array index		*/
  int nu;			/* Pointer into knot vector:
				   knt[nu-1]<=etd[kj]<knt[nu]		*/
  int ins1;			/* Number of vertices of curve in first
				   parameter direction			*/
  int iordr1;			/* Order of curve in first
				   parameter direction			*/
  int ins2;			/* Number of vertices of curve in second
				   parameter direction			*/
  int iordr2;			/* Order of curve in second
				   parameter direction			*/
  double *knt1 = SISL_NULL;		/* Original knot-vector of surface.     */
  double *knt2 = SISL_NULL;		/* Original knot vector of surface.	*/
  double *coef = SISL_NULL;		/* Pointer to array of coefficients of
				   the surface				*/
  int idim;			/* Dimension of space where the
				   curve lies				*/
  double sum;			/* Used to store vertices of
				   new curve				*/
  double sarray[MAX_SIZE];
  int alloc_needed=FALSE;
  double *alfa = SISL_NULL;		/* Array needed in subroutine
				   s1937 (Oslo-algorithm)		*/
  double *ktsurf = SISL_NULL;	/* Array for internal use only		*/

  *jstat = 0;


  /* Initialization. */

  knt1 = srf->et1;
  knt2 = srf->et2;
  ins1 = srf->in1;
  ins2 = srf->in2;
  iordr1 = srf->ik1;
  iordr2 = srf->ik2;
  if (srf->ikind == 2 || srf->ikind == 4)
    {
      idim = srf->idim + 1;
      coef = srf->rcoef;
    }
  else
    {
      idim = srf->idim;
      coef = srf->ecoef;
    }
     
  /* Test if legal input. */

  if (iordr1 < 1 || iordr2 < 1)
    goto err115;
  if (ins1 < iordr1 || ins2 < iordr2)
    goto err116;
  if (idim < 1)
    goto err102;


  /* Allocate array for internal use only. */

  if (MAX(iordr1,iordr2) > MAX_SIZE)
    {
      if ((alfa = newarray(MAX(iordr1,iordr2),DOUBLE)) == SISL_NULL)
	goto err101;
      alloc_needed = TRUE;
    }
  else
    alfa = sarray;
  
  ktsurf = newarray (inr1 * inr2 * idim, DOUBLE);
  if (ktsurf == SISL_NULL)
    goto err101;


  /* Allocate array surfr for output. */

  *surfr = newarray (inr1 * inr2 * idim, DOUBLE);
  if (*surfr == SISL_NULL)
    goto err101;


  /* Find if etr1 is a refinement of the original
     knot vector knt1 (srf->et1). */

  kj = iordr1 - 1;

  for (ki = 0; kj < ins1; ki++)
    {
      if (ki >= inr1)
	goto err116;
      if (knt1[kj] > etr1[ki])
	continue;
      if (knt1[kj] < etr1[ki])
	goto err117;
      kj++;
    }

  /* etr1 is a refinement of original knot vector knt1
   * Produce surface refined in one direction. */

  nu = 1;
  for (kj = 0; kj < inr1; kj++)
    {
      /* We want to find  knt1[nu-1] <= etr1[kj] < knt1[nu]
	 The copying of knots guarantees the nu-value to be found.
	 Since kj is increasing, the nu-values will be increasing
	 due to copying of knots. */

      for (; (((knt1[nu - 1] > etr1[kj]) || (etr1[kj] >= knt1[nu]))
	      && (nu != ins1)); nu++) ;


      /* Now we have  knt1[nu-1] <= etr1[kj] < knt1[nu],
	 so the discrete B-splines can be calculated. */

      s1937 (knt1, iordr1, kj + 1, nu, alfa, etr1);


      /* Compute the temporary surface. */

      low = nu - iordr1 + 1;
      for (ki = 0; ki < ins2; ki++)
	{
	  rem = idim * kj + idim * inr1 * ki;

	  for (kl = 0; kl < idim; kl++)
	    {
	      sum = (double) 0.0;
	      start = nu - iordr1 + 1;
	      if (start < 1)
		start = 1;

	      for (kr = start; kr <= nu; kr++)
		sum += alfa[kr - low] * coef[ki * ins1 * idim + (kr - 1) * idim + kl];
	      ktsurf[rem + kl] = sum;
	    }
	}
    }

  /* Find if etr2 is a refinement of the original
   * knot vector knt2 (srf->et2). */

  kj = iordr2 - 1;

  for (ki = 0; kj < ins2; ki++)
    {
      if (ki >= inr2)
	goto err116;
      if (knt2[kj] > etr2[ki])
	continue;
      if (knt2[kj] < etr2[ki])
	goto err117;
      kj++;
    }

  /* etr2 is a refinement of original knot vector knt2
   * Produce surface refined in one direction. */

  nu = 1;
  for (ki = 0; ki < inr2; ki++)
    {
      /* We want to find  knt2[nu-1] <= etr2[ki] < knt2[nu]
	 The copying of knots guarantees the nu-value to be found.
	 Since kj is increasing, the nu-values will be increasing
	 due to copying of knots. */

      for (; (((knt2[nu - 1] > etr2[ki]) || (etr2[ki] >= knt2[nu]))
	      && (nu != ins2)); nu++) ;


      /* Now we have  knt2[nu-1] <= etr2[kj] < knt2[nu],
	 so the discrete B-splines can be calculated */

      s1937 (knt2, iordr2, ki + 1, nu, alfa, etr2);


      /* Compute the temporary surface. */

      low = nu - iordr2 + 1;
      for (kj = 0; kj < inr1; kj++)
	for (kl = 0; kl < idim; kl++)
	  {
	    sum = (double) 0.0;
	    start = nu - iordr2 + 1;
	    if (start < 1)
	      start = 1;
	    for (kr = start; kr <= nu; kr++)
	      sum += alfa[kr - low] * ktsurf[kj * idim + (kr - 1) * idim * inr1 + kl];
	    (*surfr)[ki * idim * inr1 + kj * idim + kl] = sum;
	  }
    }


  /* OK. */

  goto out;


  /* Memory error */

err101:
  *jstat = -101;
  s6err ("s1938", *jstat, kpos);
  goto out;

  /* Error in B-spline surface description:

     Dimension less than 1. */


err102:
  *jstat = -102;
  s6err ("s1938", *jstat, kpos);
  goto out;

  /* Order less than 1. */

err115:
  *jstat = -115;
  s6err ("s1938", *jstat, kpos);
  goto out;

  /* No. of vertices less than order. */

err116:
  *jstat = -116;
  s6err ("s1938", *jstat, kpos);
  goto out;

  /* Error in knot vector */

err117:
  *jstat = -117;
  s6err ("s1938", *jstat, kpos);
  goto out;

out:
  if (alloc_needed)
    freearray (alfa);
  if (ktsurf != SISL_NULL)
    freearray (ktsurf);

  return;
}
#undef MAX_SIZE
