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
 * $Id: s1894.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1894

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1894 (double oknots[], int oik, int oin, int der1, int der2, double earray[],
       int dimp1, int narr, double *nknots[], int *nik, int *nin, int *jstat)
#else
void
s1894 (oknots, oik, oin, der1, der2, earray, dimp1, narr, nknots,
       nik, nin, jstat)
     double oknots[];
     int oik;
     int oin;
     int der1;
     int der2;
     double earray[];
     int dimp1;
     int narr;
     double *nknots[];
     int *nik;
     int *nin;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :  To produce a knot vector for the dot product of two derivates
*		of a B-spline curve having the input knot vector.
*
* INPUT      :  oknots	- The original knot vector.
*		oik	- The order of the original knot vector.
*		oin	- The number of degrees of freedom in the original
*			  knot vector.
*		der1	- The first of the derivatives.
*		der2	- The secondof the derivatives.
*		earray	- Description of the array to be used.
*		dimp1	- The dimension of the first matrix plane.
*		narr	- Number of parallel matrix planes.
*
* OUTPUT     :  nknots	- The new knot vector.
*		nik	- The order of the new knot vector.
*		nin	- The number of degrees of freedom in the new
*			  knot vector.
*		jstat    - Status variable:
*                                               > 0     : warning
*                                               = 0     : ok
*                                               < 0     : error
*
* METHOD     :  The order nk of the new new basis is determined.
*		The multiplicity of the knots are counted and the right number
*		of knots inserted according to the order of the original basis
*		and the derivates involved. At et[ik] and et[in+1] nk knots
*		are inserted.
*
* REFERENCES :  Fortran version:
*		Tor Dokken, SI, 1982-10
*
* CALLS      : s6err.
*
* WRITTEN BY :  Trond Vidar Stensby, SI, 1991-06
* REVISED BY :  Johannes Kaasa, SI, May 92. (Taking the number of derivatives
*               into account when calculating the multiplicity of the interior
*               knots, to reduce the continuity in accordance with the
*               derivatives. I also changed minimum allowed order from 1 to 2,
*               to avoid errors in s1890).
*
**************************************************************** */
{
  int size;			/* The total size of earray. */
  int mult;			/* Multiplicity of knots */
  int numb;			/* Number of new knots. */
  int kdim;			/* dimp1 -1  (sub-matrix dimension) */
  int empty;			/* Used to check if sub-matrix of earray
				   is zero. */
  int kl;			/* Loop control varibles. */
  int count1;
  int count2;
  int count3;
  int start;
  int stop;

  double eps;			/* Resolution. */
  double maximum;		/* The maximum value in et. */
  double prev;			/* Knot value. (extracted from orig) */
  double curr;			/* Knot value. (extracted from orig) */
  int kpos = 0;
  int der = max(der1, der2);

  *jstat = 0;
  

  /* Test if legal input. */

  if (oik <= 1 || oin < oik)
    goto err112;


  /* Test if knot vector degenerate. */

  if (oknots[oik - 1] >= oknots[oin])
    goto err112;


  /* The maximal number of knots to be produced at a specified knot value
   * is the order of the B-spline basis produced. */

  /* Allocate space for new knot vector */

  (*nknots) = newarray ((oin + oik) * oik, DOUBLE);
  if (*nknots == SISL_NULL)
    goto err101;


  /* Check if sub-matrix is zero. */

  kdim = dimp1 - 1;
  size = dimp1 * dimp1;
  empty = TRUE;

  for (count1 = 0; count1 < narr && empty; count1++)
    for (count2 = 0; count2 < kdim && empty; count2++)
      for (count3 = 0; count3 < kdim && empty; count3++)
	if (earray[count1 * size + count2 * dimp1 + count3] != (double)0.)
	  empty = FALSE;


  /* Assign value to nk. */

  if (empty)
    (*nik) = oik - min (der1, der2);
  else
    (*nik) = 2 * oik - der1 - der2 - 1;
  if ((*nik) < 2)
    (*nik) = 2;
  *nin = 0;


  /* Make resolution to be used for testing of knot value equalness. */

  eps = fabs (oknots[oin] - oknots[oik - 1]) * 1.0e-11;


  /* Production of knots. Initiate for calculation of knots.
     Find first knot not equal to start of curve. */

  maximum = oknots[oin];
  prev = oknots[oik - 1];
  for (kl = oik; prev >= oknots[kl]; kl++) ;

  curr = oknots[kl];
  for (mult = oik; curr < maximum; mult++)
    {
      if (curr < prev)
	goto err112;

      if (prev > curr || curr > prev + eps)
	{

	  /* New knot value found. Fill in old value. */

	   /* numb = (*nik) - oik + mult; */
	  numb = (*nik) - oik + mult + der;
	  if (numb > (*nik))
	    numb = (*nik);


	  /* If numb >= nik, test if all the numb knots are equal
	     or if they only are equal within the resolution eps.
	     If not totally equal knumb=nik-1. */

	  if (numb == (*nik))
	    {
	       /* start = max (kl - oik, 1);
	      stop = kl - 2;
	      for (count1 = start; count1 <= stop; count1++)
		if (oknots[count1 - 1] != oknots[count1])
		  numb = (*nik) - 1; */

	      start = kl - oik + der;
	      stop = kl - 2;
	      for (count1 = start; count1 <= stop; count1++)
		if (oknots[count1] != oknots[count1 + 1])
		  numb = (*nik) - 1;
	    }

	  if (prev == oknots[oik - 1])
	    numb = (*nik);
	  for (count1 = 1; count1 <= numb; count1++)
	    (*nknots)[(*nin)++] = prev;


	  /* Initialize multiplicity. */

	  mult = 0;
	  prev = curr;
	}
      kl++;
      curr = oknots[kl];
    }

  /* Knot for the next last knot value not produced. */

  /* numb = min ((*nik) - oik + mult, (*nik)); */
  numb = min ((*nik) - oik + mult + der, (*nik));


  /* If numb >= nik, test if all the numb knots are equal or if they
   * only are equal within the resolution eps. */

  /* I not totally equal numb=nik-1. */

  if (numb >= (*nik))
    {
       /* start = max (kl - oik, 1);
      stop = kl - 2;
      for (count1 = start; count1 <= stop; count1++)
	if (oknots[count1 - 1] != oknots[count1])
	  numb = (*nik) - 1; */

      start = kl - oik + der;
      stop = kl - 2;
      for (count1 = start; count1 <= stop; count1++)
	if (oknots[count1] != oknots[count1 + 1])
	  numb = (*nik) - 1;
    }

  for (count1 = 1; count1 <= numb; count1++)
    (*nknots)[(*nin)++] = prev;


  /* Knot at et[oin+1] not produced. */

  for (count1 = 1; count1 <= (*nik); count1++)
    (*nknots)[(*nin)++] = maximum;


  /* Knots produced. Correct nin and length of nknots. */

  (*nin) -= (*nik);
  *nknots = increasearray (*nknots, (*nik) + (*nin), DOUBLE);
  if (*nknots == SISL_NULL)
    goto err101;

  goto out;

  /* Not enough memory. */

err101:
  *jstat = -101;
  s6err ("s1894", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1894", *jstat, kpos);
  goto out;

out:
  return;
}
