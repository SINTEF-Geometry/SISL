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
 * $Id: s1916.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1916

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1916 (int inbcrv, double et2[], double ecoef[], int in2, int iord, int idim,
       int iopen, double par[], int der[], int *jstat)
#else
void
s1916 (inbcrv, et2, ecoef, in2, iord, idim, iopen, par, der, jstat)
     int inbcrv;
     double et2[];
     double ecoef[];
     int in2;
     int iord;
     int idim;
     int iopen;
     double par[];
     int der[];
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To adjust the length of the derivative curves such that the
*		resulting product of the derivative of the lofting surface
*		corresponds better to the curve length.
*
* INPUT      :	inbcrv	- The number of input curves given to the interpolation
*			  problem.
*		et2	- The knot vector of the B-spline curves.
*		ecoef	- Array containing the B-spline vertices of the input
*			  curves.
*		in2	- Number of vertices in the B-spline curves.
*		iord1	- The order of the B-spline lofting.
*		idim	- The dimension of the space in which the curves lie.
*		eptyp	- The description of each curve.
*			  1 - Ordinary curve.
*			  2 - Knuckle curve, is treated as ordinary curve.
*			  3 - Tangent to next curve.
*			  4 - Tangent to prior curve.
*			  5 - Double derivative to prior curve.
*			  6 - Double derivative to next curve.
*		iopen	- Open/closed surface in lofting direction
*		astpar	- Start parameter value.
*
* OUTPUT     :	ecoef	- The adjusted B-spline vertices.
*		etyp	- The type indicators.
*		par	- The parametrization.
*		der	- The derivative indicators.
*		jstat    - Status variable:
*                                               > 0     : warning
*                                               = 0     : ok
*                                               < 0     : error
*
* METHOD     :
*
* REFERENCES :	Fortran version by Tor Dokken, SI, 1984-10
*
* CALLS      :	s1919,s6err.
*
* WRITTEN BY :	Trond Vidar Stensby, SI, 1991-07
*
*********************************************************************
*/
{
  int kpekc;			/* Pointer to current curve. */
  int kpekp;			/* Pointer to previous curve. */
  int kpekf;			/* Pointer to following curve. */
  int kstart;			/* Used for controllong loops. */
  int kstop;
  int kip;			/* Indicates if previous curve exists. */
  int kif;			/* Indicates if following curve exists. */
  int ki, kj;			/* Loop control variable. */

  double tc;			/* Parameter value for current curve. */
  double tp;			/* Parameter value for previous curve. */
  double tf;			/* Parameter value for following curve. */
  double *prev = SISL_NULL;		/* Vertices of previous curve. */
  double *curr = SISL_NULL;		/* Vertices of current curve. */
  double *deriv = SISL_NULL;		/* Vertices of derivative curve. */
  double *follow = SISL_NULL;	/* Vertices of following curve. */
  int kstat = 0;
  int kpos = 0;

  *jstat = 0;


  /* Test if legal input. */

  if (in2 < iord || iord < 1)
    goto err112;

  /* Allocate space for the curves. */

  prev = newarray (in2 * idim, DOUBLE);
  if (prev == SISL_NULL)
    goto err101;
  curr = newarray (in2 * idim, DOUBLE);
  if (curr == SISL_NULL)
    goto err101;
  deriv = newarray (in2 * idim, DOUBLE);
  if (deriv == SISL_NULL)
    goto err101;
  follow = newarray (in2 * idim, DOUBLE);
  if (follow == SISL_NULL)
    goto err101;

  for (ki = 0; ki < inbcrv; ki++)
    {
      if (der[ki] == 1)
	{
	  /* Derivative curve found.
	     kpekp is to be the pointer to the previous curve.
	     kpekc is to be the pointer to the current curve.
	     kpekf is to be the pointer to the following curve. */

	  kpekc = ki;
	  tc = par[ki];

	  /* Find previous curve if any. */

	  kstop = ki;
	  kip = 0;
	  for (kj = 1; kj <= kstop && kip == 0; kj++)
	    {
	      kpekp = ki - kj;
	      tp = par[kpekp];
	      if (der[kpekp] == 0 && tp == tc)
		kpekc = kpekp;
	      if (der[kpekp] == 0 && tp < tc)
		kip = 1;
	    }

	  /* Find following curve if any. */

	  kstart = ki + 1;
	  kif = 0;
	  for (kj = kstart; kj < inbcrv && kif == 0; kj++)
	    {
	      kpekf = kj;
	      tf = par[kpekf];
	      if (der[kpekf] == 0 && tf == tc)
		kpekc = kpekf;
	      if (der[kpekf] == 0 && tf > tc)
		kif = 1;
	    }

	  /* Previous, current and next position curve found if any.
	     Check that at least current and one of previous or next is found. */

	  if (kpekc == ki || (kip == 0 && kif == 0))
	    goto err186;


	  /* Copy the current curves into temporary arrays. */

	  if (kip == 1)
	    memcopy (prev, &ecoef[kpekp * in2 * idim], in2 * idim, DOUBLE);

	  memcopy (curr, &ecoef[kpekc * in2 * idim], in2 * idim, DOUBLE);
	  memcopy (deriv, &ecoef[ki * in2 * idim], in2 * idim, DOUBLE);

	  if (kif == 1)
	    memcopy (follow, &ecoef[kpekf * in2 * idim],
		     in2 * idim, DOUBLE);

	  s1919 (et2, prev, curr, deriv, follow, in2, iord, idim, kip,
		 kif, tp, tc, tf, &kstat);
	  if (kstat < 0)
	    goto error;


	  /* Copy the adjusted derivative curve back. */

	  memcopy (&ecoef[ki * in2 * idim], deriv, in2 * idim, DOUBLE);
	}
    }

  /* OK */

  goto out;

  /* Error in allocation. */

err101:
  *jstat = -101;
  s6err ("s1916", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1916", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1916", *jstat, kpos);
  goto out;

  /* Special error. */

err186:
  *jstat = -186;
  s6err ("s1916", *jstat, kpos);
  goto out;

out:
  if (prev != SISL_NULL)
    freearray (prev);
  if (curr != SISL_NULL)
    freearray (curr);
  if (deriv != SISL_NULL)
    freearray (deriv);
  if (follow != SISL_NULL)
    freearray (follow);

  return;
}
