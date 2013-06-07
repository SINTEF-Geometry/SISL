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
 * $Id: s1917.c,v 1.3 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1917
#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1917 (int inbcrv, double ecoef[], int in2,
       int idim, int eptyp[], double astpar, int iopen,
       double *par[], int *der[], int *inumb, int *jstat)
#else
void
s1917 (inbcrv, ecoef, in2, idim, eptyp, astpar, iopen,
       par, der, inumb, jstat)
     int inbcrv;
     double ecoef[];
     int in2;
     int idim;
     int eptyp[];
     double astpar;
     int iopen;
     double *par[];
     int *der[];
     int *inumb;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To calculate the parametrization for a spline lofted curve
*		interpolation. At the same time the illegal interpolation
*		conditions are discarded and the legal ones along with the end
*		conditions are returned in the input coefficient array.
*
* INPUT      :	inbcrv	- The number of input curves given to the interpolation
*			  problem.
*		ecoef	- Array containing the B-spline verices of the input
*			  curves.
*		in2	- Number of vertices in the B-spline curves.
*		idim	- The dimension of the space in which the curve lie.
*		eptyp	- The decription of each curve.
*			  1 - Ordinary curve.
*			  2 - Knuckle curve.
*			  3 - Tangent to next curve.
*			  4 - Tangent to prior curve.
*			  5 - Double derivative to prior curve.
*			  6 - Double derivative to next curve.
*		astpar	- Start value of parametrization.
*		iopen	- Open/closed surface in lofting direction.
*
* OUTPUT     :	ecoef	- The checked B-spline vertices.
*		eptyp	- The checked type indicators.
*		par	- The parametrization array.
*		der	- The derivative indicarors.
*		inumb	- Number of updated interpolation conditions.
*               jstat    - Status variable:
*                                               > 0     : warning
*                                               = 0     : ok
*                                               < 0     : error
* METHOD     :	The start end-condition, the interpolation conditions and
*		the end end-condition are checked and the legal interpolation
*		conditions produced. If the curve is closed an extra parameter
*		value is produced. At last the sequence of the end
*		end-condition is checked to ensure that the equation system can
*		be solved without pivotation.
*
* REFERENCES :	Fortran version by Tor Dokken, SI, 81-12
*		Line correction by Steinar Meen, KV
*
* CALLS      : s6err.
*
* WRITTEN BY :  Trond Vidar Stensby, SI, 1991-07
* Revised by : Paal Fugelli, SINTEF, Oslo 02/08-1994.  Changed order of
*              occurence of test conditions to avoid overrunning array bounds.
*
*********************************************************************
*/
{
  int kpos = 0;
  int knext;			/* True if next point is being processed. */
  int ktang;			/* True if tangent has been given. */
  int kcurv;			/* True if double derivative has been given. */
  int kerr;			/* True if an error has been found.
				   The error does not cause a program stop, but a
				   default action is taken. */
  int kprev;			/* True if working with previus point. */
  int knumb;			/* Number of exepted interpolation conditions. */
  int kpopar;
  int kkkksm;
  int kant;			/* Number of sucessfully processed interpolation
				   conditions. */
  int ktyp = 0;
  int ki, kj, kl, kr;		/* Loop control variables. */
  int legal;			/* Dummy variables. */
  int dummy1;
  int dummy2;

  double tsum;			/* Variables used when calculating the */
  double tdum;			/* parametrization. */
  int    kdum;
  double tdiff;

  *jstat = 0;


  /* Test if legal input. */

  if (idim < 1)
    goto err102;
  if (inbcrv < 2)
    goto err179;


  /* Allocate space for output parameters. */

  *par = newarray (inbcrv + 1, DOUBLE);
  if (*par == SISL_NULL)
    goto err101;

  *der = newarray (inbcrv, INT);
  if (*der == SISL_NULL)
    goto err101;


  /* Calculate parametrization. */

  kerr = 0;
  knumb = 0;
  kpopar = 0;
  kkkksm = 0;


  /* Initate logical variables taking care of the state of the
     curves, tangents and double derivatives. */

  knext = TRUE;
  ktang = FALSE;
  kcurv = FALSE;


  /* Test curve type for input curves. If internal curves are
     knuckle curves or curves are assigned two derivatives the
     illegal information is discarded. */

  kant = inbcrv;
  ki = 1;

  while (ki <= kant)
    {
      legal = TRUE;

      ktyp = eptyp[ki - 1];

      if (ktyp == 1 || ktyp == 2)
	{
	  if (ktyp == 2 && ki != 1 && ki != inbcrv)
	    kerr = 1;

	  if (knext == FALSE)
	    {
	      ktang = FALSE;
	      kcurv = FALSE;
	    }
	  knext = FALSE;


	  /* Calculate paramerization. The first curve shall
	     have parameter value astpar. */

	  if (kkkksm != 0)
	    {
	      tsum = (double) 0.0;
	      for (kl = 0; kl < in2; kl++)
		{
		  tdum = (double) 0.0;
		  dummy1 = (ki - 1) * in2 * idim + kl * idim;
		  dummy2 = (kkkksm - 1) * in2 * idim + kl * idim;
		  for (kj = 0; kj < idim; kj++)
		    {
		      tdiff = ecoef[dummy1 + kj] - ecoef[dummy2 + kj];
		      tdum += tdiff * tdiff;
		    }
		  tsum += sqrt (tdum);
		}
	      if (tsum == (double) 0.0)
		legal = FALSE;
	      else
		astpar += tsum / (double)in2;
	    }

	  if (legal)
	    {
	      (*der)[knumb] = 0;
	      for (kj = kpopar; kj <= knumb; kj++)
		(*par)[kj] = astpar;

	      kprev = ki;
	      kkkksm = ki;
	      kpopar = knumb + 1;
	    }
	}
      else if (ktyp == 3 || ktyp == 13)
	{
	  /* TANGENT TO NEXT CURVE. */

	  /* Test that tangent not already given. */

	  if (knext && ktang)
	    legal = FALSE;
	  else
	    {
	      if (knext == FALSE)
		kcurv = FALSE;
	      ktang = TRUE;
	      knext = TRUE;

	      /* Legal tanent. */

	      (*der)[knumb] = 1;
	    }
	}
      else if (ktyp == 4 || ktyp == 14)
	{
	  /* TANGENT TO PRIOR CURVE */

	  /* Test that tangent not already given or that we are
	     already working with the next curve. */

	  if (knext || (ktang && !knext))
	    legal = FALSE;
	  else
	    {
	      ktang = TRUE;

	      /* Legal tangent. */

	      (*der)[knumb] = 1;
	      (*par)[knumb] = astpar;
	      kpopar = knumb + 1;
	    }
	}
      else if (ktyp == 5)
	{
	  /* SECOND DERIVATIVE TO NEXT CURVE. */

	  /* Test that double derivative not already given. */

	  if (knext && kcurv)
	    legal = FALSE;
	  else
	    {
	      if (knext == FALSE)
		ktang = FALSE;
	      kcurv = TRUE;
	      knext = TRUE;

	      /* Legal double derivative. */

	      (*der)[knumb] = 2;
	    }
	}
      else if (ktyp == 6)
	{
	  /* SECOND DERIVATIVE TO PRIOR CURVE. */

	  /* Test that double derivative not already given or that we are
	     already working with the next curve. */

	  if (knext || (kcurv && !knext))
	    kcurv = TRUE;

	  /* Legal double derivative. */

	  (*der)[knumb] = 2;
	  (*par)[knumb] = astpar;
	  kpopar = knumb + 1;
	}
      else
	legal = FALSE;


      if (legal)
	{
	  /* Legal curve or derivative. */

	  /* Since curve is legal, the curve is already in the right position. */

	  ki++;
	  knumb++;
	}
      else
	{
	  kerr = 1;

	  /* Decrease the number of conditions by one and copy remaining
	     pointers. */

	  kant--;

	  if (ki <= kant)
	    {
	      memcopy (&eptyp[ki - 1], &eptyp[ki], kant - ki, INT);
	      memcopy (&ecoef[(ki - 1) * in2 * idim], &ecoef[ki * in2 * idim],
		       (kant - ki) * in2 * idim, DOUBLE);
	    }
	}
    }

  /* If tangent or curvature to next specified after end curve, remove this. */

  if (knext)
    {
      /* Remove tangent and/or curvature. */

      if (ktang)
	knumb--;
      if (kcurv)
	knumb--;

      ktang = FALSE;
      kcurv = FALSE;
      kerr = 1;
    }

  if (kpopar > knumb)
    {
      /* Error at end of curve. */

      kerr = 1;
      knumb = kpopar;
    }

  /* If closed surface in lofting direction, find parametrization from
     last to first curve. */

  if (iopen != SISL_CRV_OPEN)
    {
      for (ki = 1; eptyp[ki] != 1 && eptyp[ki] != 2; ki++)
	if (ki > inbcrv)
	  goto err179;


      /* Calculate distance */

      tsum = (double) 0.0;

      for (kl = 0; kl < in2; kl++)
	{
	  tdum = (double) 0.0;
	  dummy1 = (ki - 1) * idim * in2 + kl * idim;
	  dummy2 = (kprev - 1) * idim * in2 + kl * idim;
	  for (kj = 0; kj < idim; kj++)
	    {
	      tdiff = ecoef[dummy1 + kj] - ecoef[dummy2 + kj];
	      tdum += tdiff * tdiff;
	    }
	  tsum += sqrt (tdum);
	}

      /* If tsum is zero the curve must be discarded. */

      if (tsum == (double) 0.0)
	{
	  knumb--;
	  kerr = 1;
	}
      else
	{
	  astpar += tsum / (double)in2;
	  (*par)[knumb] = astpar;
	}
    }

  /* Parametrization completed. */

  /* To get an interpolation problem solvable without pivotation,
     last interpolation condition must be interpolation of position. */

  for (kl = knumb - 2; kl >= 1 && (*par)[kl] >= (*par)[knumb - 1] &&
       (*der)[kl] != 0; kl--) ;

  if (kl >= 1 && (*par)[kl] >= (*par)[knumb - 1])
    {
      /* Interpolation of position at end not last condition.
	 Interchange conditions. */

      (*der)[kl] = (*der)[knumb - 1];
      (*der)[knumb - 1] = 0;
      kdum = eptyp[knumb - 1];
      eptyp[knumb - 1] = eptyp[kl];
      eptyp[kl] = kdum - 1;


      for (kr = kl + 1; kr < knumb - 1; kr++)
	{
	  ktyp = eptyp[kr];
	  if (ktyp == 4 || ktyp == 6)
	    (eptyp[kr])--;
	}

      for (kr = 0; kr < in2; kr++)
	{
	  dummy1 = (knumb - 1) * in2 * idim + kr * idim;
	  dummy2 = kl * in2 * idim + kr * idim;
	  for (kj = 0; kj < idim; kj++)
	    {
	      tdum = ecoef[dummy1 + kj];
	      ecoef[dummy1 + kj] = ecoef[dummy2 + kj];
	      ecoef[dummy2 + kj] = tdum;

	      ktyp = eptyp[kl];
	      if (ktyp >= 10)
		ecoef[dummy2 + kj] = (double) 2.0 *ecoef[dummy1 + kj]
		- ecoef[dummy2 + kj];
	    }
	}
    }

  /* To get an interpolation problem solvable without pivotation, the
     next last interpolation condition should be interpolation of
     tangent if such interpolation is specified. */

  for (kl = knumb - 3; kl >= 1 && (*par)[kl] >= (*par)[knumb - 2] &&
       (*der)[kl] != 1; kl--) ;

  if (kl >= 1 && (*par)[kl] >= (*par)[knumb - 2])
    {
      /* More than one interpolation condition at the end, and interpolation
	 of tangent at end not last but one condition. */

      (*der)[kl] = (*der)[knumb - 2];
      (*der)[knumb - 2] = 1;
      kdum = eptyp[knumb - 2];
      eptyp[knumb - 2] = eptyp[kl];
      eptyp[kl] = kdum;

      for (kr = 0; kr < in2; kr++)
	{
	  dummy1 = (knumb - 2) * in2 * idim + kr * idim;
	  dummy2 = kl * in2 * idim + kr * idim;
	  for (kj = 0; kj < idim; kj++)
	    {
	      tdum = ecoef[dummy1 + kj];
	      ecoef[dummy1 + kj] = ecoef[dummy2 + kj];
	      ecoef[dummy2 + kj] = tdum;
	    }
	}
    }

  if (knumb < 2)
    goto err187;

  /* OK */

  *inumb = knumb;

  if (kerr == 1)
    *jstat = 102;
  else
    goto out;


  /*  Not able to allocate memory. */

err101:
  *jstat = -101;
  s6err ("s1917", *jstat, kpos);
  goto out;

  /* Error in description of spline. */

err102:
  *jstat = -102;
  s6err ("s1917", *jstat, kpos);
  goto out;

  /* Error in description if interpolation problem. */

err179:
  *jstat = -179;
  s6err ("s1917", *jstat, kpos);
  goto out;

  /* Error in description of input curves. */

err187:
  *jstat = -187;
  s6err ("s1917", *jstat, kpos);
  goto out;

out:
  return;
}
