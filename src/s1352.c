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
 * $Id: s1352.c,v 1.5 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1352

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1352(double t[],int n,int k,double inteps[],double lefteps[],
	   double righteps[],int dim,int leftfix,int rightfix,
	   double eps[],int *stat)
#else
void s1352(t, n, k, inteps, lefteps, righteps, dim, leftfix, rightfix,
	   eps, stat)
     double	t[];
     int	n;
     int	k;
     double 	inteps[];
     double	lefteps[];
     double	righteps[];
     int	dim;
     int	leftfix;
     int	rightfix;
     double	eps[];
     int	*stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* Purpose: To compute a variable tolerance, eps, along a knot vector. This
*          tolerance will be used for knot removal of tensor product B-spline
*          surfaces. The tolerance is to be lefteps at the left end of the
*          knot vector, then increase to inteps at the midpoint and then
*          decrease again to righteps at the right end of the knot vector
*          (this description applies to each of the dim components of the
*          tolerance.) The way the tolerance increases from the left end to
*          the midpoint, is determined by the variable leftfix.
*          The tolerance will increase according to the first half
*          of a perfect B-spline of order leftfix+1 (if leftfix=0 the
*          tolerance will be inteps along the first half of the knot vector
*          irrespective of the value of lefteps) translated and scaled
*          appropriately. The value of eps(.,i)
*          will be the value of this spline at the 'de Boor point',
*
*                           (t[i+1]+...+t[i+k-1])/(k-1)),
*
*          assuming that this point lies to the left of the midpoint.
*          To the right of the midpoint the tolerance will decrease in a
*	   similar way depending on the integer irend.
*
*
*
*
* Input:
*          t           - the knot vector along which we are to compute
*			 tolerances.
*
*          n           - the number of B-splines associated by t.
*
*          k           - the polynomial order of the B-splines associated
*                        with t.
*
*          inteps      - double array of length dim giving the maximum
*			 tolerance in the interior of t.
*
*	   lefteps     - double array of length dim giving the tolerance
*                        at the beginning of t.
*
*          righteps    - double array of length dim giving the tolerance
*                        at the end of t.
*          dim         - the dimension of the tolerance vectors.
*
*          leftfix     - the number of derivatives to be kept fixed during
*                        knot removal at the edge corresponding to the
*                        beginning of t. The greater the value of leftfix
*                        the more slowly the tolerance is allowed to increase
*			 towards the middle of t.
*          rightfix    - similar to leftfix.
*
*
* Output:
*          eps         - double array of length dim containing the computed
*                        tolerance vector.
*
* Method:
*     As indicated above the main ingredient in this routine is the perfect
*     B-spline. First the knots for the perfect B-spline of order leftfix+1
*     are computed and the tolerance over the first half of the knot vector
*     is calculated. Then the knots for the perfect B-spline of order
*     rightfix+1 are computed and the tolerance over the second half of
*     the knot vector is calculated.
*
*
* References:
*     Larry Schumaker, Spline Functions: Basic Theory, Wiley.
*
*
*-
*
* Calls: s1221
*
* WRITTEN BY : Knut Moerken, University of Oslo, July 1992, based
*              on an earlier Fortran version.
* Changed by: Paal Fugelli, SINTEF, 1994-07.
*             Added code at end to to fix memory leakage problem.
*
*********************************************************************
*/
{
  int ih = 3*MAX(leftfix, rightfix) + 2; /* Size of knot vector.       */
  double *th = SISL_NULL;			 /* Knot vector.               */
  double *hcoef = SISL_NULL;			 /* B-spline coefficients.     */
  double *w = SISL_NULL;			 /* Work array for s1221.      */

  SISLCurve *hspline = SISL_NULL;             /* Object that will be used
					    for storing the perfect
					    B-spline.		       */

                                         /* The usage of the rest of
					    the variables should be
					    clear from the code.       */

  int left=0, i, lstat, j, jh, pos=0;
  double k1inv, ta, tb, tmid, hh, maxval, sum, cnst, ch, val;

  /* A useful constant. */

  if (k == 1) k1inv = 1.0;
  else k1inv = 1.0 / (float) (k-1);

  /* Determine the left and right ends and the midpoint of the knot
     interval over which the tolerance vector is to be computed,
     hh will also be useful. */

  ta = t[k-1]; tb = t[n];
  tmid = ta + (tb-ta)/2.0;
  hh = 1.0 / (tb-ta);

  /* Allocate space for the knot vector of the perfect B-spline. For
     simplicity we use s1221 for computing values on this B-spline and
     then it is easiest if it is a B-spline on a knot vector with
     multiple knots at the ends. If the order of the perfect B-spline is
     r then we use array of length 3r-1, so that the perfect B-spline
     becomes the middle one on this knot vector. Since we do not
     bother to allocate space twice we make sure we have enough space
     for largest of these knot vectors. */

  th = newarray(ih, double);
  if (th == SISL_NULL) goto err101;

  /* The length of the coefficient vector should be clear from
     the argument above. We set all coefficients to zero except the one
     that multiplies the perfect B-spline which of course must be one. */

  hcoef = new0array(2*MAX(leftfix, rightfix)+1, double);
  if (hcoef == SISL_NULL) goto err101;

  hcoef[leftfix] = 1.0;

  /* Determine the knot vector. */

  if (leftfix == 0)
    {
      /* For order one we only need two knots which we could really pick
	 quite arbitrarily as long as we avoid the discontinuity. */

      th[0] = MIN(-1.1, (2*t[0]-ta-tb)*hh);
      th[1] = MAX(1.1, (2*t[n+k-1]-ta-tb)*hh);
    }
  else
    {

      /* Compute the knots of the perfect B-spline and give the first and
	 last knot maximal multiplicity. */

      cnst = PI / (float) (leftfix+1);

      /* First the "real knots". */

      for (j=0; j<=leftfix+1; j++)
	th[leftfix+j] = cos((leftfix+1-j)*cnst);

      /* Increase the multiplicity at the ends (-1 and 1 if the arithmetic
	 is exact). */

      for (j=0; j<leftfix; j++)
	{
	  th[j] = th[leftfix];
	  th[j+2*leftfix+2] = th[j+2*leftfix+1];
	}
    }

  /* Allocate space for a work array to be used by s1221 which will also
     be used for returning results. */

  w = newarray(MAX((leftfix+1),(rightfix+1)), double);
  if (w == SISL_NULL) goto err101;

  /* Create the spline object that represents the perfect B-spline. */

  hspline = newCurve(2*leftfix+1, leftfix+1, th, hcoef, 1, 1, 0);
  if (hspline == SISL_NULL) goto err101;

  /* Evaluate the perfect B-spline at its maximum, the midpoint 0,
     and save it in maxval. */

  s1221(hspline, 0, 0.0, &left, w, &lstat);
  if (lstat < 0) goto err;

  maxval = w[0];

  /* Now we are ready to compute tolerances. To do this we stretch the
     perfect B-spline so that it starts at ta and ends at tb and therefore
     has its maximum at tmid. To compute tolerances in the left half
     of [ta, tmid] we use the first half of the perfect B-spline, adjusted
     so that its value at ta is lefteps and its value at tmid is inteps.
     The right half is treated similarly. The tolerance for coefficient
     i will then from the adjusted perfect B-spline at the de Boor point
     (t[i+1]+...+t[i+k-1])/(k-1). */

  /* We start by computing the tolerance at the left end which is lefteps. */

  if (leftfix > 0)
    {
      for (ih=0; ih<dim; ih++) eps[ih] = lefteps[ih];
    }
  else
    {
      for (ih=0; ih<dim; ih++) eps[ih] = inteps[ih];
    }

  /* We then compute the second de Boor point (sum). */

  if (k == 1) sum = t[1];
  else
    {
      sum = 0.0;

      for (j=2; j<k+1; j++) sum += t[j];

      sum *= k1inv;
    }

  /* Then we enter a loop to compute the tolerances up to the midpoint. */

  i = 1;  ih = dim;

  while (sum < tmid)
    {

      /* Determine the correct argument to the perfect B-spline. */

      ch = (2*sum-ta-tb)*hh;

      /* Evaluate the perfect B-spline. */

      s1221(hspline, 0, ch, &left, w, &lstat);
      if (lstat < 0) goto err;

      /* Scale the value. */

      val = w[0];
      val /= maxval;

      /* Compute the tolerance at this de Boor point. */

      for (jh=0; jh<dim; jh++, ih++)
	eps[ih] = MIN((inteps[jh]-lefteps[jh])*val + lefteps[jh],inteps[jh]);

      ++i;

      /* Compute the next de Boor point. */

     if (k == 1) sum = t[i];
     else
       {
	 sum = 0.0;

	 for (j=i+1; j<i+k; j++) sum += t[j];

	 sum *= k1inv;
       }
    }

  /* Now we have to prepare for computing tolerances in [tmid,ta].
     If the number of derivatives to be kept fixed is deifferent at the
     two ends we must compute a new knot vector. */

  if (leftfix != rightfix)
    {
      if (rightfix == 0)
	{
	  th[0] = MIN(-1.1, (2*t[0]-ta-tb)*hh);
	  th[1] = MAX(1.1, (2*t[n+k-1]-ta-tb)*hh);
	}
      else
	{
	  cnst = PI / (float) (rightfix+1);

	  for (j=0; j<=rightfix+1; j++)
	    th[rightfix+j] = cos((rightfix+1-j)*cnst);
	  for (j=0; j<rightfix; j++)
	    {
	      th[j] = th[rightfix];
	      th[j+2*rightfix+2] = th[j+2*rightfix+1];
	    }
	}

      /* For simplicity we always update the coefficient vector. */

      hcoef[leftfix] = 0.0;
      hcoef[rightfix] = 1.0;

      /* For simplicity we also free the old spline object (has icopy==0). */

      freeCurve(hspline);
      hspline = newCurve(2*rightfix+1, rightfix+1, th, hcoef, 1, 1, 0);
      if (hspline == SISL_NULL) goto err101;

      /* Find the maximum of the perfect B-spline. */

      s1221(hspline, 0, 0.0, &left, w, &lstat);
      if (lstat < 0) goto err;

      maxval = w[0];
    }

  /* Compute the remaining tolerances as above (the next de Boor point is
     already known). */

  while (i<n-1)
    {
      ch = (2*sum-ta-tb)*hh;

      s1221(hspline, 0, ch, &left, w, &lstat);
      if (lstat < 0) goto err;

      val = w[0];
      val /= maxval;

      for (jh=0; jh<dim; jh++, ih++)
	eps[ih] = MIN((inteps[jh]-righteps[jh])*val + righteps[jh],inteps[jh]);

      ++i;

     if (k == 1) sum = t[i];
     else
       {
	 sum = 0.0;

	 for (j=i+1; j<i+k; j++) sum += t[j];

	 sum *= k1inv;
       }
    }

  /* Compute the final tolerance which is righteps. */

  if (rightfix > 0)
    {
      for (jh=0; jh<dim; jh++, ih++) eps[ih] = righteps[jh];
    }
  else
    {
      for (jh=0; jh<dim; jh++, ih++) eps[ih] = inteps[jh];
    }

  /* Everything went well. */

  *stat = 0;
  goto out;

  /* Error in memory allocation. */

 err101:
  *stat = -101;
  s6err("s1352", *stat, pos);

  goto out;

  /* Error in lower level routine. */

 err:
  *stat = lstat;
  s6err("s1352", *stat, pos);

  /* Free memory. */

 out:
  if (w != SISL_NULL) freearray(w);
  if (hspline != SISL_NULL) freeCurve(hspline);

  /* Must also remove knots and coefs because hspline had icopy==0. */

  if (th != SISL_NULL) freearray(th);
  if (hcoef != SISL_NULL) freearray(hcoef);

  return;

}
