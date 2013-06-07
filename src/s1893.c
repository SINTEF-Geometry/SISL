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
 * $Id: s1893.c,v 1.3 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1893

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1893 (SISLCurve * orig, double earray[], int dimp1, int narr, int der1,
       int der2, SISLCurve ** ncurve, int *jstat)

#else
void
s1893 (orig, earray, dimp1, narr, der1, der2, ncurve, jstat)
     SISLCurve *orig;
     double earray[];
     int dimp1;
     int narr;
     int der1;
     int der2;
     SISLCurve **ncurve;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :  To make the function
*		F(t) = (D**der1)P(t) earray (D**der2)P(t),
*		If earray is the identity matrix then we make the dot
*		product of ider1-derivate and ider2-derivate of the B-spline
*		curve P(t) described by orig.
*
* INPUT      :  orig	- The original B-spline curve, P(t).
*		earray	- The array to be used.
*		dimp1	- The dimension orig.idim+1. This is due to that the
*			  matrix earray is used for calculation with homogenous
*			  coordinates.
*		narr	- Number of parallel arrays.
*		der1	- The order of the first derivate to be used.
*		der2	- The order of the second derivate to be used.
*
* OUTPUT     :  ncurve	- The new B-spline curve.
*		jstat	- Status variable.
*                        > 0     : warning
*                        = 0     : ok
*                        < 0     : error
*
* METHOD     :  The product of the der1 derivate of a B-spline curve a matrix
*		and the der2 derivate of the same B-spline curve is a B-spline
*		function of order 2*ik-der1-der2, if the order of the B-spline
*		curve is ik.
*		The subroutine has the following steps.
*		1. Calculate knot vector of B-spline function.
*		2. Calculate parameter value to be used for calculation of
*		   points on the B-spline function.
*		3. Interpolate the calculated values on the B-spline function.
*		Note: the parameter values of the B-spline function must be
*		calculated to ensure reproduction of the B-spline function.
*
* REFERENCES :  Fortran version
*		Tor Dokken, SI, 1982-02
*
* CALLS      :  s1894, s1890, s1221, s1891, s6err.
*
* WRITTEN BY : Trond Vidar Stensby, SI, 1991-06
* REVISED BY : Johannes Kaasa, SI, May 1992 (Corrected the sequence of
*              indexes in the product earray*val1*val2)
* REVISED BY : Michael Floater, SI, June 92. The rational stuff
*              was completely messed up. But it works now.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, September 1994.  Size of
*              'tau' must be according to parameters passed to s1891().
*********************************************************************
*/
{
  int nik;			/* The order of the new basis. */
  int nin;			/* The number of verices in the new basis. */
  int mder;			/* max(der1,der2) */
  int left;			/* Interval indicator. */
  int pos;			/* Used to index earray */
  int pos1;			/* Position of the first derivatives in the
				 * array deriv. (returned form s1221); */
  int pos2;			/* Position of the second derivatives in the
				 * array deriv. (returned form s1221); */
  int count1, count2;		/* Loop control variables. */
  int count3=0;
  int kr, kl, kp;
  int *der = SISL_NULL;		/* The derivative indicators. (0) */

  double *nknots = SISL_NULL;	/* The new knot vector. */
  double *coef = SISL_NULL;		/* Coefficients of the new B-spline curve. */
  double *par = SISL_NULL;		/* Parameter values used for interpolation. */
  double *deriv = SISL_NULL;		/* The derivates returned by s1221. */
  double *val1 = SISL_NULL;		/* Extracted values from deriv. */
  double *val2 = SISL_NULL;		/* Extracted values from deriv. */
  double *tau = SISL_NULL;		/* Interpolation points. */
  double sum;			/* Used for calculating F(t). */
  int kpos = 0;
  int kstat = 0;

  *jstat = 0;


  /* Test if legal input. */

  if (orig->ik <= 1 || orig->in <orig->ik)
    goto err112;
  if ( dimp1 < orig->idim || dimp1 > orig->idim +1 )
    goto err151;


  /* Produce a knot vector. */

  s1894 (orig->et, orig->ik, orig->in, der1, der2, earray, dimp1, narr,
	 &nknots, &nik, &nin, &kstat);
  if (kstat < 0)
    goto error;


  /* Produce parameter values and derivate indicators. */

  s1890 (nknots, nik, nin, &par, &der, &kstat);
  if (kstat < 0)
    goto error;


  /* Allocate arrays. */

  val1 = newarray (orig->idim + 1, DOUBLE);
  if (val1 == SISL_NULL)
    goto err101;
  val2 = newarray (orig->idim + 1, DOUBLE);
  if (val2 == SISL_NULL)
    goto err101;
  tau = new0array (nin * narr * narr, DOUBLE);
  /*  tau = newarray (nin * narr, DOUBLE);  (PFU 21/09-94) */
  if (tau == SISL_NULL)
    goto err101;

  mder = max (der1, der2);
  deriv = newarray ((mder + 1) * orig->idim, DOUBLE);
  if (deriv == SISL_NULL)
    goto err101;


  /* Calculate interpolation points. */

  left = 0;
  for (count1 = 0; count1 < nin; count1++)
    {
      s1221 (orig, mder, par[count1], &left, deriv, &kstat);
      if (kstat < 0)
	goto error;


      /* Extract the values/derivatives. */

      pos1 = der1 * orig->idim;
      pos2 = der2 * orig->idim;

      for (count2 = 0; count2 < orig->idim; count2++)
	{
	  val1[count2] = deriv[pos1++];
	  val2[count2] = deriv[pos2++];
	}

      if(orig->idim < dimp1)
      {
          if (der1 > 0)
            val1[orig->idim] = (double) 0.0;
          else
            val1[orig->idim] = (double) 1.0;

          if (der2 > 0)
            val2[orig->idim] = (double) 0.0;
          else
            val2[orig->idim] = (double) 1.0;
      }

      /* Calculate the functtion F(t). */

      pos = 0;
      for (kl = 0; kl < narr; kl++)
	{
	  sum = (double) 0.0;
	  for (kr = 0; kr < dimp1; kr++)
	    {
	      for (kp = 0; kp < dimp1; kp++)
		sum += earray[pos++] * val1[kr] * val2[kp];
	      /* sum += earray[pos++] * val1[kp] * val2[kr]; */
	    }
	  tau[count3++] = sum;
	}
    }

  /* Caculate new curve description */

  s1891 (par, tau, narr, nin, narr, der, TRUE, nknots, &coef, &nin,
	 nik, 0, 0, &kstat);
  if (kstat < 0)
    goto error;

  *ncurve = newCurve (nin, nik, nknots, coef, orig->ikind, narr, 2);
  if (*ncurve == SISL_NULL)
    goto err171;
  (*ncurve)->cuopen = orig->cuopen;

  /* OK */

  goto out;


  /* Memory error. */

err101:
  *jstat = -101;
  s6err ("s1893", *jstat, kpos);
  goto out;

  /* Could not create curve. */

err171:
  *jstat = -171;
  s6err ("s1893", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1893", *jstat, kpos);
  goto out;

  /* dimp1 not equal to idim+1. */

err151:
  *jstat = -151;
  s6err ("s1893", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1893", *jstat, kpos);
  goto out;

  /* Free memory. */

out:
  if (val1 != SISL_NULL)
    freearray (val1);
  if (val2 != SISL_NULL)
    freearray (val2);
  if (der != SISL_NULL)
    freearray (der);
  if (par != SISL_NULL)
    freearray (par);
  if (deriv != SISL_NULL)
    freearray (deriv);
  if (tau != SISL_NULL)
    freearray (tau);
  return;
}
