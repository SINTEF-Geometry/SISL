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
 * $Id: s1896.c,v 1.3 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1896

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1896 (SISLSurf * osurf, double earray[], int dimp1, int narr, int ders1[],
       int dert1[], int ders2[], int dert2[], SISLSurf ** nsurf, int *jstat)

#else
void
s1896 (osurf, earray, dimp1, narr, ders1, dert1, ders2, dert2, nsurf, jstat)
     SISLSurf *osurf;
     double earray[];
     int dimp1;
     int narr;
     int ders1[];
     int dert1[];
     int ders2[];
     int dert2[];
     SISLSurf **nsurf;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :  For each of the planes in the matrix we make the function
*
*		     iders1	idert1		     iders2    idert2
*		  /D \       /D \		  /D \      /D \
*	F(s,t) =   --    *    --   P(s,t) earray   --   *    --    P(s,t)
*              	  \DS/       \DT/		  \DS/      \DT/
*
*		Here:
*			iders1=ders1[i], iders2=ders2[i]
*			idert1=dert1[i], idert2=dert2[i]
*
*		If earray is the identity matrix then we make the dor product
*		of (iders1,idert1)-derivateive and (iders2,idert2)-derivative
*		of the B-spline "surface" P(s,t) described by osurf.
*
* INPUT      :  osurf	- The original B-pline surface.
*		earray	- The array to be used.
*		dimp1	- The dimension idim+1. This is due to that the matrix
*			  earray is used for calculation with homogenous
*			  coordinates.
*		narr	- Number of parallel arrays.
*		ders1	- The derivatives to be applied in the first parameter
*			  direction in the first occurence of the surface.
*		dert1	- The derivatives to be applied in the second parameter
*			  direction in the first occurence of the surface.
*		ders2	- The derivatives to be applied in the first parameter
*			  direction in the second occurence of the surface.
*		dert2	- The derivatives to be applied in the second parameter
*			  direction in the second occurence of the surface.
*
* OUTPUT     :  nsurf	- The new B-spline surface.
*		jstat	- Status variable
*                                               > 0     : warning
*                                               = 0     : ok
*                                               < 0     : error
*
* METHOD     :  The product of the (iders1,idert1)-derivative of a B-spline
*		curve, a matrix and the (iders2,idert2)-derivative of the same
*		B-spline curve is a function of order 2*ik1-ders1-ders2-1 in
*		the first parameter direction and 2*ik2-dert1-dert2-1 in the
*		second parameter direction, if the order of the B-spline
*		surface is ik1,ik2.
*		ds1 = min( ders1[i] ),	i=1,...,narr
*		dt1 = min( dert1[i] ),	i=1,...,narr
*		ds2 = min( ders2[i] ),	i=1,...,narr
*		dt2 = min( dert2[i] ),	i=1,...,narr
*
*		The subroutine have the following steps.
*		1. Calculate knot vectors of B-spline functions.
*		2. Calculate parameter value pairs to be used for calculation
*		   of points on the B-spline function.
*		3. Interpolate the calculated values on the B-spline function.
*
*		Note: The parameter values of the B-spline function values must
*		be calculated to ensure reproduction of the B-spline function.
*
* REFERENCES :  Fortran version:
*		Tor Dokken, SI, 1983-05
*
* CALLS      :  s1894,s1890,s1421,s1891,s6err.
*
* WRITTEN BY :  Trond Vidar Stensby, SI, 1991-06
* REVISED BY :  Johannes Kaasa, SI, May 92 (Changed the order of the indexes
*               in the multiplication of earray, val1 and val2. In addition I
*               changed the sequence in the second assignment of ds2 and dt2.
*               I also had to put on paranthesises in the expression for pos1
*               and pos2).
* REVISED BY : Michael Floater, SI, June 92. The rational stuff
*              was completely messed up. But it works now.
*              Dimension of output curve was wrong -- now narr.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, August 94. Fixed over-running
*              of array indices (found with Purify).
*
*********************************************************************
*/
{
  int nik1;			/* Order of new surface in
				   first parameter direction. */
  int nin1;			/* Order of new surface in
				   second parameter direction. */
  int nik2;			/* Number of vertices in first
				   parameter direction. */
  int nin2;			/* Number of vertices in second parameter direction. */
  int lfs;			/* Interval indicator. (left side) */
  int lft;			/* Interval indicator. (left side) */
  int tpos;			/* Used to index array tau. */
  int epos;			/* Used to index earray. */
  int pos1;			/* Position of values of first derivatives. */
  int pos2;			/* Position of values of second derivatives. */
  int ds1;			/* Order of derivatives. */
  int dt1;
  int ds2;
  int dt2;
  int mds1;			/* Maximum order of derivatives. */
  int mdt1;
  int mds2;
  int mdt2;
  int nder1;			/* Total order of derivatives.
				   (Both directions) */
  int nder2;
  int dim;			/* Dimension of tau. */
  int maxder;			/* Largest total order of derivatives.
				   (Both functions.) */
  int count1;			/* Loop control variables. */
  int kj, ki;
  int kl, kr, kp;
  double parval[2];
  double sum;			/* Used for calculation of P(s,t). */
  double *nknots1 = SISL_NULL;	/* New knots in first parameter direction. */
  double *nknots2 = SISL_NULL;	/* New knots in second parameter direction. */
  double *coef1 = SISL_NULL;		/* New coeficients */
  double *coef2 = SISL_NULL;		/* New coeficients */
  double *par1 = SISL_NULL;		/* Parameter values in first direction. */
  double *par2 = SISL_NULL;		/* Parameter values in second direction. */
  int *der1 = SISL_NULL;		/* Derivative indicators in first direction. */
  int *der2 = SISL_NULL;		/* Derivative indicators in second direction.*/
  double *deriv = SISL_NULL;		/* Derivatives returned by s1421. */
  double *normal = SISL_NULL;	/* Normal returned by s1421. (not used) */
  double *val1 = SISL_NULL;		/* Values extracted from deriv. */
  double *val2 = SISL_NULL;		/* Values extracted from deriv. */
  double *tau = SISL_NULL;		/* Interpolation points. */
  int kstat = 0;
  int kpos = 0;

  *jstat = 0;

  /* Test if legal input. */

  if (osurf->ik1 <= 1 || osurf->in1 < osurf->ik1)
    goto err112;
  if (osurf->ik2 <= 1 || osurf->in2 < osurf->ik2)
    goto err112;

  /* Find minimal and maximal order of derivatives */

  ds1 = mds1 = ders1[0];
  dt1 = mdt1 = dert1[0];
  ds2 = mds2 = ders2[0];
  dt2 = mdt2 = dert2[0];

  for (count1 = 1; count1 < narr; count1++)
    {
      if (ds1 > ders1[count1])	ds1 = ders1[count1];
      if (dt1 > dert1[count1])	dt1 = dert1[count1];
      if (ds2 > ders2[count1])	ds2 = ders2[count1];
      if (dt2 > dert2[count1])	dt2 = dert2[count1];

      if (mds1 < ders1[count1])	mds1 = ders1[count1];
      if (mdt1 < dert1[count1])	mdt1 = dert1[count1];
      if (mds2 < ders2[count1])	mds2 = ders2[count1];
      if (mdt2 < dert2[count1]) mdt2 = dert2[count1];
    }

  /* Produce a knot vector in the first parameter direction. */

  s1894 (osurf->et1, osurf->ik1, osurf->in1, ds1, ds2, earray, dimp1, narr,
	 &nknots1, &nik1, &nin1, &kstat);
  if (kstat < 0) goto error;

  /* Produce a knot vector in second parameter direction. */

  s1894 (osurf->et2, osurf->ik2, osurf->in2, dt1, dt2, earray, dimp1, narr,
	 &nknots2, &nik2, &nin2, &kstat);
  if (kstat < 0) goto error;

  /* Produce parameter values and derivative indicators in first
   * parameter direction. */

  s1890 (nknots1, nik1, nin1, &par1, &der1, &kstat);
  if (kstat < 0) goto error;

  /* Produce parameter values and derivative indicators in second
   * parameter direction. */

  s1890 (nknots2, nik2, nin2, &par2, &der2, &kstat);
  if (kstat < 0) goto error;

  /* Allocate memory for point calculation. */

  val1 = newarray (dimp1, DOUBLE);
  if (val1 == SISL_NULL) goto err101;
  val2 = newarray (dimp1, DOUBLE);
  if (val2 == SISL_NULL) goto err101;
  tau = newarray (narr * nin1 * nin2, DOUBLE);
  if (tau == SISL_NULL) goto err101;
  maxder = max (max (mds1, mds2), max (mdt1, mdt2));
  deriv = newarray (osurf->idim * (maxder + 1) * (maxder + 2) / 2, DOUBLE);
  if (deriv == SISL_NULL) goto err101;
  normal = newarray (osurf->idim * (maxder + 1) * (maxder + 2) / 2, DOUBLE);
  if (normal == SISL_NULL) goto err101;

  /* Calculate interpolation points. */

  lfs = 0;
  lft = 0;
  tpos = 0;
  for (kj = 0; kj < nin2; kj++)
    {
      parval[1] = par2[kj];
      for (ki = 0; ki < nin1; ki++)
	{
	  parval[0] = par1[ki];
	  epos = 0;
	  for (kl = 0; kl < narr; kl++)
	    {
	      ds1 = ders1[kl];
	      dt1 = dert1[kl];
	      ds2 = ders2[kl];
	      dt2 = dert2[kl];

	      /* ds2 = dert2[kl];
	      dt2 = ders2[kl]; */

	      maxder = max (max (ds1, ds2), max (dt1, dt2));

	      s1421 (osurf, maxder, parval, &lfs, &lft, deriv, normal, &kstat);
	      if (kstat < 0) goto error;

	      nder1 = ds1 + dt1;
	      nder2 = ds2 + dt2;
	      pos1 = osurf->idim * (nder1 * (nder1 + 1) / 2 + dt1);
	      pos2 = osurf->idim * (nder2 * (nder2 + 1) / 2 + dt2);

	      for (count1 = 0; count1 < osurf->idim; count1++)
		{
		  val1[count1] = deriv[pos1++];
		  val2[count1] = deriv[pos2++];
		}
	      if (osurf->idim < dimp1)
		{
		  val1[osurf->idim] = (double) 1.0;
		  val2[osurf->idim] = (double) 1.0;
		  if (ds1 > 0 || dt1 > 0)
		    val1[osurf->idim] = (double) 0.0;
		  if (ds2 > 0 || dt2 > 0)
		    val2[osurf->idim] = (double) 0.0;
		}

	      /* Can now calculate a interpolation point. */

	      sum = (double) 0.0;
	      for (kr = 0; kr < dimp1; kr++, epos += dimp1)
		{
		  for (kp = 0; kp < dimp1; kp++)
		    sum += earray[epos + kp] * val1[kr] * val2[kp];
		  /* sum += earray[epos + kp] * val1[kp] * val2[kr]; */
		}
	      tau[tpos++] = sum;
	    }
	}
    }

  /* Calculate new surface description. */

  /* Interpolate in second parameter direction. */

  dim = narr * nin1;

  s1891 (par2, tau, dim, nin2, 1, der2, TRUE, nknots2, &coef1, &nin2,
	 nik2, 0, 0, &kstat);
  if (kstat < 0) goto error;

  /* Interpolate in first parameter direction. */

  s1891 (par1, coef1, narr, nin1, nin2, der1, TRUE, nknots1, &coef2,
	 &nin1, nik1, 0, 0, &kstat);
  if (kstat < 0) goto error;

  /* OK */

  *nsurf = newSurf (nin1, nin2, nik1, nik2, nknots1, nknots2,
		    coef2, osurf->ikind, narr, 2);
  if (*nsurf == SISL_NULL) goto err171;

  goto out;

  /* Not enough memory. */

err101:
  *jstat = -101;
  s6err ("s1896", *jstat, kpos);
  goto out;

  /* Could not create surface, */

err171:
  *jstat = -171;
  s6err ("s1896", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1896", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1896", *jstat, kpos);
  goto out;

  /* Free pointers. */

out:
  if (coef1 != SISL_NULL)    freearray (coef1);
  if (val1 != SISL_NULL)     freearray (val1);
  if (val2 != SISL_NULL)     freearray (val2);
  if (par1 != SISL_NULL)     freearray (par1);
  if (par2 != SISL_NULL)     freearray (par2);
  if (der1 != SISL_NULL)     freearray (der1);
  if (der2 != SISL_NULL)     freearray (der2);
  if (normal != SISL_NULL)   freearray (normal);
  if (deriv != SISL_NULL)    freearray (deriv);
  if (tau != SISL_NULL)      freearray (tau);

  return;
}
