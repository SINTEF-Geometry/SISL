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
 * $Id: s1919.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1919

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1919 (double et[], double prev[], double curr[], double deriv[],
       double follow[], int in, int ik, int idim, int iip, int iif,
       double ap, double ac, double af, int *jstat)
#else
void
s1919 (et, prev, curr, deriv, follow, in, ik, idim, iip, iif, ap, ac, af, jstat)
     double et[];
     double prev[];
     double curr[];
     double deriv[];
     double follow[];
     int in;
     int ik;
     int idim;
     int iip;
     int iif;
     double ap;
     double ac;
     double af;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To make the product of the tangent legth and the parametrization
*		of the curve epd as close as possible to the curve length.
*
* INPUT      :	et	- The original knot vector.
*
* OUTPUT     :jstat    - Status variable:
*                                               > 0     : warning
*                                               = 0     : ok
*                                               < 0     : error
*
* METHOD     :
*
* REFERENCES :	Fortran version by Tor Dokken, SI, 1991-07
*
* CALLS      : s1890, s1221, s1891, s6err.
*
* WRITTEN BY :	Trond Vidar Stensby, SI, 1991-07
*
*********************************************************************
*/
{
  SISLCurve *tcurve = SISL_NULL;	/* Temporary curve. */
  int knh;			/* Local number of verticews. */
  int left;			/* Used when calling s1221 */
  int pos, pos2;		/* Used for efficent adressing of arrays. */
  int ki, kj;			/* Loop control variables. */
  double tdist1;		/* Distance between previous and current. */
  double tdist2;		/* Distance between current and following. */
  double tlength;		/* Length of tangemt. */
  double tdiff;			/* Dummy. */
  double tfak;
  double tval1;			/* Adjusted parameter values. */
  double tval2;
  double *kpar = SISL_NULL;		/* Used when calculating parametrization. */
  int *kder = SISL_NULL;
  double *kpc = SISL_NULL;		/* Points on previous curve. */
  double *kdc = SISL_NULL;		/* Points on derivative curve. */
  double *kcc = SISL_NULL;		/* Points on current curve. */
  double *kfc = SISL_NULL;		/* Points on following curve. */
  double *epd = SISL_NULL;		/* Vertices of the derivative curve. */
  int kstat = 0;
  int kpos = 0;

  *jstat = 0;


  /* Test if legal input. */

  if (ik <= 1 || in <ik)
    goto err112;


  /* Produce parameter values. */

  s1890 (et, ik, in, &kpar, &kder, &kstat);
  if (kstat < 0)
    goto error;


  /* Allocate temporary arrays. */

  kpc = newarray (idim * in, DOUBLE);
  if (kpc == SISL_NULL)
    goto err101;
  kcc = newarray (idim * in, DOUBLE);
  if (kcc == SISL_NULL)
    goto err101;
  kdc = newarray (idim * in, DOUBLE);
  if (kdc == SISL_NULL)
    goto err101;
  kfc = newarray (idim * in, DOUBLE);
  if (kfc == SISL_NULL)
    goto err101;


  if (iip == 1)
    {
      /* Caculate interpolation points on previous curve. */

      tcurve = newCurve (in, ik, et, prev, 1, idim, 1);
      if (tcurve == SISL_NULL)
	goto err101;

      left = 0;
      for (ki = 0; ki < in; ki++)
	{
	  s1221 (tcurve, 0, kpar[ki], &left, &kpc[ki * idim], &kstat);
	  if (kstat < 0)
	    goto error;
	}
      if (tcurve != SISL_NULL)
	freeCurve (tcurve);
    }

  /* Caculate interpolation points on current curve. */

  tcurve = newCurve (in, ik, et, curr, 1, idim, 1);
  if (tcurve == SISL_NULL)
    goto err101;

  left = 0;
  for (ki = 0; ki < in; ki++)
    {
      s1221 (tcurve, 0, kpar[ki], &left, &kcc[ki * idim], &kstat);
      if (kstat < 0)
	goto error;
    }
  if (tcurve != SISL_NULL)
    freeCurve (tcurve);


  /* Caculate interpolation points on derivative curve. */

  tcurve = newCurve (in, ik, et, deriv, 1, idim, 1);
  if (tcurve == SISL_NULL)
    goto err101;

  left = 0;
  for (ki = 0; ki < in; ki++)
    {
      s1221 (tcurve, 0, kpar[ki], &left, &kdc[ki * idim], &kstat);
      if (kstat < 0)
	goto error;
    }
  if (tcurve != SISL_NULL)
    freeCurve (tcurve);


  if (iif == 1)
    {
      /* Caculate interpolation points on following curve. */

      tcurve = newCurve (in, ik, et, follow, 1, idim, 1);
      if (tcurve == SISL_NULL)
	goto err101;

      left = 0;
      for (ki = 0; ki < in; ki++)
	{
	  s1221 (tcurve, 0, kpar[ki], &left, &kfc[ki * idim], &kstat);
	  if (kstat < 0)
	    goto error;
	}
      if (tcurve != SISL_NULL)
	freeCurve (tcurve);
    }

  /* Adjust the points calculated on the derivative curve according to the
     distances between previous, current and following curve. */

  if (iip == 1)
    tval1 = fabs (ac - ap);
  if (iif == 1)
    tval2 = fabs (af - ac);

  pos = 0;
  for (ki = 0; ki < in; ki++)
    {
      tdist1 = (double) 0.0;
      tdist2 = (double) 0.0;
      tlength = (double) 0.0;

      for (kj = 0; kj < idim; kj++)
	{
	  if (iip == 1)
	    {
	      tdiff = kcc[pos] - kpc[pos];
	      tdist1 += tdiff * tdiff;
	    }
	  if (iif == 1)
	    {
	      tdiff = kfc[pos] - kcc[pos];
	      tdist2 += tdiff * tdiff;
	    }
	  tlength += kdc[pos] * kdc[pos];
	  pos++;
	}

      tdist1 = sqrt (tdist1);
      tdist2 = sqrt (tdist2);
      tlength = sqrt (tlength);
      if (tlength == (double) 0.0)
	tlength = (double) 1.0;


      /* Make scaling factor of tangent/derivative. */

      tfak = (double) 1.0;
      if (iip == 1 && iif != 1)
	tfak = tdist1 / (tval1 * tlength);
      else if (iip != 1 && iif == 1)
	tfak = tdist2 / (tval2 * tlength);
      else if (iip == 1 && iif == 1)
	tfak = min (tdist1 / (tval1 * tlength),
		    tdist2 / (tval2 * tlength));


      /* Find actual length of derivative curve at the parameter value.
	 Make new derivative length. */

      pos2 = pos - idim;
      for (kj = 0; kj < idim; kj++)
	{
	  kdc[pos2 + kj] *= tfak;
	}
    }

  /* Calculate new curve description. */

  s1891 (kpar, kdc, idim, in, 1, kder, 1, et, &epd, &knh, ik, 0, 0, &kstat);
  if (kstat < 0)
    goto error;

  memcopy (deriv, epd, idim * in, DOUBLE);


  /* OK */

  goto out;


  /* Allocation error. */

err101:
  *jstat = -101;
  s6err ("s1919", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1919", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1919", *jstat, kpos);
  goto out;

out:
  if (epd != SISL_NULL)
    freearray (epd);
  if (kpc != SISL_NULL)
    freearray (kpc);
  if (kcc != SISL_NULL)
    freearray (kcc);
  if (kdc != SISL_NULL)
    freearray (kdc);
  if (kfc != SISL_NULL)
    freearray (kfc);
  if (kpar != SISL_NULL)
    freearray (kpar);
  if (kder != SISL_NULL)
    freearray (kder);

  return;
}
