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
 * $Id: s1614.c,v 1.2 1994-12-07 12:05:37 pfu Exp $
 *
 */


#define S1614

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1614(double epoint[], int inbpnt, int idim, int eptyp[],
	   double spoint[], int *jnbpnt, int sptyp[], int *jstat)

#else
void s1614(epoint, inbpnt, idim, eptyp,spoint, jnbpnt, sptyp, jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    eptyp[];
     double spoint[];
     int    *jnbpnt;
     int    sptyp[];
     int    *jstat;
#endif
/*
*************************************************************************
*
* PURPOSE: To check the point description vector to ensure that the
*          points returned from this routine can be used directly in the
*          conversion of conics to B-spline representation.
*
* INPUT:
*        Epoint - The points/tangents describing the conic.
*        Inbpnt - No. of points/tangents in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        Eptyp  - Type indicator for points/tangents :
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Tangent to next point.
*                  4 - Tangent to prior point.
*
* Output:

*        Spoint - The corrected points.
*        Jnbpnt - No. of corrected points.
*        Sptyp  - Type indicator for corrected points.
*        Jstat  - status messages:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method:
* 	The point description is run through, to check if some points are
*	 assigned two tangents. If two tangents are detected, the first is
* 	used and the second discarded. In addition, the first points
* 	interpolation-condition detected is put in the start of the output
* 	array, and the last points interpolation-condition is put in the end
* 	of the output array.
*-
* Calls: s6err.
*
* Written by: A.M. Ytrehus, SI  Oslo Oct.91.
* After FORTRAN (P1614) written by: T. Dokken  SI.
* Revised by: J. Kaasa, SI, Aug. 92 (Fixed bugs in connection with
*             type handling).
*****************************************************************
*/
{
  int ki, kp, kk, kki;
  int khelp;
  int ktell = 0;
  int kflag = 0;
  int knbpnt = 0;
  int ktyp = 1;
  int ktypm1;
  double tdum;
  int kpos = 0;

  *jstat = 0;


  /* Run through and test the type of the data points. Copy
     approved points and tangents into output variables. */

  for (ki = 0; ki < inbpnt; ki++)
    {
      ktypm1 = ktyp;
      ktyp = eptyp[ki];

      if (!((ktyp < 1 || ktyp > 4) ||
	    (knbpnt == 0 && ktyp == 4) ||
	    (ktyp == 3 && ktypm1 == 3) ||
	    (ktyp == 4 && ktypm1 == 4) ||
	    (ktyp == 4 && ktypm1 == 3)))
	{
	  /* Save this point or tangent in spoint and sptyp. */

	  sptyp[knbpnt] = ktyp;

	  kki = idim * ki;
	  kk = idim * knbpnt;

	  for (kp = 0; kp < idim; kp++)
	    spoint[kk + kp] = epoint[kki + kp];
	  knbpnt++;
	}
    }

  /* Remove last point if the type is 3 (which means tangent to next point. */

  if (ktyp == 3)
    knbpnt--;


  /* We cannot make use of more than five points.
     Choose the first five points. */

  if (knbpnt > 5)
    {
      knbpnt = 5;

      /* Check the type of (new) last point. */

      if (sptyp[knbpnt - 1] == 3)
	{
	  /* Cut out this point, and use point no. 6 instead. */

	  sptyp[knbpnt - 1] = sptyp[knbpnt];

	  kk = idim * (knbpnt - 1);

	  for (kp = 0; kp < idim; kp++)
	    spoint[kk + kp] = spoint[kk + idim + kp];
	}
    }


  /* Ensure that first entry is a point. */

  ktyp = sptyp[0];

  if (ktyp > 2)
    {
      /* Interchange the two first conditions. */

      kflag = 1;

      sptyp[0] = 1;
      sptyp[1] = 4;

      for (kp = 0; kp < idim; kp++)
	{
	  tdum = spoint[kp];
	  spoint[kp] = spoint[kp + idim];
	  spoint[kp + idim] = tdum;
	}
    }

  /* Ensure that last entry is a point. */

  ktyp = sptyp[knbpnt - 1];

  if (ktyp > 2)
    {

      /* Interchange the two last conditions. */

      kflag = 1;

      sptyp[knbpnt - 1] = 1;
      sptyp[knbpnt - 2] = 3;

      kk = idim * (knbpnt - 1);

      for (kp = 0; kp < idim; kp++)
	{
	  tdum = spoint[kk + kp];
	  spoint[kk + kp] = spoint[kk - idim + kp];
	  spoint[kk - idim + kp] = tdum;
	}
    }

  /* Count no. of points. If less than two, we cannot even
     create a straight line. */

  for (ki = 0; ki < knbpnt; ki++)
    {
      ktyp = sptyp[ki];
      if (ktyp < 3)
	ktell++;
    }

  if (ktell <= 1)
    goto err181;


  /* Run through and remove double tangents that may be introduced by
     the interchange of end-point/end-tangents. */

  ktyp = 1;

  if (kflag == 1)
    {
      khelp = knbpnt;

      ktypm1 = ktyp;

      for (ki = 0; ki < khelp; ki++)
	{
	  ktyp = sptyp[ki];

	  if ((ktyp == 3 && ktypm1 == 3) ||
	      (ktyp == 4 && ktypm1 == 4) ||
	      (ktyp == 4 && ktypm1 == 3))
	    {
	      /* Remove last of doble tangent. */

	      knbpnt--;

	      for (kk = ki; kk < knbpnt; kk++)
		{
		  sptyp[kk] = sptyp[kk + 1];

		  kki = idim * kk;

		  for (kp = 0; kp < idim; kp++)
		    spoint[kki + kp] = spoint[kki + idim + kp];
		}
	    }
	}
    }


  *jnbpnt = knbpnt;

  goto out;


  /* Error in input, less than two points given in problem formulation. */

err181:
  *jstat = -181;
  s6err ("s1614", *jstat, kpos);
  goto out;

out:

  return;
}
