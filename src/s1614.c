/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
