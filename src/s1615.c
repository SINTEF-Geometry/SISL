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
 * $Id: s1615.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1615

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1615(double epoint[], int inbpnt, int idim, int eptyp[],int *jstat)
#else
void s1615(epoint, inbpnt, idim, eptyp, jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    eptyp[];
     int    *jstat;
#endif
/*
*************************************************************************
*
* PURPOSE: To test if the points lie on two branches of a hyperbola
*
* INPUT:
*        Epoint - Points/tangents describing the conic.
*        Inbpnt - No. of points/tangents in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        Eptyp  - Type indicator for points/tangents :
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Tangent to next point.
*                  4 - Tangent to prior point.
*
* Output:
*        Jstat - status messages:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method:
* 	If 3D points, the subroutine assumes that the third coordinate
* 	is zero, and that the first entry in Epoint is a point.
*	 Derivative vectors and vectors between the points are extracted
*	 and put into a temporary array. Then the cross product between
* 	adjacent ones of these are calculated. If both positive and negative
* 	cross products are found, the points lie on different branches
*	of a hyperbola.
*-
* Calls: s6err.
*
* Written by: A.M. Ytrehus, SI Oslo, Oct.91.
* After FORTRAN, (P1615), written by: T. Dokken  SI.
*****************************************************************
*/
{
  double svector[8];
  double *spoint = SISL_NULL;
  int ki, kk, kp, kki;
  int ktyp;
  double tcross;
  int kant = 4;
  int kpluss = 0;
  int kneg = 0;
  int kdim = 2;
  int kstore = 0;
  int kpos = 0;

  *jstat = 0;


  /* Allocate local array. */

  spoint = newarray (kdim * inbpnt, DOUBLE);
  if (spoint == SISL_NULL)
    goto err101;


  /* Kant is the no. of segments between the points (that is inbpnt-1).
     We can have max. 5 points. */

  if (inbpnt < 5)
    kant = inbpnt - 1;


  /* If the no. of segments is two or less, the points cannot
     belong to two different conic branches. (Ok, go out). */

  if (kant < 3)
    goto out;


  /* Store the positions of points and tangents (next/prior point
     minus/plus tangent/vector) in the working array spoit. */


  for (ki = 0; ki < inbpnt; ki++)
    {
      ktyp = eptyp[ki];
      kki = idim * ki;

      if (ktyp == 1 || ktyp == 2)
	{
	  /* Store this point directly. */

	  kk = kdim * kstore;

	  for (kp = 0; kp < kdim; kp++)
	    spoint[kk + kp] = epoint[kki + kp];

	  kstore++;
	}
      else if (ktyp == 3)
	{
	  /* Tangent to next point. Store the coordinates of
	     next point minus this one. */

	  kk = kdim * kstore;

	  for (kp = 0; kp < kdim; kp++)
	    spoint[kk + kp] = epoint[kki + idim + kp] - epoint[kki + kp];

	  kstore++;

	}
      else if (ktyp == 4)
	{
	  /* Tangent to prior point. Store the coordinates of
	     last point plus this one. */

	  kk = kdim * kstore;

	  for (kp = 0; kp < kdim; kp++)
	    spoint[kk + kp] = epoint[kki - idim + kp] + epoint[kki + kp];

	  kstore++;
	}
    }


  /* Run through the points/derivatives and put difference
     vectors into the svector-array. */


  kstore = 0;

  for (ki = 1; ki < inbpnt; ki++)
    {
      kk = kdim * kstore;
      kki = kdim * ki;

      for (kp = 0; kp < kdim; kp++)
	svector[kk + kp] = spoint[kki + kp] - spoint[kki - kdim + kp];

      kstore++;
    }


  /* If the cross product of adjacent vectors all have the same
     direction, then the points all lie on the same conic branch. */

  for (ki = 0; ki < kant - 1; ki++)
    {
      kk = kdim * ki;

      tcross = svector[kk] * svector[kk + 3]
	- svector[kk + 1] * svector[kk + 2];

      if (tcross > 0)
	kpluss++;
      if (tcross < 0)
	kneg++;
    }

  /* If both kpluss and kneg are nonzero, then the points
     lie on different branches. */

  if (kpluss > 0 && kneg > 0)
    *jstat = 1;

  goto out;

  /* Allocation error. */

err101:
  *jstat = -101;
  s6err ("s1615", *jstat, kpos);
  goto out;

out:

  /* Free local arrays. */

  if (spoint != SISL_NULL)
    freearray (spoint);

  return;
}
