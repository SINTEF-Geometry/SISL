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
 * $Id: s1619.c,v 1.2 1995-01-26 08:44:48 pfu Exp $
 *
 */


#define S1619

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1619(double epoint[], int inbpnt, int idim, int eptyp[],
	   double econic[], int ityp, double etang[],
	   double *ashape, int *jstat)
#else
void s1619(epoint, inbpnt, idim, eptyp, econic, ityp,etang, ashape, jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    eptyp[];
     double econic[];
     int    ityp;
     double etang[];
     double *ashape;
     int    *jstat;
#endif
/*
*************************************************************************
*
* PURPOSE: To Produce the intersection point of the start- and end-tangent,
*          the midpoint of the line from the start to the endpoint
*          and a "shape factor".
* INPUT:
*        Epoint - The points/tangents describing the conic.
*        Inbpnt - No. of points/tangents in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        Eptyp  - Type indicator for the points/tangents :
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Tangent to next point.
*                  4 - Tangent to prior point.
*        Econic - The conic coefficients of the points in the Epoint array.
*        Ityp   - Type of conic to be produced.
*                  1 - Straight line.
*                  2 - Ellipse.
*                  3 - Parabola.
*                  4 - Hyperbola.
*
* Output:
*        Etang   -The intersection point of start and end tangents.
*        Ashape - Shape factor of the conic.
*        Jstat  - status variable:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method:
*        The points are calculated using the equation of lines and
*        the equation of the conic.
*-
* Calls: No.
*
* Written by: A.M. Ytrehus, si Oslo, Oct.91.
* After FORTRAN, (P1619), written by: T. Dokken  SI.
*****************************************************************
*/
{
  double ta11, ta12, ta13, ta22, ta23, ta33;
  double tx, ty, tx1, ty1, tx2, ty2, axt, ayt;
  double tv1, tv2, tvx, tvy, tnx, tny, tix, tiy;
  double ta, tb, tc, te1, te2, te3;
  double tys1, tys2, txs1, txs2;
  double as, ts1, ts2;
  double tlong, tlength;
  double tdum;
  int ktyp;

  int ki, kk;
  int kki = 0;

  int krem = 0;
  int ksam = 0;

  *jstat = 0;


  /* Initiate variables. */

  kk = idim * (inbpnt - 1);

  ta11 = econic[0];
  ta12 = econic[1];
  ta13 = econic[3];		/* ? */
  ta22 = econic[2];		/* ? */
  ta23 = econic[4];
  ta33 = econic[5];

  tx = (epoint[0] + epoint[kk]) / ((double) 2.0);
  ty = (epoint[1] + epoint[kk + 1]) / ((double) 2.0);


  /* Calculate tangent vectors in start- and end-point. */

  tx1 = -ta22 * epoint[1] - ta12 * epoint[0] - ta23;
  ty1 = ta11 * epoint[0] + ta12 * epoint[1] + ta13;
  tlong = sqrt (tx1 * tx1 + ty1 * ty1);
  tx1 = tx1 / tlong;
  ty1 = ty1 / tlong;

  tx2 = -ta22 * epoint[kk + 1] - ta12 * epoint[kk] - ta23;
  ty2 = ta11 * epoint[kk] + ta12 * epoint[kk + 1] + ta13;
  tlong = sqrt (tx2 * tx2 + ty2 * ty2);
  tx2 = tx2 / tlong;
  ty2 = ty2 / tlong;


  /* Calculate intersection point of tangents (axt,ayt). */

  tdum = ty1 * tx2 - ty2 * tx1;


  /* Remember if parallel tangents. */

  krem = 0;

  if (fabs (tdum) <= REL_PAR_RES)
    {
      /* Parallel tangents. */

      krem = 1;
      tv1 = tx1;
      tv2 = ty1;
    }
  else
    {
      axt = (epoint[0] * ty1 * tx2 - epoint[kk] * ty2 * tx1 +
	     tx1 * tx2 * (epoint[kk + 1] - epoint[1])) / tdum;
      ayt = (ty1 * ty2 * (epoint[0] - epoint[kk]) -
	     epoint[1] * tx1 * ty2 + epoint[kk + 1] * tx2 * ty1) / tdum;


      /* Calculate the intersection point(s) between the conic
	 and the line from (axt,ayt) to (tx,ty). */

      tv1 = tx - axt;
      tv2 = ty - ayt;
      tlong = sqrt (tv1 * tv1 + tv2 * tv2);
      tv1 = tv1 / tlong;
      tv2 = tv2 / tlong;
    }

  ta = -tv2;
  tb = tv1;
  tc = -tx * ta - ty * tb;

  if (fabs (ta) >= fabs (tb))
    {
      /* We put the straight line into the equation
	 of the conic and get: te1*y*y + te2*y + te3 = 0. */

      te1 = ta11 * tb * tb / (ta * ta) - ((double) 2.0) * ta12 * tb / ta + ta22;
      te2 = ((double) 2.0) * tb * tc * ta11 / (ta * ta) - ((double) 2.0) * ta12 * tc / ta
	- ((double) 2.0) * ta13 * tb / ta + ((double) 2.0) * ta23;
      te3 = ta11 * tc * tc / (ta * ta) - ((double) 2.0) * ta13 * tc / ta + ta33;
      tdum = te2 * te2 - ((double) 4.0) * te1 * te3;


      /* If tdum < 0.0 no intersection point. Produce straight line. */

      if (tdum < 0.0)
	{
	  /* Produce straight line; */

	  *jstat = 1;
	  goto out;
	}

      tdum = sqrt (tdum);
      tys1 = (-te2 - tdum) / (((double) 2.0) * te1);
      tys2 = (-te2 + tdum) / (((double) 2.0) * te1);
      txs1 = -(tb * tys1 + tc) / ta;
      txs2 = -(tb * tys2 + tc) / ta;
    }
  else
    {
      /* ta equal to zero. We put the straight line into the conic
	 equation and get: te1*x*x + te2*x + te3 = 0; */

      te1 = ta11 - ((double) 2.0) * ta12 * ta / tb + ta22 * ta * ta / (tb * tb);
      te2 = ((double) 2.0) * ta * tc * ta22 / (tb * tb) - ((double) 2.0) * ta12 * tc / tb
	- ((double) 2.0) * ta23 * ta / tb + ((double) 2.0) * ta13;
      te3 = ta22 * tc * tc / (tb * tb) - ((double) 2.0) * ta23 * tc / tb + ta33;
      tdum = te2 * te2 - ((double) 4.0) * te1 * te3;


      /* If tdum < 0.0, no intersection point. Produce straight line. */

      if (tdum < (double) 0.0)
	{
	  /* Produce straight line. */

	  *jstat = 1;
	  goto out;
	}

      tdum = sqrt (tdum);
      txs1 = (-te2 - tdum) / (((double) 2.0) * te1);
      txs2 = (-te2 + tdum) / (((double) 2.0) * te1);
      tys1 = -(ta * txs1 + tc) / tb;
      tys2 = -(ta * txs2 + tc) / tb;
    }

  /* Calculate shape factors, and find the shape factor to be used.
     Special treatment if parallel tangents. */

  if (krem == 1)
    {
      /* Parallel tangents. One of the intersection points is to be
	 used tshape tangent intersection point. We try the first. */

      axt = txs1;
      ayt = tys1;
    }

  tlength = (tx - axt) * (tx - axt) + (ty - ayt) * (ty - ayt);


  /* If the tangent point and the midpoint
     are the same, produce straight line. */

  if (tlength < (double) 0.0)
    {
      /* Produce straight line. */

      *jstat = 1;
      goto out;
    }


  ts1 = ((txs1 - tx) * (axt - tx) + (tys1 - ty) * (ayt - ty)) / tlength;
  ts2 = ((txs2 - tx) * (axt - tx) + (tys2 - ty) * (ayt - ty)) / tlength;

  if (ts1 >= (double) 1.0 && ts2 >= (double) 1.0)
    {
      /* Produce straight line. */

      *jstat = 1;
      goto out;
    }

  /* Treat parabolas and hyperbolas. */

  if (ityp >= 3)
    {
      /* Hyperbola or parabola. */

      as = ts1;

      if (as >= (double) 1.0)
	as = ts2;

      if (as >= (double) 1.0)
	{
	  /* Produce straight line. */

	  *jstat = 1;
	  goto out;
	}

      ts1 = as;
      ts2 = as;
    }

  /* The variable Ksam is assigned the value 1 if the elliptic arc
     search is on the same side as the tangent point, else 0.
     Search for internal point. */

  for (ki = 1; ki < inbpnt - 1; ki++)
    {
      ktyp = eptyp[ki];
      if (ktyp < 3)
	break;
    }

  /* If we are here and ktyp < 3, internal point found.
     Create normal-vector for line from start-point to endpoint. */

  if (ktyp < 3)
    {
      /* Internal point. */

      tnx = -epoint[1] + epoint[kk + 1];
      tny = epoint[0] - epoint[kk];


      /* Vector from end-point to internal point. */

      kki = idim * ki;
      tix = epoint[kki] - epoint[0];
      tiy = epoint[kki + 1] - epoint[1];


      /* Vector from start-point to tangent point. */

      tvx = axt - epoint[0];
      tvy = ayt - epoint[1];

      ksam = 1;
      if ((tnx * tix + tny * tiy) * (tnx * tvx + tny * tvy) < (double) 0.0)
	ksam = 0;
    }
  else if (ktyp > 3)
    {
      /* No internal points. Look for start-vector. */

      ktyp = eptyp[1];

      if (ktyp == 4)
	{
	  /* Start-tangent. Vector from start-point to tangent point. */

	  tvx = axt - epoint[0];
	  tvy = ayt - epoint[1];
	  ksam = 1;
	  if ((epoint[idim] * tvx + epoint[idim + 1] * tvy) < (double) 0.0)
	    ksam = 0;
	}
      else
	{
	  /* No internal point or start-tangent. Look for end-tangent.
             If no end-tangent, produce a straight line. */

	  ktyp = eptyp[inbpnt - 2];

	  if (ktyp != 3)
	    {
	      /* Produce straight line. */

	      *jstat = 1;
	      goto out;
	    }

	  /* Vector from tangent intersection point to end-point. */

	  tvx = epoint[kk] - axt;
	  tvy = epoint[kk + 1] - ayt;
	  ksam = 1;
	  if ((epoint[kki - idim] * tvx + epoint[kki - idim + 1] * tvy)
	      < (double) 0.0)
	    ksam = 0;
	}
    }

  /* Find which branch of the ellipse to use. */

  if (ksam == 1)
    {
      /* Special treatment if parallel tangents. We will use
	 the found tangent intersection point. */

      if (krem == 0)
	{

	  /* Use branch of same side as tangent intersection point.
	     Positive shape factor. */

	  as = ts1;
	  if (as < (double) 0.0)
	    as = ts2;
	  if (as < (double) 0.0)
	    {
	      /* Create straight line. */

	      *jstat = 0;
	      goto out;
	    }
	  else
	    {
	      /* Special treatment if parallel tangents. Use
		 (txs2,tys2) as tangent intersection point. */

	      if (krem == 1)
		{
		  /* Parallel tangents. */

		  axt = txs2;
		  ayt = tys2;
		}
	      else
		{
		  /* use branch on opposite side as tangent intersection
		     point. Negative shape factor. */

		  as = ts1;
		  if (as > (double) 0.0)
		    as = ts2;
		  if (as > (double) 0.0)
		    {
		      /* Create straight line. */

		      *jstat = 1;
		      goto out;
		    }
		}
	    }
	}
    }

  /* Shape factor for ellipse found. */

  /* If parallel tangent, subtract the middle of the line from the
     start-point to the end-point from the tangent-point. */

  if (krem == 1)
    {
      axt = axt - (epoint[0] + epoint[kk]) / ((double) 2.0);
      ayt = ayt - (epoint[1] + epoint[kk + 1]) / ((double) 2.0);
      as = (double) 0.0;
    }

  *ashape = as;
  etang[0] = axt;
  etang[1] = ayt;
  if (idim == 3)
    etang[2] = (double) 0.0;

  goto out;

out:

  return;
}
