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
 * $Id: s1712.c,v 1.3 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1712

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1712 (SISLCurve * pc, double abeg, double aend, SISLCurve ** rcnew, int *jstat)
#else
void
s1712 (pc, abeg, aend, rcnew, jstat)
     SISLCurve *pc;
     double abeg;
     double aend;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE     :To take one part of a open B-spline curve and make a new
*              curve of the part, if aend < abeg  the new curve is turned.
*
*
*
* INPUT      : pc     - SISLCurve to take a part of.
*              abeg     - Start parameter-value of the curve part picked.
*              aend     - End parameter-value of the curve part picked.
*
*
*
* OUTPUT     : rcnew      -The new curve that is a part of the orginal curve.
*              jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : freeCurve - Free space occupied by given curve-object.
*              s1710     - Divide a curve into two parts.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* Revised by : Tor Dokken, SI, 26-feb-1988.
*              Allow for converting from closed to open bases
* MODIFIED BY : Arne Laksaa, SI, 92-09. Using D(N)EQUAL() insted of (!=)/==.
*
**********************************************************************/
{
  int kstat;			/* Local status variable.           */
  int kpos = 0;			/* Position of error.               */
  double tbeg, tend;		/* The less and greater point.      */
  SISLCurve *q1 = SISL_NULL;		/* Pointer to new curve-object.     */
  SISLCurve *q2 = SISL_NULL;		/* Pointer to new curve-object.     */
  SISLCurve *q3 = SISL_NULL;		/* Pointer to new curve-object.     */
  int kturn = 0;

  /* Check that we have a curve to pick a part of. */

  if (!pc)
    goto err150;

  /* Check that the intersection points are interior points. */

  if ((abeg < pc->et[0] && DNEQUAL(abeg,pc->et[0])) ||
      (abeg > pc->et[pc->in+pc->ik-1] && DNEQUAL(abeg,pc->et[pc->in+pc->ik-1])))
    goto err151;
  if ((aend < pc->et[0] && DNEQUAL(aend,pc->et[0])) ||
      (aend > pc->et[pc->in+pc->ik-1] && DNEQUAL(aend,pc->et[pc->in+pc->ik-1])))
    goto err151;

  /* Check that abeg is not like aend. */

  if (DEQUAL(abeg,aend))
    goto err151;

  if (pc->cuopen == SISL_CRV_PERIODIC)
  {
      double delta = pc->et[pc->in] - pc->et[pc->ik - 1];

      if (abeg > aend)	kturn = 1;

      if (abeg < pc->et[pc->ik - 1] && DNEQUAL(abeg, pc->et[pc->ik - 1]))
	abeg += delta;
      if (abeg > pc->et[pc->in] || DEQUAL(abeg, pc->et[pc->in]))
	abeg -= delta;

      if (aend < pc->et[pc->ik - 1] && DNEQUAL(aend, pc->et[pc->ik - 1]))
	aend += delta;
      if (aend > pc->et[pc->in] && DNEQUAL(aend, pc->et[pc->in]))
	aend -= delta;

      if ((abeg > aend && !kturn) || (abeg < aend && kturn))
	 kturn = 1;
      else
	 kturn = 0;
  }

  /* Find the smaller and greater of the intersection points. */

  if (abeg < aend)
    {
      tbeg = abeg;
      tend = aend;
    }
  else
    {
      tbeg = aend;
      tend = abeg;
    }

  /* Divide into two at each point. The new curve is
     the middelmost curve.*/

  /* We start with dividing at the first point, and free the curv at
     left. */

  s1710 (pc, tbeg, &q1, &q2, &kstat);
  if (kstat < 0)
    goto err153;
  /* UJK, periodicity */
  if (kstat && q1 && !q2)
    {
      q2 = q1;
      q1 = SISL_NULL;
    }
  else if (q1)
    {
      freeCurve (q1);
      q1 = SISL_NULL;
    }

  /* Then we divide at the last point. The curve to left is the new curve.*/

  s1710 (q2, tend, &q1, &q3, &kstat);
  if (kstat < 0)
    goto err153;
  /* BOH, periodicity */
  if (kstat && !q1 && q3)
  {
    q1 = q3;
    q3 = SISL_NULL;
  }

  /* The curve is turned if nessesary. */

  if ((abeg > aend && !kturn) || (abeg < aend && kturn))
    s1706 (q1);

  /* Updating output. */

  *rcnew = q1;
  *jstat = 0;
  goto out;


/* Error. Error at low level function. */

err153:
  *jstat = kstat;
  goto outfree;


/* Error. No curve to pick a part of.  */

err150:
  *jstat = -150;
  s6err ("s1712", *jstat, kpos);
  goto out;


/* Error. No part, abeg and aend has illegal values.  */

err151:
  *jstat = -151;
  s6err ("s1712", *jstat, kpos);
  goto out;


/* Error in output. */

outfree:
  if (q1)
    freeCurve (q1);


/* Free local used memory. */

out:
   if (q2)
      freeCurve (q2);
   if (q3)
      freeCurve (q3);
   return;
}
