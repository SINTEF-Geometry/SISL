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
 * $Id: s1714.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1714

#define SISL_CRV_PERIODIC -1
#define SISL_CRV_OPEN 1
#define SISL_CRV_CLOSED 0

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1714 (SISLCurve * pc, double apar1, double apar2, SISLCurve ** rcnew1, SISLCurve ** rcnew2, int *jstat)
#else
void 
s1714 (pc, apar1, apar2, rcnew1, rcnew2, jstat)
     SISLCurve *pc;
     double apar1;
     double apar2;
     SISLCurve **rcnew1;
     SISLCurve **rcnew2;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE   :Divide a closed B-spline curve in two parts at two
*            spesified values. The first curve start at apar1.
*            If the curve is open the last part of the curve
*            that is going through the gap between end and start
*            of the orginal curve is translating .
*            NOTE: For periodic curves when apar2 < apar1,
*                  the second curve will hav startparameter
*                  value shifted one period.
*
*
* INPUT      : pc        - SISLCurve to divide
*              apar1     - Start parameter-value of the first new curve.
*              apar2     - Start parameter-value of the second new curve.
*
*
*
* OUTPUT     : rcnew1    - The first new curve.
*              rcnew2    - The second new curve.
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
*              s1712     - Pick a part of an open curve.
*              s1713     - Pick a part of a closed curve.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* MODIFIED BY : Ulf J. Krystad, SI, 92-01. Periodic crvs.
* MODIFIED BY : Arne Laksaa, SI, 92-09. Using D(N)EQUAL() insted of (!=)/==.
*
**********************************************************************/
{
  int kstat;			/* Local status variable.        */
  int kpos = 0;			/* Position of error.            */
  SISLCurve *q1 = SISL_NULL;		/* Pointer to new curve-object.  */
  SISLCurve *q2 = SISL_NULL;		/* Pointer to new curve-object.  */

  /* Check that we have a curve to devide. */

  if (!pc)
    goto err150;

  /* Check that apar1 is not equal apar2. */

  if (DEQUAL (apar1, apar2))
    goto err151;

  /* Treating periodicity UJK, jan.92 and later ALA, sep.92 ------- */
  if (pc->cuopen == SISL_CRV_PERIODIC)
    {
      double delta = pc->et[pc->in] - pc->et[pc->ik - 1];
      
      while(apar1 < pc->et[pc->ik - 1] && DNEQUAL(apar1, pc->et[pc->ik - 1]))
	 apar1 += delta;
      while(apar1 > pc->et[pc->in] || DEQUAL(apar1, pc->et[pc->in]))
	 apar1 -= delta;

      while (apar2 < apar1 || DEQUAL(apar2, apar1))
	apar2 += delta;
      
      while (apar2 > (apar1+delta) && DNEQUAL(apar2, (apar1+delta)))
	apar2 -= delta;

      /* Shift startpoint of curve (and make it ordinary closed )*/
      s1710 (pc, apar1, &q1, &q2, &kstat);
      if (kstat < 0)
	goto err153;

      if (q2)
	freeCurve (q2);
      q2 = SISL_NULL;

      /* Split into two */
      s1710 (q1, apar2, rcnew1, rcnew2, &kstat);
      if (kstat < 0)
	goto err153;

	if (q1)
	freeCurve (q1);
      q1 = SISL_NULL;


      *jstat = 0;
      goto out;
    }
  /* End of treating periodicity UJK, jan.92 ------- */


    /* Divide the curve into two at each point.
     Join the two end curves at each end.*/

if (apar1 < apar2)
  {
    s1712 (pc, apar1, apar2, &q1, &kstat);
    if (kstat)
      goto err153;

    s1713 (pc, apar2, apar1, &q2, &kstat);
    if (kstat)
      goto err153;
  }

else
  {
    s1712 (pc, apar2, apar1, &q2, &kstat);
    if (kstat)
      goto err153;

    s1713 (pc, apar1, apar2, &q1, &kstat);
    if (kstat)
      goto err153;
  }

 /* Updating output. */

*rcnew1 = q1;
*rcnew2 = q2;
*jstat = 0;
goto out;


 /* Error. Subrutine error. */

err153:
*jstat = kstat;
goto outfree;


 /* Error. No curve to pick a part of.  */

err150:
*jstat = -150;
s6err ("s1714", *jstat, kpos);
goto out;


 /* Error. No part, apar1 and apar2 has illegal values.  */

err151:
*jstat = -151;
s6err ("s1714", *jstat, kpos);
goto out;


 /* Error in output. */

outfree:
if (q1)
  freeCurve (q1);
if (q2)
  freeCurve (q2);

out:
return;
}
