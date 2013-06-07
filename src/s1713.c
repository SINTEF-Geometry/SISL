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
 * $Id: s1713.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1713

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1713(SISLCurve *pc,double abeg,double aend,SISLCurve **rcnew,int *jstat)
#else
void s1713(pc,abeg,aend,rcnew,jstat)
     SISLCurve  *pc;
     double abeg;
     double aend;
     SISLCurve  **rcnew;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE     :To take one part of a closed B-spline curve and make a new
*              curve of the part. If the routine is used on an open curve
*              and aend<=abeg then the portion between aend and abeg is
*              removed by translating the last part of the new curve.
*              Works also over the seem for periodic curves.
*
*
* INPUT      : pc       - SISLCurve to take a part of.
*              abeg     - Start parameter-value of the curve part picked.
*              aend     - End parameter-value of the curve part picked.
*
*
*
* OUTPUT     : rcnew    - The new curve that is a part of the orginal curve.
*              jstat    - status messages
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
*              s1714     - Subdivide closed (periodic) crvs in two.
*              s1715     - Join two curves into one.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* MODIFIED BY : Ulf J. Krystad, SI, 92-01. Periodic crvs.
* MODIFIED BY : Arne Laksaa, SI, 92-09. Using D(N)EQUAL() insted of (!=)/==.
*
**********************************************************************/
{
  int kstat;          /* Local status variable.          */
  int kpos=0;         /* Position of error.              */
  double tbeg,tend;   /* The smaller and greater point.  */
  SISLCurve *q1=SISL_NULL; /* Pointer to new curve-object.    */
  SISLCurve *q2=SISL_NULL; /* Pointer to new curve-object.    */
  SISLCurve *q3=SISL_NULL; /* Pointer to new curve-object.    */
  SISLCurve *q4=SISL_NULL; /* Pointer to new curve-object.    */
  
  /* Check that we have a curve to pick a part of. */
  
  if (!pc) goto err150;
  
  /* Treating periodicity UJK, jan.92 ------- */
  if (pc->cuopen == SISL_CRV_PERIODIC)
  {
     s1714 (pc, abeg, aend, rcnew, &q1, jstat);
     if (q1) freeCurve(q1);q1=SISL_NULL;
     goto out;
  }

  /* Check that the intersection points is interior points. */
  
  if ((abeg < pc->et[0] && DNEQUAL(abeg,pc->et[0])) || 
      (abeg > pc->et[pc->in+pc->ik-1] && DNEQUAL(abeg,pc->et[pc->in+pc->ik-1])))
     goto err151;
  if ((aend < pc->et[0] && DNEQUAL(aend,pc->et[0])) || 
      (aend > pc->et[pc->in+pc->ik-1] && DNEQUAL(aend,pc->et[pc->in+pc->ik-1])))
    goto err151;
      
  /* Find the smaller and greater of the intersection points. */
  
  if (abeg<aend)
    {
      tbeg = abeg;
      tend = aend;
    } 
  else
    if (abeg>aend)
      {
	tbeg = aend;
	tend = abeg;
      }
  
  if (DEQUAL(abeg,aend))
  {
     /* In this case we have just one point to
	devide at. The result is two curves, q1 q1. */
     
     s1710(pc,abeg,&q1,&q3,&kstat);
     if (kstat<0 || kstat==2) goto err153;
  } 
  else
  {
     /* Devide into two at each point,
	we than have tree curves, q1 q2 q3.*/
     
     s1710(pc,tbeg,&q1,&q4,&kstat);
     if (kstat<0 || kstat==2) goto err153;
     
     s1710(q4,tend,&q2,&q3,&kstat);
     if (kstat<0 || kstat==2) goto err153;
     
     freeCurve(q4);  q4 = SISL_NULL;
  }
  
  /* If nessesary we have to join curve q3 and q1 to get the new curve.*/
  
  if (abeg > aend || DEQUAL(abeg,aend))
  {
     if (q2) 
     {
	freeCurve(q2);
	q2 = SISL_NULL;
     }
     if (!q1)
     {
	q2 = q3;
	q3 = SISL_NULL;
     }
     else if (!q3)
     {
	q2 = q1;
	q1 = SISL_NULL;
     }
     else
     {
	s1715(q3,q1,1,0,&q2,&kstat);
	if (kstat) goto err153;
     }
  }
  
  /* Updating output. */
  
  *rcnew = q2;
  *jstat = 0;
  goto out;
  
  /* Error. Subrutine error. */
  
 err153:
  *jstat = kstat;
  goto outfree;
  
  /* Error. No curve to pick a part of.  */
  
 err150:
  *jstat = -150;
  s6err("s1713",*jstat,kpos);
  goto out;
  
  /* Error. No part, abeg and aend has illegal values.  */
  
 err151:
  *jstat = -151;
  s6err("s1713",*jstat,kpos);
  goto out;
  
  /* Error in output. */
  
 outfree:
  if(q2) freeCurve(q2);
  
  /* Free local used memory. */
  
 out:
  if(q1) freeCurve(q1);
  if(q3) freeCurve(q3);
  if(q4) freeCurve(q4);
  return;
}

