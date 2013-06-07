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
 * $Id: s1716.c,v 1.3 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1716

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1716(SISLCurve *pc1,SISLCurve *pc2,
	   double aeps,SISLCurve **rcnew,int *jstat)
#else
void s1716(pc1,pc2,aeps,rcnew,jstat)
     SISLCurve  *pc1;
     SISLCurve  *pc2;
     double aeps;
     SISLCurve  **rcnew;
     int    *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To join the closest of one B-spline curve with the closest
*              end of another B-spline curve by translating the second curve
*              if the curves are closer to each other than aeps.
*              If pc1 is to be joined at the start the direction of the
*              curve is turned, and if pc2 is to be joined at the end
*              the direction of this curve is turned. This means that
*              pc1 always is at the beginning at the new curve.
*              If aeps is to small to any joining a SISL_NULL pointer is returned.
*
*
* INPUT      : pc1     - First curve to join.
*              pc2     - Second curv to join.
*              aeps    - The curves is to be joined if aeps is greater or
*                        like the distance between the closest end.
*                        If aeps is negativ the curve automaticaly is joined.
*
*
*
* OUTPUT     : rcnew   - The new joined curve.
*              jstat   - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : First we are finding the smallest distens between
*              curves. If this distens is smaller than aeps the curves
*              are joining at the closest ends
*
*
*
* REFERENCES :
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              s1715.c   - Join two curves at specified ends.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* Revised by : Paal Fugelli. SINTEF, Oslo, Norway, Sept. 1994.  Fixed array
*              bounds over-running.
*
**********************************************************************/
{
  int kstat;              /* Local status variable.                  */
  int kpos=0;             /* Position of error.                      */
  int knr;                /* Number to mark type of junction.        */
  int km11=0,km12=0;      /* Knot mutiplicety at the ends of
			     the first curve.                        */
  int km21=0,km22=0;      /* Knot mutiplicety at the ends of
			     the second curve.                       */
  int kdim;               /* Dimensjon of the space in whice curves
			     lies.                                   */
  int kk1=pc1->ik;        /* The order of the first old curve.       */
  int kk2=pc2->ik;        /* The order of the second old curve.      */
  int kn1=pc1->in;        /* Number of vertices in the old curves.   */
  int kn2=pc2->in;        /* Number of vertices in the old curves.   */
  int ki,kj1,kj2;      /* Control variable in loop, and others.   */
  double t1,tdel,tdelmin; /* The translation of the knots to the
			     second curve.                           */
  SISLCurve *qc=SISL_NULL;         /* Pointer to the new curve-object.        */

  /* Check that we have curves to join. */

  if (!pc1 || !pc2) goto err150;

  /* Check that The curves is in the same room, have the same kdim. */

  if (pc1->idim != pc2->idim) goto err106;
  else kdim = pc1->idim;

  /* Finding the knot multiplicity at the ends. */

  while (pc1->et[km11] == *pc1->et) km11++;
  while (pc1->et[kn1+kk1-1-km12] == pc1->et[kn1+kk1-1]) km12++;
  while (pc2->et[km21] == *pc2->et) km21++;
  while (pc2->et[kn2+kk2-1-km22] == pc2->et[kn2+kk2-1]) km22++;

  /* Finding the smallest distance between the ends. */

  /* First we compute the square distance between the start of both
     curves, Then we mark this to be the shortest. */

  for (tdel=DZERO,ki=0; ki<kdim; ki++)
    {
      if (km11<kk1)  t1 = DZERO;
      else           t1 = pc1->ecoef[kdim*(km11-kk1)+ki];
      if (km21>=kk2) t1 -= pc2->ecoef[kdim*(km21-kk2)+ki];
      tdel += t1*t1;
    }
  tdelmin = tdel;
  knr = 0;

  /* The start of the first curve and the end of the second curve. */

  for (tdel=DZERO,ki=0; ki<kdim; ki++)
    {
      if (km11<kk1)  t1 = DZERO;
      else           t1 = pc1->ecoef[kdim*(km11-kk1)+ki];
      if (km22>=kk2) t1 -= pc2->ecoef[kdim*(kn2-1-km22+kk2)+ki];
      tdel += t1*t1;
    }
  if (tdel<tdelmin)
    {
      tdelmin = tdel;
      knr = 1;
    }

  /* The end of the first curve and the start of the second curve. */

  for (tdel=DZERO,ki=0; ki<kdim; ki++)
    {
      if (km12<kk1)  t1 = DZERO;
      else           t1 = pc1->ecoef[kdim*(kn1-1-km12+kk1)+ki];
      if (km21>=kk2) t1 -= pc2->ecoef[kdim*(km21-kk2)+ki];
      tdel += t1*t1;
    }
  if (tdel<tdelmin)
    {
      tdelmin = tdel;
      knr = 2;
    }

  /* The end of the first curve and the end of the second curve. */

  for (tdel=DZERO,ki=0; ki<kdim; ki++)
    {
      if (km12<kk1)  t1 = DZERO;
      else           t1 = pc1->ecoef[kdim*(kn1-1-km12+kk1)+ki];
      if (km22>=kk2) t1 -= pc2->ecoef[kdim*(kn2-1-km22+kk2)+ki];
      tdel += t1*t1;
    }
  if (tdel<tdelmin)
    {
      tdelmin = tdel;
      knr = 3;
    }


  if (aeps < DZERO || aeps >= sqrt(tdelmin))
    {
      /* We mark what ends we are going to use in the junction.
	 and then call a function to join the curves. */

      if (knr<2) kj1 = 0;
      else        kj1 = 1;

      if (knr==0 || knr==2) kj2 = 0;
      else                kj2 = 1;


      s1715(pc1,pc2,kj1,kj2,&qc,&kstat);
      if (kstat) goto err153;
    } else
      {
	/* Aeps was to small We just have to return SISL_NULL. */

	qc = SISL_NULL;
      }


  /* Updating output. */

  *rcnew = qc;
  *jstat = 0;
  goto out;


  /* Error. Subrutine error. */

 err153:
  *jstat = kstat;
  goto outfree;


  /* Error. No curve to subdevice.  */

 err150:
  *jstat = -150;
  s6err("s1716",*jstat,kpos);
  goto out;


  /* Error. Different dimensjon of the room.  */

 err106:
  *jstat = -106;
  s6err("s1716",*jstat,kpos);
  goto out;


 outfree:
  if(qc) freeCurve(qc);


 out:
  return;
}
