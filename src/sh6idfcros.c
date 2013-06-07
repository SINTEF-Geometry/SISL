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
 * $Id: sh6idfcros.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6IDFCROSS

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6idfcross(SISLIntdat *pintdat, SISLIntpt *vcross[], int *jncross,
		  int ipar1, int ipar2, int *jstat)
#else
void sh6idfcross(pintdat, vcross, jncross, ipar1, ipar2, jstat)
   SISLIntdat *pintdat;
   SISLIntpt  *vcross[];
   int        *jncross;
   int        ipar1;
   int        ipar2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given intersection data. Check if there exist 4
*              intersection points where where the parameters
*              corresponding to each object in the intersection
*              are pairwise equal.
*
*
* INPUT      : pintdat  - Intersection data.
*              vcross   - Intersection points found so far.
*              jncross  - Number of intersection points found so far.
*              ipar1    - Number of parameter directions of first object.
*              ipar2    - Number of parameter directions of second object.
*
* OUTPUT     : vcross   - Intersection points found.
*              jncross  - Number of intersection points found.
*              jstat    - Status
*                         jstat = 0   => No set of points is found.
*                         jstat = 1   => Successful. 4 points are found.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*
* CALLS      : s6dist    -  Distance between points.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 12.92.
*
*********************************************************************
*/
{ 
   int ki,kj;       /* Counters.                               */
   int kpt;         /* Index of last intersection point found. */
   int kpar1;       /* Start index of current parameter set.   */
   int kpar2;       /* Number of parameter in current set.     */
   double tdist;    /* Distance between parameter points.      */
   SISLIntpt *pt;   /* Current intersection point.             */
   SISLIntpt *qnext; /* Next point to find.                    */
   
 
   /* Test if there is 4 points in pintdat.  */
   
   if (pintdat->ipoint < 4)
   {
      /* No possibility of cross intersections. */
      
      *jstat = 0;
      return;
   }
   
   /* Test if a set of cross intersections is found. */
   
   if (*jncross == 4)
   {
      /* Test if the second parameter set of the last intersection point
	 found is equal to that of the first point.         */
      
      tdist = s6dist(vcross[0]->epar+ipar1,vcross[3]->epar+ipar1,ipar2);
      if (DEQUAL(tdist+(double)1.0,(double)1.0))
	 /* The set of points is found.  */
	 
	 *jstat = 1;
      else
	 *jstat = 0;
      
      return;
   }
   
   /* Prepare for a search for the next point in the set.  */
   
   kpt = (*jncross) - 1;
   pt = vcross[kpt];
   kpar1 = (kpt % 2 == 0) ? 0 : ipar1;
   kpar2 = (kpt % 2 == 0) ? ipar1 : ipar2;
   
   /* Traverse the intersection points to find a point that has got
      one parameter set equal to the current one.  */
   
   for (ki=0; ki<pintdat->ipoint; ki++)
   {
      qnext = pintdat->vpoint[ki];
      
      /* Check if the point is found already.  */
      
      for (kj=0; kj<=kpt; kj++)
	 if (qnext == vcross[kj]) break;
      if (kj <= kpt) continue;
      
      /* Check if the next point belongs to the wanted set. */
      
      tdist = s6dist(qnext->epar+kpar1,pt->epar+kpar1,kpar2);
      if (DEQUAL(tdist+(double)1.0,(double)1.0))
      {
	 /* A point is found.  */
	 
	 kpt++;
	 vcross[kpt] = qnext;
	 (*jncross)++;
	 
	 /* Find next point.  */
	 
	 sh6idfcross(pintdat,vcross,jncross,ipar1,ipar2,jstat);
	 if (*jstat == 1) return;  /* The entire set is found.  */
	 
	 (*jncross)--;
	 kpt--;
      }
   }
   
   /* No set of cross intersections exist.  */
   
   *jstat = 0;
   return;
}
