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
