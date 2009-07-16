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
 * $Id: sh6iscnect.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6ISCONNECT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6isconnect(SISLIntpt *pt0, SISLIntpt *pt1, SISLIntpt *pt2)
#else
int sh6isconnect(pt0, pt1, pt2)
    SISLIntpt *pt0;
    SISLIntpt *pt1;
    SISLIntpt *pt2;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Boolean. Test if there exist a connection between
*              the two intersection points pt1 and pt2.
*
*
* INPUT      : pt0       - Pointer to previous Intpt, possibly SISL_NULL.
*              pt1       - Pointer to first Intpt.
*              pt2       - Pointer to second Intpt.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, 93-02.
*
*********************************************************************
*/
{
   int kstat = 0;   /* Status on wether a connection is found.  */
   int ki;          /* Counter.                                 */
   SISLIntpt *qt;   /* Intersection point.                      */
   int been_here = -199;

   /* Test if the points are equal.  */

   if (pt1 == pt2) return 1;
   
   /* UJK, aug 93, oo-loop in sh6isconn. */
   if (pt1->marker == been_here) return 0;
   pt1->marker = been_here;

   /* Traverse the intersection points connected to pt1.  */
   
   for (ki=0; ki<pt1->no_of_curves; ki++)
   {
      qt = sh6getnext(pt1,ki);
      if (qt == pt0) continue;
      
      kstat = sh6isconnect(pt1,qt,pt2);
      
      if (kstat) return 1;
   }
   
   /* No connection is found.  */
   
   return 0;
}
