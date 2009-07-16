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
 * $Id: sh6setdir.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6SETDIR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6setdir(SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6setdir(pt1,pt2,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set the direction of the curve through pt1 and pt2 such
*              that pt1 points to pt2.
*              If they are not connected, give error.
*
*
* INPUT      : pt1      - Pointer to first Intpt.
*              pt2      - Pointer to second Intpt.
*              jstat    - Error flag.
*                        jstat =  0  => Successful
*                        jstat = -1  => Points are not connected.
*                        jstat = -2  => Error in subfunction.
*
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. July 91.
*********************************************************************
*/
{
   int kstat;         /* error flag. */
   int index1,index2; /* dummy indices.           */

   *jstat = 0;

   /* Check if pt1 and pt2 are already connected. */

   sh6getlist(pt1,pt2,&index1,&index2,&kstat);
   if(kstat < 0) goto err2;
   if(kstat > 1) goto err1; /* Not connected. */

   /* Set direction from pt1 to pt2. */

   pt1->curve_dir[index1] |= 1;
/*   pt2->curve_dir[index2]  = (-1 ^ 33); */
   pt2->curve_dir[index2] = -31;
   pt2->curve_dir[index2] |= pt1->curve_dir[index1];


   goto out;

   /* Points are not connected. */
err1:

   *jstat = -1;
   s6err("sh6setdir",*jstat,0);
   goto out;

   /* Error in subfuction. */
err2:

   *jstat = -2;
   s6err("sh6setdir",*jstat,0);
   goto out;

   out :
      return;
}

