/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: sh6setcnsd.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6SETCNSDIR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6setcnsdir(SISLIntpt *pt1,SISLIntpt *pt2,int ipar,int *jstat)
#else
void sh6setcnsdir(pt1,pt2,ipar,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       ipar;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set the direction of the curve through pt1 and pt2 such
*              that pt1 points to pt2 along a constant parameter line.
*              If they are not connected, give error.
*              NOTE: if the points are given singular status
*                    when they are of type SI_ORD.
*
* INPUT      : pt1      - Pointer to first Intpt.
*              pt2      - Pointer to second Intpt.
*              ipar     - parameter index: 0,1,2,3
*              jstat    - Error flag.
*                        jstat =  0  => Successful
*                        jstat = -1  => Points are not connected.
*                        jstat = -2  => Error in subfunction.
*                        jstat = -3  => Error in subfunction.
*
*
*
* REFERENCES : Same as sh6setdir, different values used in curve_dir
*              
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
   int kstat;         /* error flag. */
   int index1,index2; /* dummy indices.           */
   int sgn=1;
   
   *jstat = 0;
   /* Legal value on ipar ? */
   if (ipar < 0 || ipar > 3) goto err0;
			     
   /* Check if pt1 and pt2 are already connected. */

   sh6getlist(pt1,pt2,&index1,&index2,&kstat);
   if(kstat < 0) goto err2;
   if(kstat > 1) goto err1; /* Not connected. */
		 /*
		 if(pt1->iinter == SI_ORD)       pt1->iinter =  SI_SING;
		 else if(pt1->iinter == -SI_ORD) pt1->iinter = -SI_SING;
		 
		 if(pt2->iinter == SI_ORD)       pt2->iinter =  SI_SING;
		 else if(pt2->iinter == -SI_ORD) pt2->iinter = -SI_SING;
		 */
 /* Set constant direction between pt1 and pt2. */
   pt1->curve_dir[index1] |= (1<<(ipar+1));
   pt2->curve_dir[index2] |= (1<<(ipar+1));

   goto out;

   /* Wrong value on ipar. */
err0:

   *jstat = -3;
   s6err("sh6setcnsdir",*jstat,0);
   goto out;

   /* Points are not connected. */
err1:

   *jstat = -1;
   s6err("sh6setcnsdir",*jstat,0);
   goto out;

   /* Error in subfuction. */
err2:

   *jstat = -2;
   s6err("sh6setcnsdir",*jstat,0);
   goto out;

   out :
      return;
}

