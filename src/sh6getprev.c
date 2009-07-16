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
 * $Id: sh6getprev.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETPREV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6getprev(SISLIntpt *pt1,SISLIntpt *pt2)
#else
int sh6getprev(pt1,pt2)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt pt1 and a pointer to another Intpt pt2,
*              fetch the index of the pt1 array corresponding
*              to pt2. If no such index exists return -1.
*
*
* INPUT      : pt1       - Pointer to the Intpt.
*              pt2     - Pointer to another Intpt.
*
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. May 91.
*
*********************************************************************
*/
{
   int       ncurv;   /* number of curves pt1 is connected to       */
   int       index;   /* index number for pnext array              */

   index = -1;

   if(pt1 == SISL_NULL || pt2 == SISL_NULL) goto out;

   ncurv = pt1->no_of_curves;  /* note ncurv can be zero */

   index=0;
   while(index < ncurv && pt1->pnext[index] != pt2) index++;
   if(index == ncurv) index = -1;  /* no index found */

   goto out;

   out :
      return index;
}
