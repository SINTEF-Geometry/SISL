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
 * $Id: sh6rempnt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6REMOVEPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6removept(SISLIntpt *pt1,SISLIntpt *pt2,SISLIntpt *ptold,int *jstat)
#else
void sh6removept(pt1,pt2,ptold,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   SISLIntpt *ptold;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove the pt ptold from the list which contains
*              pt1,ptold,pt2 consecutively.
*              Error if there is no such list.
*              The list is repaired afterwards.
*
* INPUT      : pt1      - First point.
*              pt2      - Second point.
*              ptold    - Point to be removed between pt1 and pt2.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1,ptold,pt2 are not connected
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
* CALLS      : s6err      - Gives error message.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
  int kstat;      /* Local status variables.          */
  
   *jstat = 0;
  


  sh6disconnect(pt1,ptold,&kstat);
  if(kstat < 0) goto error;

  sh6disconnect(pt2,ptold,&kstat);
  if(kstat < 0) goto error;

  sh6connect(pt1,pt2,&kstat);
  if(kstat < 0) goto error;

  
  goto out;
  


/* Error in subfunction. */

error:  *jstat = kstat;
        s6err("sh6removept",*jstat,0);
        goto out;


   out:
      return;
}
