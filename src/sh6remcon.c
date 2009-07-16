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
 * $Id: sh6remcon.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6REMCON

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6remcon(SISLIntdat **pintdat,SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6remcon(pintdat,pt1,pt2,jstat)
   SISLIntdat **pintdat;
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove connection between pt1 and pt2.
*
* INPUT      : pintdat  - Pointer to pointer to the SISLIntdat data.
*              pt1      - First point.
*              pt2      - Second point.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat = -2  => error in space allocation
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
  int kstat;                 /* Local status variable.                     */
  
   *jstat = 0;
  

  sh6disconnect(pt1,pt2,&kstat);
  if (kstat < 0) goto error;

  /* Do we want to remove either point from pintdat? */

  
  goto out;  
  

  /* Error in lower level routine. */

  error : *jstat = kstat;
  s6err("sh6remcon",*jstat,0);
  goto out;                       
  
  
  out: ;
}
