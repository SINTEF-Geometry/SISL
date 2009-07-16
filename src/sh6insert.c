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
 * $Id: sh6insert.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6INSERT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6insert(SISLIntdat **pintdat,SISLIntpt *pt1,SISLIntpt *pt2,SISLIntpt **ptnew,int *jstat)
#else
void sh6insert(pintdat,pt1,pt2,ptnew,jstat)
   SISLIntdat **pintdat;
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   SISLIntpt **ptnew;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Place a new intersection point between two (connected)
*              points. 
*              NOTE that if there exists a point in pintdat which is
*              close to ptnew, then the input point is killed and
*              no insertion is done.
*              pt1 and pt2 must lie contiguously in a list.
*
* INPUT      : pintdat  - Pointer to pointer to the SISLIntdat data.
*              pt1      - First point.
*              pt2      - Second point.
*              ptnew    - Point to be placed between pt1 and pt2.
*              jstat    - Error flag.
*                         jstat =  1  => Warning, point existing in
*                                        pintdat, nothing done.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat = -2  => error in space allocation
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
* CALLS      : s6err      - Gives error message.
*              sh6idnpt    - Insert a new intpt structure.
*              copyIntpt  - Copy an intpt structure.
*              newIntdat  - Create new intdat structure.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* MODYFIED BY: UJK, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
  int kstat;                /* Local status variable.                     */
  
   *jstat = 0;
  
  /* First we have to be sure that pintdat contains ptnew. */
  
  sh6idnpt(pintdat,ptnew,1,&kstat);
  if (kstat < 0) goto error;
  if (kstat > 0) 
    { 
       /* Point already existing in data structure, point killed, 
	  no insertion */
       *jstat = 1;
       goto out;
    }
  
  /* UJK, aug. 92 insert always mainpts if one of the neighbour is a main */
  /* if (sh6ismain(pt1) && sh6ismain(pt2)) */
   if (sh6ismain(pt1) || sh6ismain(pt2))
     sh6tomain(*ptnew,&kstat);
  else
     sh6tohelp(*ptnew,&kstat);
  if (kstat < 0) goto error;

  /* Then insert the point. */
  sh6insertpt(pt1,pt2,*ptnew,&kstat);
  if (kstat < 0) goto error;
  
  
  goto out;
  

/* Error. pt1 and pt2 are not properly connected.  */


/* Error in sub function.  */

error:  *jstat = kstat;
        s6err("sh6insert",*jstat,0);
        goto out;

   out:
      return;
}
