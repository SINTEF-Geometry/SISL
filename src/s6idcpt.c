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
 * $Id: s6idcpt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6IDCPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idcpt(SISLIntdat *pintdat,SISLIntpt *pintpt,SISLIntpt **rintpt)
#else
void s6idcpt(pintdat,pintpt,rintpt)
     SISLIntdat *pintdat;
     SISLIntpt  *pintpt;
     SISLIntpt  **rintpt;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To find the point which is closest to pintpt
*              in the parametric space. If pintpt is the only
*              point in pintdat *rintpt is NULL.
*
*
*
* INPUT       :pintpt   - Pointer to an intersection point.
*              pintdat  - Pointer to intersection data.
*
*
* OUTPUT     : rintpt   - Pointer to a pointer to a point closest
*                         to pintpt.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  if (pintdat == NULL)
    *rintpt = NULL;
  else
    {
      int ki,knr;                /* Counters.          */
      double tdist,td;           /* To store distanse. */
      
      if (pintpt == pintdat->vpoint[0])
        tdist = HUGE;
      else
        tdist = s6dist(pintdat->vpoint[0]->epar,pintpt->epar,pintpt->ipar);
      
      for (knr=0,ki=1; ki<pintdat->ipoint; ki++)
        {
	  if (pintpt == pintdat->vpoint[ki])
	    td = HUGE;
	  else
	    td = s6dist(pintdat->vpoint[ki]->epar,pintpt->epar,pintpt->ipar);
	  
	  if (td < tdist)
	    {
	      knr = ki;
	      tdist = td;
	    }
        }
      
      if (tdist == HUGE)
        *rintpt = NULL;
      else
        *rintpt = pintdat->vpoint[knr];
    }
}

