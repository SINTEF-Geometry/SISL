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
 * $Id: sh6getnext.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETNEXT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt* 
      sh6getnext(SISLIntpt *pt,int index)
#else
SISLIntpt* sh6getnext(pt,index)
   SISLIntpt *pt;
   int       index;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt and an index, fetch the next point
*              given by index.
*              If error, return SISL_NULL.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              index    - Index of link at pt.
*
*
* OUTPUT     : 
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{

   SISLIntpt *nextpt = SISL_NULL;

   /* check if index is within range */

   if(pt != SISL_NULL &&
      index >= 0 &&
      index < pt->no_of_curves) nextpt = pt->pnext[index];

   goto out;

   
   out :
      return nextpt;
}
