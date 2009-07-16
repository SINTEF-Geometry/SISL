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
 * $Id: sh6nmbhelp.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6NMBHELP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6nmbhelp(SISLIntpt *pt,int *jstat)
#else
int sh6nmbhelp(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt, return the number of help points
*              it is linked to.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*********************************************************************
*/
{
   int num; /* Number of lists. */
   int ki; /* Loop variable.  */

   num=0;

   /* Count number of main lists pt lies in. */

   for(ki=0; ki<pt->no_of_curves; ki++)
   {
       if(pt->pnext[ki] == SISL_NULL) goto err1;
       if(sh6ishelp(pt->pnext[ki])) num++;
   }

   goto out;
   

err1:
   /* Error in data structure. */
   
   *jstat = -1;
   s6err("sh6nmbhelp",*jstat,0);
   goto out;
   
   
   out :
      return num;
}






