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
 * $Id: sh6nmbmain.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6NMBMAIN

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6nmbmain(SISLIntpt *pt,int *jstat)
#else
int sh6nmbmain(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt, return the number of main points
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
       if(pt->pnext[ki] == NULL) goto err1;
       if(sh6ismain(pt->pnext[ki])) num++;
   }

   goto out;
   

err1:
   /* Error in data structure. */
   
   *jstat = -1;
   s6err("sh6nmbmain",*jstat,0);
   goto out;
   
   
   out :
      return num;
}






