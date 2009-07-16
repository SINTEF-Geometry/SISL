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
 * $Id: sh6tomain.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TOMAIN

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6tomain(SISLIntpt *pt,int *jstat)
#else
void sh6tomain(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if pt is a help point. If it is, transform to
*              main point, if not give a message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => pt was a help point, now main.
*                         jstat =  1  => pt is not a help point.
*                         jstat = -1  => Error in pt.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* MODYFIED BY: UJK, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
   int ki; /* Loop variable. */
   int num; 
   int kstat;
   
   *jstat=0;

   if(pt == SISL_NULL) goto err1;

   if(sh6ishelp(pt))  /* If pt is a help point. */
   {
       pt->iinter = -pt->iinter;  /* Convert status to main point. */

       /* Go through all neighbours and keep invariant:
	  not more than one mainpoint connected to a help point. */
       for(ki=0; ki<pt->no_of_curves; ki++) 
       {
	   if(sh6ishelp(pt->pnext[ki]))
	   {
	      /* UJK, change all NON-terminators to main */
	      /* num=sh6nmbmain(pt->pnext[ki],&kstat); */
	       num = pt->pnext[ki]->no_of_curves;
	       if(num > 1) sh6tomain(pt->pnext[ki],&kstat);
	   }
       }

   }
   else
   {
       *jstat=1;
   }

   goto out;
   

err1:
   /* Error in input. pt is null. */
   
   *jstat = -1;
   s6err("sh6tomain",*jstat,0);
   goto out;
   
   
   out :
      return;
}
