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
 * $Id: sh6tohelp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TOHELP

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6tohelp(SISLIntpt *pt,int *jstat)
#else
void sh6tohelp(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if pt is a mai point. If it is, transform to
*              help point, if not give a message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => pt was a main point, now help.
*                         jstat =  1  => pt is not a main point.
*                         jstat = -1  => Error in pt.
*                         jstat = -2  => Illegal to convert status
*                         jstat = -3  => Error in data structure.
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
   int kstat; /* Local status */
   int num; 

   *jstat=0;

   if(pt == SISL_NULL) goto err1;

   if(sh6ismain(pt))  /* If pt is a help point. */
   {
      /* ??????????? */
      /* if(pt->no_of_curves > 2) goto err2; */
      
      num=sh6nmbmain(pt,&kstat);
      /* Problem in sh6edgred when starting reduction */
      /* if(num > 1) goto err2; */

       pt->iinter = -pt->iinter;  /* Convert status to main point. */
   }
   else
   {
       *jstat=1;
   }

   goto out;
   

err1:
   /* Error in input. pt is null. */
   
   *jstat = -1;
   s6err("sh6tohelp",*jstat,0);
   goto out;
   
   /* Error, Illegal to change status. */
   
   /* err2:
    *jstat = -2;
   s6err("sh6tohelp",*jstat,0);
   goto out; */
   
   
   out :
      return;
}
