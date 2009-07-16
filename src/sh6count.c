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
 * $Id: sh6count.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6COUNT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
int 
      sh6count(SISLIntpt *pt,int *jstat)
#else
int sh6count(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given a main Intpt 
*              return the number of points in the list (of
*              main points)
*              and determine whether it is open or closed
*              through *jstat. If the Intpt is not
*              in a unique list, say so.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              index    - Index specifying a list containing pt.
*              jstat    - Error flag.
*                         jstat = 0   => Successful. List is open.
*                         jstat = 1   => Successful. List is closed.
*                         jstat = 2   => pt is a junction pt.
*                         jstat = 3   => pt is isolated.
*                         jstat = -1  => Error in pt.
*                         jstat = -2  => List is inconsistent.
*
*
*
* METHOD     : 1. Traverse forwards from pt to one end of
*                 the list.
*              2. Then (unless list is closed) go to the
*                 other end of the list.
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. Sept. 91.
*
*********************************************************************
*/
{
   int       ki;      /* Counter.                        */
   SISLIntpt *pt1, *pt2;  /* Neighbours of pt.           */
   SISLIntpt *lastpt, *nowpt, *nextpt;
			 /* 3 adjacents pts in the list. */
   int       kstat;   /* Status variable.                */

   *jstat = 0;
   ki=1;

   if(pt == SISL_NULL) goto err1;
   if(!sh6ismain(pt)) goto err1;

   sh6getnhbrs(pt,&pt1,&pt2,&kstat);
   if(kstat < 0) goto err2; /* Bad list. */
   if(kstat == 2)
   {
       *jstat = 2;
       goto out;
   }
   if(kstat == 3)
   {
       *jstat = 3;
       goto out;
   }


   /* Traverse in the first direction unless at the end. */

   nextpt=pt1;
   nowpt=pt;
   while(nextpt != SISL_NULL && nextpt != pt)
   {
       ki++;
       lastpt=nowpt;
       nowpt=nextpt;
       sh6getother(nowpt,lastpt,&nextpt,&kstat);
       if(kstat < 0) goto err2; /* Bad list. */
   }

   /* Now if nextpt == pt the list is closed and we're finished. */

   if(nextpt == pt)
   {
       *jstat=1;
       goto out;
   }

   /* Otherwise traverse in the second direction. */

   nextpt=pt2;
   nowpt=pt;
   while(nextpt != SISL_NULL)
   {
       ki++;
       lastpt=nowpt;
       nowpt=nextpt;
       sh6getother(nowpt,lastpt,&nextpt,&kstat);
       if(kstat < 0) goto err2; /* Bad list. */
   }

   goto out;

err1:
   /* Error. Error in pt or index. */
   *jstat = -1;
   s6err("sh6count",*jstat,0);
   goto out;

err2:
   /* Error. List is inconsistent. */
   *jstat = -2;
   s6err("sh6count",*jstat,0);
   goto out;

   out :
      return ki;
}
