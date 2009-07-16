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
 * $Id: s6newbox.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6NEWBOX

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s6newbox(SISLbox *pbox,int inum,int itype,
	      double aepsge,int *jstat)
#else
void s6newbox(pbox,inum,itype,aepsge,jstat)
     SISLbox *pbox;
     int    inum;
     int    itype;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Create a particular box exist within an existing
*              box instance.
*
*
*
* INPUT      : pbox   - Box to modify.
*              inum   - Number of elements in min- and max-arrays.
*              itype  - Kind of box to create.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*              aepsge - Geometry resolution.
*
* OUTPUT     :  jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
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
* WRITTEN BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/                                     
{
   int knum = (inum == 1) ? inum : 2*inum;  /* If the geometry space has
					       dimension larger than 1,
					       a double set of min- and
					       max-arrays is to be made. */

   if (itype < 0 || itype > 2) goto err126;
   
   /* Test no such box exist, create the necessary arrays.  */
   
   if (pbox->e2min[itype] == SISL_NULL)
   {
      if ((pbox->e2min[itype] = newarray(knum,DOUBLE)) == SISL_NULL) goto err101;
      if ((pbox->e2max[itype] = newarray(knum,DOUBLE)) == SISL_NULL) goto err101;
   }
  
   /* Set the tolerance. */
   
   if (itype != 0) pbox->etol[itype] = aepsge;
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in input.  Kind of box do not exist.  */
   
   err126 : *jstat = -126;
   goto out;
   
   out :
      return;
}
