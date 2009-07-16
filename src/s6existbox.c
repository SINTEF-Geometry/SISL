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
 * $Id: s6existbox.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6EXISTBOX

#include "sislP.h" 

#if defined(SISLNEEDPROTOTYPES)
int s6existbox(SISLbox *pbox,int itype,double aepsge)
#else
int s6existbox(pbox,itype,aepsge)
     SISLbox *pbox;
     int    itype;
     double aepsge;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if a particular box exist within an existing
*              box instance.
*
*
*
* INPUT      : pbox   - Box to test.
*              itype  - Kind of box to test existance of.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*              aepsge - Geometry resolution.
*
* OUTPUT     : s6existbox -  Status.
*                            -1 : Kind of box exist, but is expanded
*                                 with another tolerance.
*                             0 : Kind of box do not exist.
*                             1 : Requested box exist.
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
   if (pbox->e2min[itype] == SISL_NULL) return(0);  /* No box is made. */
   
   if (itype != 0 && DNEQUAL(pbox->etol[itype],aepsge))
      return(-1);  /* Box exist, but with another size of the expansion. */
   
   return(1);
}
   
