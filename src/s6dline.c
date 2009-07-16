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
 * $Id: s6dline.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DLINE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
  s6dline(double estart[],double eend[],double epoint[],
	  int idim,int *jstat)
#else
double s6dline(estart,eend,epoint,idim,jstat)
   double estart[];
   double eend[];
   double epoint[];
   int    idim;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between a line segment and a point.
*
*
*
* INPUT      : estart - Start point of line segment. Dimension is idim.
*              eend   - End point of line segment. Dimension is idim.
*              epoint - Point. Dimension is idim.
*              idim   - Dimension of geometry space.
*
*
*
* OUTPUT     : s6dline - Distance between line and point.
*              jstat   - status messages 
*                        = 2      : Zero line segment. Closest point
*                                   is outside segment.
*                        = 1      : Closest point on the extended
*                                   line through estart and end is
*                                    outside the segment.
*                        = 0      : ok
*                        < 0      : error
*
*
* METHOD     : 
*              
*
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   -  Scalar product between two vectors.
*              s6diff   -  Difference vector between two vectors.  
*              s6length -  Length of vector.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
   int kstat = 0;         /* Local status varaible.           */
   int ki;                /* Counter.                         */
   double tpar;           /* Parameter of closest point.      */
   double tdist;          /* Distance between point and line. */
   double t1;             /* Scalar product.                  */
   double *sline = SISL_NULL;  /* Line vector.                     */
   double *sdiff = SISL_NULL;  /* Difference vector.               */
   
   /* Allocate scratch for local vectors.  */
   
   if ((sline = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   if ((sdiff = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Compute help vectors.  */
   
   s6diff(eend,estart,idim,sline);
   s6diff(epoint,estart,idim,sdiff);
   
   /* Compute parameter of closest point. */
   
   t1 = s6scpr(sline,sline,idim);
   if (t1 <= REL_COMP_RES) 
   {
      /* Compute distance between point and first endpoint of line. */
      
      tdist = s6dist(estart,epoint,idim);
       
      /* Set a warning.  */
      
      *jstat = 2;
      goto out;
   }
   
   tpar = s6scpr(sline,sdiff,idim)/t1;
   
   /* Compute vector between input point and closest point on
      line.      */
   
   for (ki=0; ki<idim; ki++)
      sdiff[ki] = estart[ki] + tpar*sline[ki] - epoint[ki];
   
   /* Compute length of vector.  */
   
   tdist = s6length(sdiff,idim,&kstat);
   
   /* Set status.  */
   
   *jstat = (tpar < 0 || tpar > 1) ? 1 : 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   out :
      /* Free space occupied by local arrays.  */
      
      if (sline != SISL_NULL) freearray(sline); 
      if (sdiff != SISL_NULL) freearray(sdiff);
			 
      return tdist;
 }
