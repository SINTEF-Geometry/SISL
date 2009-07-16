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
 * $Id: s1775.c,v 1.2 2001-03-19 15:58:53 afr Exp $
 *
 */
#define S1775

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
s1775(SISLSurf *surf, double point[], int dim, double epsge,
	   double start[], double end[], double guess[], double clpar[],
	   int *stat)
#else
void s1775(surf, point, dim, epsge, start, end, guess, clpar, stat)
     SISLSurf   *surf;
     double point[];
     int dim;
     double epsge;
     double start[];
     double end[];
     double guess[];
     double clpar[];
     int    *stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a surface and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameters is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : surf    - The surface in the closest point problem.
*              point   - The point in the closest point problem.
*              dim     - Dimension of the geometry.
*              epsge   - Geometry resolution.
*              start   - Surface parameters giving the start of the search
*                        area (umin, vmin).
*              end     - Surface parameters giving the end of the search
*                        area (umax, vmax).
*              guess   - Surface guess parameters for the closest point
*                        iteration.
*
*
*
* OUTPUT     : clpar   - Resulting surface parameters from the iteration.
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SINTEF Oslo, August 1995.
*
*********************************************************************
*/                       
{ 
   int kpos = 0;             /* Error indicator. */
   SISLPoint* ppoint = SISL_NULL; /* SISL point.      */

   /* Generate a SISL point. */
   
   ppoint = newPoint(point, dim, 0);
   
   /* Call s1773. */
   
   s1773(ppoint, surf, epsge, start, end, guess, clpar, stat);
   if (*stat < 0) goto error;
  
   goto out;
  
   /* Error in lower level routine.  */

   error :
   s6err("s1775", *stat, kpos);
   goto out;

 out:    
    if (ppoint != SISL_NULL) freePoint(ppoint);
}

