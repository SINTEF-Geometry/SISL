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
 * $Id: s1858.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1858

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1858(SISLSurf *ps1,SISLCurve *pc1,double aepsco,double aepsge,
	   int *jpt,double **gpar1,double **gpar2,int *jcrv,
	   SISLIntcurve ***wcurve,int *jstat)
#else
void s1858(ps1,pc1,aepsco,aepsge,jpt,gpar1,gpar2,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     SISLCurve    *pc1;    
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar1;
     double   **gpar2;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find all intersections between a B-spline surface
*              and a B-spline curve.
*
*
*
* INPUT      : ps1    - Pointer to the surface.
*              pc1    - Pointer to the curve.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single intersection points.
*              gpar1  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the first curve. The points lie 
*                       continuous. Intersection curves are stored in wcurve.
*              gpar2  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the second curve.
*              jcrv   - Number of intersection curves.
*              wcurve - Array containing descriptions of the intersection
*                       curves. The curves are only described by points
*                       in the parameter plane. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*                       If the curves given as input are degnenerate an
*                       intersection point can be returned as an intersection
*                       curve. Use s1327 to decide if an intersection curve
*                       is a point on one of the curves.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1858, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
*
*********************************************************************
*/
{  
  int kstat = 0;           /* Local status variable.                       */
  int kpos = 0;            /* Position of error.                           */
  int trackflag = 0;
  int jtrack;
  int *pretop=SISL_NULL;
  SISLTrack **wtrack=SISL_NULL;
          
  sh1858(ps1,pc1,aepsco,aepsge,trackflag,&jtrack,&wtrack,jpt,
	 gpar1,gpar2,&pretop,jcrv,wcurve,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);
     
  /* 
   * Intersections found.  
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error : 
    *jstat = kstat;
    s6err("s1858",*jstat,kpos);
    goto out;

  out:
    return;
}                                               
