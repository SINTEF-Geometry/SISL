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
 * $Id: s1859.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1859

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1859(SISLSurf *ps1,SISLSurf *ps2,double aepsco,double aepsge,
	   int *jpt,double **gpar1,double **gpar2,int *jcrv,
	   SISLIntcurve ***wcurve,int *jstat)
#else
void s1859(ps1,ps2,aepsco,aepsge,jpt,gpar1,gpar2,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     SISLSurf     *ps2;    
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
* PURPOSE    : Find all intersections between two B-spline surfaces.
*
*
*
* INPUT      : ps1    - Pointer to first surface.
*              ps2    - Pointer to second surface.
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
* METHOD     : The curves are subdivided until we have a simple problem
*              where the intersection can be found by iteration. Then
*              Newton iteration is used.
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1859, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
*
*********************************************************************
*/
{            
  int kstat = 0;               /* Local status variable.                   */
  int kpos = 0;                /* Position of error.                       */
  int i;

  int trackflag = 0;
  int jtrack;
  SISLTrack **wtrack=SISL_NULL;
  int *pretop=SISL_NULL;
  int jsurf;
  SISLIntsurf **wsurf=SISL_NULL;

  sh1859 (ps1, ps2, aepsco, aepsge, trackflag, &jtrack, &wtrack, jpt, gpar1, 
	  gpar2, &pretop, jcrv, wcurve,&jsurf,&wsurf,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  for(i=0; i<jsurf; i++)
    freeIntsurf(wsurf[i]);
  if(wsurf != SISL_NULL) freearray(wsurf);

  if(jsurf > 0) 
    *jstat=10;
  else 
    *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error : 
    *jstat = kstat;
    s6err("s1859",*jstat,kpos);
    goto out;

  out:

  /*
   * Exit s1859.
   * -----------
   */

    return;
}                                               
