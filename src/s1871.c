/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1871

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1871(SISLCurve *pc1, double *pt1, int idim, double aepsge,
	 int *jpt,double **gpar1,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1871(pc1,pt1,idim,aepsge,jpt,gpar1,jcrv,wcurve,jstat)
     SISLCurve     *pc1;
     double    *pt1;
     int	idim;
     double   aepsge;
     int      *jpt;
     double   **gpar1;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find all intersections between a B-spline curve
*              and a point.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              pt1    - coordinates of the point.
*	       idim   - number of coordinates in pt1.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single intersection points.
*              gpar1  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the curve. The points lie
*                       continuous. Intersection curves are stored in wcurve.
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
* REFERENCES :
*
*-
* CALLS      : sh1871      - Perform the actual intersection.
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 9403.
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
  double aepsco = REL_COMP_RES;

  sh1871(pc1, pt1, idim, aepsco, aepsge, trackflag, &jtrack, &wtrack,
	 jpt, gpar1, &pretop, jcrv, wcurve, &kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  /*
   * Intersections found.
   * --------------------
   */

  *jstat = kstat;
  goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1871",*jstat,kpos);
    goto out;

  out:
    return;
}
