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
 * $Id: sh1858.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */


#define SH1858

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
sh1858(SISLSurf *ps1,SISLCurve *pc1,double aepsco,double aepsge,
	int trackflag, int *jtrack, SISLTrack *** wtrack,
	   int *jpt,double **gpar1,double **gpar2,int **pretop,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void sh1858(ps1,pc1,aepsco,aepsge,trackflag,jtrack,wtrack,jpt,
	    gpar1,gpar2,pretop,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     SISLCurve    *pc1;    
     double   aepsco;
     double   aepsge;
     int       trackflag;
     int       *jtrack;
     SISLTrack ***wtrack;
     int      *jpt;
     double   **gpar1;
     double   **gpar2;
     int      **pretop;
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
*              trackflag - If true, create tracks.
*
*
*
* OUTPUT     : jtrack - Number of tracks created
*              wtrack - Array of pointers to tracks
*              jpt    - Number of single intersection points.
*              gpar1  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the first curve. The points lie 
*                       continuous. Intersection curves are stored in wcurve.
*              gpar2  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the second curve.
*              pretop - Topology info. for single intersection points.
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
* CALLS      : sh1761      - Perform object/object-intersection.
*              hp_s1880      - Put intersections on output format.
*              newObject  - Create new object.
*              freeObject - Free space occupied by an object.
*              freeIntdat - Free space occupied by the intdat structure.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-10.
* REWRITTEN BY : B.O. Hoset, SI, 89-06.
*
*********************************************************************
*/
{            
  double *nullp = SISL_NULL;
  int kstat = 0;                 /* Local status variable.                   */
  int kpos = 0;                  /* Position of error.                       */
  SISLObject *qo1 = SISL_NULL;            /* Object containing first curve in 
				    the intersection.                        */
  SISLObject *qo2 = SISL_NULL;            /* Object containing second curve in 
				    the intersection.                        */
  SISLIntdat *qintdat = SISL_NULL;        /* Structure holding the intersection data. */
  int      ksurf=0;         /* Dummy number of Intsurfs. */
  SISLIntsurf **wsurf=SISL_NULL;    /* Dummy array of Intsurfs. */
  int kdeg=0;
  
  /* 
   * Check dimensions.  
   * -----------------
   */

  *jpt  = 0;
  *jcrv = 0;
  *jtrack = 0;

  if (ps1 -> idim != pc1 -> idim) goto err106;

  /* 
   * Create objects and connect surface/curve to the objects.  
   * --------------------------------------------------------
   */

  if ((qo1 = newObject(SISLSURFACE)) == SISL_NULL) goto err101;
  qo1 -> s1 = ps1;
  qo1 -> o1 = qo1;
  
  if ((qo2 = newObject(SISLCURVE)) == SISL_NULL) goto err101;
  qo2 -> c1 = pc1;
  qo2 -> o1 = qo2;
  
  /* 
   * Find intersections.  
   * -------------------
   */

  sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error;

  /* Represent degenerated intersection curves as one point.  */

  sh6degen(qo1,qo2,&qintdat,aepsge,&kstat);
  if (kstat < 0) goto error;

  /* Join periodic curves */
  int_join_per( &qintdat,qo1,qo2,nullp,kdeg=0,aepsge,&kstat);
  if (kstat < 0)
    goto error;

  /* Create tracks */
  if (trackflag && qintdat)
    {
      make_tracks (qo1, qo2, 0, nullp,
		   qintdat->ilist, qintdat->vlist, 
		   jtrack, wtrack, aepsge, &kstat);
      if (kstat < 0)
	goto error;

    }

  /* 
   * Express intersections on output format.  
   * ---------------------------------------
   */

  if (qintdat)/* Only if there were intersections found */
    {
      hp_s1880(qo1, qo2, kdeg,
	       2,1,qintdat,jpt,gpar1,gpar2,pretop,jcrv,wcurve,&ksurf,&wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* 
   * Intersections found.  
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* 
   * Error in space allocation.  
   * --------------------------
   */

 err101: *jstat = -101;                
        s6err("sh1858",*jstat,kpos);
        goto out;

  /* Dimensions conflicting.  */

 err106: *jstat = -106;
  s6err("sh1858",*jstat,kpos);
        goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
        s6err("sh1858",*jstat,kpos);
        goto out;

 out:

  /* 
   * Free allocated space.  
   * ---------------------
   */

  if (qo1) 
    {
      qo1 -> s1 = SISL_NULL;  freeObject(qo1);
    }
  if (qo2) 
    {
      qo2 -> c1 = SISL_NULL;  freeObject(qo2);
    }
  if (qintdat) freeIntdat(qintdat);

  /*
   * Exit sh1858.
   * -----------
   */
                                        
  return;
}                                               

