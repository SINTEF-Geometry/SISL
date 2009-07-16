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
 * $Id: sh1510.c,v 1.3 2001-03-19 15:59:04 afr Exp $
 *
 */


#define SH1510

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1510(SISLSurf * ps,double eyepoint[],int idim,double aepsco,
	    double aepsge,int trackflag,int *jtrack,SISLTrack ***wtrack,
	    int *jpt,double **gpar,int **pretop,int *jcrv,
	    SISLIntcurve ***wcurve,int *jsurf,SISLIntsurf ***wsurf,int *jstat)
#else
void sh1510(ps, eyepoint, idim, aepsco, aepsge,	trackflag,jtrack,wtrack,
	    jpt,gpar,pretop,jcrv,wcurve,jsurf,wsurf,jstat)
     SISLSurf *ps;
     double eyepoint[];
     int idim;
     double aepsco;
     double aepsge;
     int       trackflag;
     int       *jtrack;
     SISLTrack ***wtrack;
     int *jpt;
     double **gpar;
     int      **pretop;
     int *jcrv;
     SISLIntcurve ***wcurve;
     int      *jsurf;
     SISLIntsurf ***wsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find the silhouette curves and points of a surface when
*              the surface is viewed perspectively from a specific eye point. In
*              addition to the points and curves found by this routine,
*              break curves and edge-curves might be silhouette curves.
*
*
*
* INPUT      : ps  -       Pointer to the surface.
*              eyepoint -  The eye point vector.
*              idim   -    Dimension of the space in which eyepoint lies.
*              aepsco -    Computational resolution.
*              aepsge -    Geometry resolution.
*              trackflag - If true, create tracks.
*
*
*
* OUTPUT     : jtrack - Number of tracks created
*              wtrack - Array of pointers to tracks
*              jpt    - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single silhouette points in the parameter
*                       plane of the surface. The points lie continuous.
*                       Silhouette curves are stored in wcurve.
*              pretop - Topology info. for single intersection points.
*              jcrv   - Number of silhouette curves.
*              wcurve - Array containing descriptions of the silhouette
*                       curves. The curves are only described by points
*                       in the parameter plane. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : sh1761       - Perform point object-intersection.
*              s1512       - Make a function whose zeroes describe the
*                            silhouette lines of a surface.
*              hp_s1880       - Put intersections on output format.
*              make_sf_kreg   - Ensure k-regularity of surface.
*              newObject   - Create new object.
*              newPoint    - Create new point.
*              freeObject  - Free space occupied by an object.
*              freeIntdat  - Free space occupied by an intersection data.
*
* WRITTEN BY : Mike Floater, SI, 90-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994. Added
*              initialization of 'jsurf'.
*
*********************************************************************
*/
{
  int kstat = 0;		/* Local status varible.                        */
  int kpos = 0;			/* Position of error.                           */
  int kdim = 1;			/* Dimension of space in which the point in the
				   intersect point and surface problem lies.    */
  double *spar = SISL_NULL;		/* Dummy array containing parameter values of
				   second object of single intersection points. */
  double spoint[1];		/* SISLPoint to intersect with object.              */
  SISLSurf *qs = SISL_NULL;		/* Surface whose zeroes are the silhouette
			       curves/points of the original surface.     */
  SISLPoint *qp = SISL_NULL;		/* Pointer to point in
			       surface/point intersection.  */
  SISLObject *qo1 = SISL_NULL;	/* Pointer to surface in
			       object/point intersection. */
  SISLObject *qo2 = SISL_NULL;	/* Pointer to point in
			       object/point intersection    */
  SISLIntdat *qintdat = SISL_NULL;	/* Intersection result */
  int kdeg = 1004;              /* The degree of the implicit equation */
  double simpli[16];             /* Array containing the implicit description */

  SISLObject *track_obj=SISL_NULL;
  SISLSurf *qkreg=SISL_NULL; /* Input surface ensured k-regularity. */

  /* -------------------------------------------------------- */

  if (ps->cuopen_1 == SISL_SURF_PERIODIC ||
      ps->cuopen_2 == SISL_SURF_PERIODIC)
  {
     /* Cyclic surface. */

     make_sf_kreg(ps,&qkreg,&kstat);
     if (kstat < 0) goto error;
   }
  else
    qkreg = ps;

  /*
  * Create new object and connect surface to object.
  * ------------------------------------------------
  */

  if (!(track_obj = newObject (SISLSURFACE)))
    goto err101;
  track_obj->s1 = ps;

  simpli[0] = eyepoint[0];
  simpli[1] = eyepoint[1];
  simpli[2] = eyepoint[2];






  /*
   * Check dimension.
   * ----------------
   */

  *jpt = 0;
  *jcrv = 0;
  *jtrack = 0;
  *jsurf = 0;

  if (idim != qkreg->idim)
    goto err106;

  /*
   * Make surface whose zeroes describe silhouette lines.
   * ----------------------------------------------------
   */

  s1512 (qkreg, eyepoint, idim, &qs, &kstat);
  if (kstat < 0)
    goto error;

  /*
   * Create new object and connect surface to object.
   * ------------------------------------------------
   */

  if (!(qo1 = newObject (SISLSURFACE)))
    goto err101;
  qo1->s1 = qs;
  qo1->o1 = qo1;

  /*
   * Create new object and connect point to object.
   * ----------------------------------------------
   */

  if (!(qo2 = newObject (SISLPOINT)))
    goto err101;
  spoint[0] = DZERO;
  if (!(qp = newPoint (spoint, kdim, 1)))
    goto err101;
  qo2->p1 = qp;

  /*
   * Find intersections.
   * -------------------
   */

  sh1761 (qo1, qo2, aepsge, &qintdat, &kstat);
  if (kstat < 0)
    goto error;

  /* Represent degenerated intersection curves as one point.  */

  sh6degen(track_obj,track_obj,&qintdat,aepsge,&kstat);
  if (kstat < 0) goto error;

  /* Split curves at knots where the surface is only C1. */

  spli_silh(&qintdat,track_obj,&kstat);
  if (kstat < 0) goto error;

  /* Create tracks */
  if (trackflag && qintdat)
    {

      refine_all (&qintdat, track_obj, track_obj, simpli, kdeg, aepsge, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* Join periodic curves */
  int_join_per( &qintdat,track_obj,track_obj,simpli,kdeg,aepsge,&kstat);
  if (kstat < 0)
    goto error;

  if (trackflag && qintdat)
    {
      make_tracks (track_obj, track_obj, kdeg, simpli,
		   qintdat->ilist, qintdat->vlist,
		   jtrack, wtrack, aepsge, &kstat);
      if (kstat < 0)
	goto error;

    }

  /*
   * Express silhouette curves and points on output format.
   * -----------------------------------------------------
   */

  if (qintdat)			/* Only if there were intersections found */
    {
      hp_s1880(track_obj, track_obj, kdeg,
	       2,0,qintdat,jpt,gpar,&spar,pretop,jcrv,wcurve,jsurf,wsurf,&kstat);
      if (kstat < 0)
	goto error;
    }


  /*
   * Silhouette curves and points found.
   * ----------------------------------
   */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh1510", *jstat, kpos);
  goto out;

  /* Dimensions conflicting.  */

err106:*jstat = -106;
  s6err ("sh1510", *jstat, kpos);
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  s6err ("sh1510", *jstat, kpos);
  goto out;

out:

  /* Free allocated space.  */

  if (spar)
    freearray (spar);
  if (qo1)
    freeObject (qo1);
  if (qo2)
    freeObject (qo2);
  if (qintdat)
    freeIntdat (qintdat);
  if (track_obj)
    {
       track_obj->s1 = SISL_NULL;
       freeObject(track_obj);
    }

  /* Free local surface.  */
    if (qkreg != SISL_NULL && qkreg != ps) freeSurf(qkreg);

  return;
}
