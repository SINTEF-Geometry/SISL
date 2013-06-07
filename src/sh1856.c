/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: sh1856.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */


#define SH1856

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1856(SISLSurf *ps1,double epoint[],double edir[],int idim,
	    double aepsco,double aepsge,
	    int trackflag, int *jtrack, SISLTrack *** wtrack,
	    int *jpt,double **gpar,int **pretop,int *jcrv,
	    SISLIntcurve ***wcurve,int *jstat)
#else
void sh1856(ps1,epoint,edir,idim,aepsco,aepsge,
	    trackflag,jtrack,wtrack,jpt,gpar,pretop,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   epoint[];
     double   edir[];
     int      idim;
     double   aepsco;
     double   aepsge;
     int       trackflag;
     int       *jtrack;
     SISLTrack ***wtrack;
     int      *jpt;
     double   **gpar;
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
* PURPOSE    : Find all intersections between a tensor-product surface
*              and an infinite straight line.
*
*
*
* INPUT      : ps1    - Pointer to surface.
*              epoint - SISLPoint on the line.
*              edir   - Direction vector of the line.
*              idim   - Dimension of the space in which the line lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              trackflag - If true, create tracks.
*
*
*
* OUTPUT     : jtrack - Number of tracks created
*              wtrack - Array of pointers to tracks
*              jpt    - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       plane of the surface. The points lie continuous. 
*                       Intersection curves are stored in wcurve.
*              pretop - Topology info. for single intersection points.
*              *jcrv  - Number of intersection curves.
*              wcurve  - Array containing descriptions of the intersection
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
* METHOD     : The line is described as the intersection between two
*              planes. The vertices of the surface are put into the equation 
*              of this planes achieving a surface in the two-dimentional 
*              space. Then the zeroes of this surface is found.
*
*
* REFERENCES :
*
*-
* CALLS      : sh1761 - Perform point object-intersection.
*              s1328 - Equation of surface into equations of two planes.
*              s1329 - Equation of surface into equation of plane.
*              make_sf_kreg   - Ensure k-regularity of surface.
*              hp_s1880 - Put intersections on output format.
*              s6twonorm - Make two vectors of length one that are normal
*                          to a 3d input vector.
*              newPoint    - Create new point.
*              newObject - Create new object.
*              freeObject - Free space occupied by an object.
*              freeIntdat  - Free space occupied by an intersection data.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
* REWRITTEN BY : Bjoern Olav Hoset, SI, 89-06.
*
*********************************************************************
*/
{            
  double *nullp = SISL_NULL;
  int kstat = 0;           /* Local status varible.                        */
  int kpos = 0;            /* Position of error.                           */
  int kdim;                /* Dimension of space in which the point in the
			      intersect point and surface problem lies.    */
  double *spar = SISL_NULL;     /* Dummy array containing parameter values of
			      second object of single intersection points. */
  double spoint[2];        /* SISLPoint to intersect with object.              */
  double *snorm1 = SISL_NULL;   /* Normal to direction vector of line.          */
  double *snorm2 = SISL_NULL;   /* Normal to direction vector of line and snorm1.*/
  SISLSurf *qs = SISL_NULL;         /* Pointer to surface in 
			      surface/point intersection.*/
  SISLPoint *qp = SISL_NULL;        /* Pointer to point in 
			      surface/point intersection.  */
  SISLObject *qo1 = SISL_NULL;      /* Pointer to surface in 
			      object/point intersection. */
  SISLObject *qo2 = SISL_NULL;      /* Pointer to point in 
			      object/point intersection    */
  SISLIntdat *qintdat = SISL_NULL;  /* Intersection result */
  int      ksurf=0;         /* Dummy number of Intsurfs. */
  SISLIntsurf **wsurf=SISL_NULL;    /* Dummy array of Intsurfs. */
  int      kdeg=2000;       /* input to int_join_per. */
  SISLObject *track_obj=SISL_NULL;
  SISLSurf *qkreg=SISL_NULL; /* Input surface ensured k-regularity. */

  /* -------------------------------------------------------- */  

  if (ps1->cuopen_1 == SISL_SURF_PERIODIC ||
      ps1->cuopen_2 == SISL_SURF_PERIODIC)
  {
     /* Cyclic surface. */

     make_sf_kreg(ps1,&qkreg,&kstat);
     if (kstat < 0) goto error;
   }
  else
    qkreg = ps1;
  
  /*
  * Create new object and connect surface to object.
  * ------------------------------------------------
  */
  
  if (!(track_obj = newObject (SISLSURFACE)))
    goto err101;
  track_obj->s1 = ps1;

  /* 
   * Check dimension.  
   * ----------------
   */

  *jpt  = 0;
  *jcrv = 0;
  *jtrack = 0;

  if (idim != 2 && idim != 3) goto err105;
  if (idim != qkreg -> idim) goto err106;

  /* 
   * Allocate space for normal vectors.  
   * ----------------------------------
   */

  snorm1 = newarray(idim,double);
  snorm2 = newarray(idim,double);
  if (snorm1 == SISL_NULL || snorm2 == SISL_NULL) goto err101;

  if (idim == 3)
    {

      /* 
       * Find two planes that intersect in the given line.  
       * -------------------------------------------------
       */

      s6twonorm(edir,snorm1,snorm2,&kstat);
      if (kstat < 0) goto error;

      /* 
       * Put the surface into the plane equations.  
       * -----------------------------------------
       */

      s1328(qkreg,epoint,snorm1,snorm2,idim,&qs,&kstat);
      if (kstat < 0) goto error;

      /*
       * Create new object and connect point to object.
       * ----------------------------------------------
       */

      kdim      = 2;
      spoint[0] = spoint[1] = DZERO;
      if (!(qo2  = newObject(SISLPOINT))) goto err101;
      if (!(qp   = newPoint(spoint,kdim,1))) goto err101;
      qo2 -> p1 = qp;
    }
  else if (idim == 2)
    {

      /* 
       * Find normal vector of line.  
       * ---------------------------
       */

      snorm1[0] = edir[1];
      snorm1[1] = (-1)*edir[0];

      /* 
       * Put surface into line-equation.  
       * -------------------------------
       */

      s1329(qkreg,epoint,snorm1,idim,&qs,&kstat);
      if (kstat < 0) goto error;

      /*
       * Create new object and connect point to object.
       * ----------------------------------------------
       */

      kdim      = 1;
      spoint[0] = DZERO;
      if (!(qo2  = newObject(SISLPOINT))) goto err101;
      if (!(qp   = newPoint(spoint,kdim,1))) goto err101;
      qo2 -> p1 = qp;
    }

  /* 
   * Create new object and connect surface to object.  
   * ------------------------------------------------
   */

  if(!(qo1 = newObject(SISLSURFACE))) goto err101;
  qo1 -> s1 = qs;
  qo1 -> o1 = qo1;

  /* 
   * Find intersections.  
   * -------------------
   */

  sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error;

  /* Represent degenerated intersection curves as one point.  */

  sh6degen(track_obj,track_obj,&qintdat,aepsge,&kstat);
  if (kstat < 0) goto error;

  /* Join periodic curves */
  int_join_per( &qintdat,track_obj,track_obj,nullp,kdeg,aepsge,&kstat);
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
      hp_s1880(track_obj, track_obj, kdeg,
	       2,0,qintdat,jpt,gpar,&spar,pretop,jcrv,wcurve,&ksurf,&wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* 
   * Intersections found.  
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
        s6err("sh1856",*jstat,kpos);
        goto out;

  /* Error in input. Dimension different from two or three.  */

 err105: *jstat = -105;
        s6err("sh1856",*jstat,kpos);
        goto out;

  /* Dimensions conflicting.  */

 err106: *jstat = -106;
        s6err("sh1856",*jstat,kpos);
        goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
        s6err("sh1856",*jstat,kpos);
        goto out;

 out:

  /* Free allocated space.  */

  if (snorm1)  freearray(snorm1);
  if (snorm2)  freearray(snorm2);
  if (spar)    freearray(spar);
  if (qo1)     freeObject(qo1);
  if (qo2)     freeObject(qo2);
  if (qintdat) freeIntdat(qintdat);
  if (track_obj)
    {
       track_obj->s1 = SISL_NULL;
       freeObject(track_obj);
    }

  /* Free local surface.  */
    if (qkreg != SISL_NULL && qkreg != ps1) freeSurf(qkreg);

return;
}                                               

