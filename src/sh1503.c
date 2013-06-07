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
 * $Id: sh1503.c,v 1.3 2001-03-19 15:59:04 afr Exp $
 *
 */


#define SH1503

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh1503(SISLSurf *ps1,double base[],double norm[],double axisA[],double alpha,
       double ratio,int idim,double aepsco,double aepsge,
       int trackflag,int *jtrack,SISLTrack ***wtrack,
       int *jpt,double **gpar,int **pretop,int *jcrv,SISLIntcurve ***wcurve,
       int *jsurf,SISLIntsurf ***wsurf,int *jstat)
#else
void sh1503(ps1,base,norm,axisA,alpha,ratio,idim,aepsco,aepsge,trackflag,
	    jtrack,wtrack,jpt,gpar,pretop,jcrv,wcurve,jsurf,wsurf,jstat)
     SISLSurf     *ps1;
     double   base[];
     double   norm[];
     double   axisA[];
     double   alpha;
     double   ratio;
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
     int      *jsurf;
     SISLIntsurf ***wsurf;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find all intersections between a tensor-product surface
*              and a cone.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              base   - Base point of cone
*              norm   - Direction of cone axis
*              axisA  - One of the two ellipse axis vectors
*              alpha  - The opening angle of the cone at axisA
*              ratio  - The ratio of axisA to axisB
*              idim     - Dimension of the space in which the cone lies.
*              aepsco   - Computational resolution.
*              aepsge   - Geometry resolution.
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
* METHOD     : The vertices of the surface are put into the equation of
*              the cone resulting in a surface in one-dimensional space.
*              Then the zeroes of this surface are found.
*
*
* REFERENCES :
*
*-
* CALLS      : sh1761       - Perform point object-intersection.
*              s1320       - Put equation of surface into equation of implicit
*                            surface.
*              s1500       - Represent cone as implicit function.
*              hp_s1880       - Put intersections on output format.
*              make_sf_kreg   - Ensure k-regularity of surface.
*              newObject   - Create new object.
*              newPoint    - Create new point.
*              freeObject  - Free space occupied by an object.
*              freeIntdat  - Free space occupied by an intersection data.
*
* WRITTEN BY : Mike Floater, SI, 90-10.
* Changed by : Per OEyvind Hvidsten, SINTEF, 1194.
*              jcrv and jsurf set to zero when no intersection are found.
*
*********************************************************************
*/
{
  int kstat = 0;           /* Local status varible.                        */
  int kpos = 0;            /* Position of error.                           */
  int kdim = 1;            /* Dimension of space in which the point in the
			      intersect point/surface problem lies.        */
  double *spar = SISL_NULL;     /* Dummy array containing parameter values of
			      second object of single intersection points. */
  double spoint[1];        /* SISLPoint to intersect with object.              */
  double *scone = SISL_NULL;    /* Description of a cone as implicit surface. */
  SISLSurf *qs = SISL_NULL;         /* Pointer to surface in
			      surface/point intersection.*/
  SISLPoint *qp = SISL_NULL;        /* Pointer to point in
			      surface/point intersection.  */
  SISLObject *qo1 = SISL_NULL;      /* Pointer to surface in
			      object/point intersection. */
  SISLObject *qo2 = SISL_NULL;      /* Pointer to point in
			      object/point intersection    */
  SISLIntdat *qintdat = SISL_NULL;  /* Intersection result */
  int kdeg=2;         /* The degree of the implicit equation   */
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

  if (idim != qkreg -> idim) goto err106;

  /*
   * Allocate space for matrix describing a cone.
   * --------------------------------------------
   */

  if ((scone = newarray((idim+1)*(idim+1),double)) == SISL_NULL) goto err101;

  /*
   * Make a matrix of dimension (idim+1)x(idim+1) describing a
   * cone as an implicit function.
   * ---------------------------------------------------------
   */

  s1500(base,norm,axisA,alpha,ratio,idim,1,scone,&kstat);
  if (kstat < 0) goto error;

  /*
   * Put the description of the input surface into the implicit
   * equation for the cone.
   * ----------------------------------------------------------
   */

  s1320(qkreg,scone,1,0,&qs,&kstat);
  if (kstat < 0) goto error;

  /*
   * Create new object and connect surface to object.
   * ------------------------------------------------
   */

  if(!(qo1 = newObject(SISLSURFACE))) goto err101;
  qo1 -> s1 = qs;
  qo1 -> o1 = qo1;

  /*
   * Create new object and connect point to object.
   * ----------------------------------------------
   */

  if (!(qo2 = newObject(SISLPOINT))) goto err101;
  spoint[0] = DZERO;
  if (!(qp = newPoint(spoint,kdim,1))) goto err101;
  qo2 -> p1 = qp;

  /*
   * Find intersections.
   * -------------------
   */

  sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error;

  /* Represent degenerated intersection curves as one point.  */

  sh6degen(track_obj,track_obj,&qintdat,aepsge,&kstat);
  if (kstat < 0) goto error;

  /* Create tracks */
  if (trackflag && qintdat)
    {

      refine_all (&qintdat, track_obj, track_obj, scone, kdeg, aepsge, &kstat);
      if (kstat < 0)
	goto error;
    }


  /* Join periodic curves */
  int_join_per( &qintdat,track_obj,track_obj,scone,kdeg,aepsge,&kstat);
  if (kstat < 0)
    goto error;

  if (trackflag && qintdat)
    {
       make_tracks (track_obj, track_obj, kdeg, scone,
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
	       2,0,qintdat,jpt,gpar,&spar,pretop,jcrv,wcurve,jsurf,wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  else
  {
    *jcrv = 0;
    *jsurf = 0;
  }

  /*
   * Intersections found.
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
        s6err("sh1503",*jstat,kpos);
        goto out;

  /* Dimensions conflicting.  */

 err106: *jstat = -106;
        s6err("sh1503",*jstat,kpos);
        goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
        s6err("sh1503",*jstat,kpos);
        goto out;

 out:

  /* Free allocated space.  */

  if (spar)    freearray(spar);
  if (scone)   freearray(scone);
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
