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
 * $Id: s1954.c,v 1.3 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1954

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1954(SISLSurf *psurf,double epoint[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1954(psurf,epoint,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
	   SISLSurf     *psurf;
	   double   epoint[];
	   int      idim;
	   double   aepsco;
	   double   aepsge;
	   int      *jpt;
	   double   **gpar;
	   int      *jcrv;
	   SISLIntcurve ***wcurve;
	   int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find the closest point between the surface psurf and the
*              point epoint.
*
*
*
* INPUT      : psurf  - Pointer to the surface in the closest point problem.
*              epoint - The point in the closest point problem.
*              idim   - Dimension of the space in which epoint lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single closest points.
*              gpar   - Array containing the parameter values of the
*                       single closest points in the parameter area
*                       of the surface. The points lie continuous.
*                       Closest curves are stored in wcurve.
*              jcrv   - Number of closest curves.
*              wcurve - Array containing descriptions of the closest
*                       curves. The curves are only described by points
*                       in the parameter area. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Put the equation of the surface into the equation of a
*              sphere with center in the point epoint and radius zero,
*              achieving a one-dimensional surface. Find the minima of
*              this surface.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1321 - Create matrix describing the sphere as an
*                      implicit function.
*              s1320 - Put the description of the surface into the
*                      equation of the sphere.
*              s1921 - Find minimum of surface.
*              freeSurf - Free space occupied by a surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input surface with periodic basis correctly.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                    */
  int kpos = 0;             /* Position of error.                        */
  int kdim = 1;             /* Dimension of curve in extremal problem.   */
  double tradius = 0;       /* Radius of hyper-sphere describing point.  */
  double tdir = -1;         /* Direction of extremal value.              */
  double *sarray = SISL_NULL;    /* Matrix describing hyper-sphere.           */
  SISLSurf *qs = SISL_NULL;      /* Surface of which to find extremal points. */
  SISLSurf *qkreg = SISL_NULL;   /* Input surface with k-regularity ensured.  */
  int ratflag = 0;          /* Flag to indicate if surface is rational.  */
  int ki;                   /* Counter.                                  */


  *jstat = 0;

  /* Ensure k-regular basis */

  if ( psurf -> cuopen_1 == SISL_SURF_PERIODIC ||
       psurf -> cuopen_2 == SISL_SURF_PERIODIC )
  {
    /* Cyclic (i.e. periodic) surface */

    make_sf_kreg(psurf, &qkreg, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg = psurf;


  /* Test input.  */

  if ( qkreg -> idim != idim )  goto err106;

  if ( qkreg -> ikind == 2  ||  qkreg -> ikind == 4)  ratflag = 1;

  /* Allocate space for array describing a hyper-sphere.  */

  if ( (sarray = newarray((idim+1)*(idim+1), DOUBLE)) == SISL_NULL )  goto err101;

  /* Make a matrix of dimension (idim+1)*(idim+1) to describe
     the hyper-shpere.                                        */

  s1321(epoint, tradius, idim, kdim, sarray, &kstat);
  if ( kstat < 0 )  goto error;

  /* Put surface into equation of hyper-sphere.  */

  s1320(qkreg, sarray, kdim, ratflag, &qs, &kstat);
  if ( kstat < 0 )  goto error;

  /* Find minimum points of the new surface.  */

  s1921(qs, &tdir, kdim, aepsco, aepsge, jpt, gpar, jcrv, wcurve, &kstat);
  if ( kstat < 0 )  goto error;

  /* Handle periodicity (remove extraneous points) */
  if ( (*jpt) > 1  &&  idim > 1  && (psurf -> cuopen_1 == SISL_SURF_PERIODIC ||
				     psurf -> cuopen_2 == SISL_SURF_PERIODIC) )
  {
    for ( ki=0;  ki < (*jpt);  ki++ )
    {
      if ( (psurf -> cuopen_1 == SISL_SURF_PERIODIC &&
	    (*gpar)[2*ki]     == psurf -> et1[psurf->in1]) ||
	   (psurf -> cuopen_2 == SISL_SURF_PERIODIC &&
	    (*gpar)[2*ki+1]     == psurf -> et2[psurf->in2]) )
      {
	(*jpt)--;
	(*gpar)[2*ki]   = (*gpar)[2*(*jpt)];
	(*gpar)[2*ki+1] = (*gpar)[2*(*jpt)+1];
	ki--;
      }
    }
  }

  /* Closest points/intervals found.  */

  *jstat = 0;
  goto out;


  /* Error in space allocation.  */

 err101:
  *jstat = -101;
  s6err("s1954",*jstat,kpos);
  goto out;

  /* Dimensions conflicting.  */

 err106:
  *jstat = -106;
  s6err("s1954",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1954",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space.  */

  if ( qkreg  &&  qkreg != psurf )  freeSurf(qkreg);
  if (sarray)  freearray(sarray);
  if (qs)  freeSurf(qs);

  return;
}
