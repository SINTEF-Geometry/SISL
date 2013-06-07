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
