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
 * $Id: s1921.c,v 1.3 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1921

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1921(SISLSurf *ps1,double edir[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1921(ps1,edir,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   edir[];
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
* PURPOSE    : Find the extremal points/curves of the surface ps1 in
*              the direction edir.
*
*
* INPUT      : ps1    - Pointer to the surface.
*              edir   - The direction in which the extremal point(s)
*                       and/or interval(s) are to be calculated. If
*                       idim=1 a positive value indicates the maximum
*                       of the B-spline function and a negative value
*                       the minimum. If the dimension is greater that
*                       1 the array contains the coordinates of the
*                       direction vector.
*              idim   - Dimension of the space in which the vector edir
*                       lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single extremal points.
*              gpar   - Array containing the parameter values of the
*                       single extremal points in the parameter
*                       area of the surface. The points lie continuous.
*                       Extremal curves are stored in wcurve.
*              jcrv   - Number of extremal curves.
*              wcurve - Array containing descriptions of the extremal
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
* METHOD     : The scalar-product of the coefficients of the surface with
*              the vector edir is calculated giving a surface with dimension
*              1. Then the maxima of this surface is calculated.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1880 - Put extremal points/intervals into output format.
*              s1161 - Find maxima of one-dimensional object.
*              s6scpr - Scalar-product between two vectors.
*              newSurf    - Create new surface
*              newObject  - Create new object.
*              newIntdat  - Create new max point .
*              freeObject - Free the space occupied by a given object.
*              freeIntdat  - Free space occupied by an extremal data.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-10.
* REVISED BY : Mike Floater, SI, 91-09.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input surface with periodic basis correctly.
*
*********************************************************************
*/
{
  int ikind;               /* Type of surface ps1 is.                     */
  int kstat = 0;           /* Local status variable.                      */
  int kpos = 0;            /* Position of error.                          */
  int ki;                  /* Counter.                                    */
  int kn1,kn2;             /* Number of vertices of surface.              */
  int kk1,kk2;             /* Order of surface.                           */
  double tmax;             /* Estimate of maximal value of 1-dim. surface.*/
  double *st1,*st2;        /* Pointer to knotvectors of surface.          */
  double *scoef;           /* Pointer to vertices of surface.             */
  double *sc = SISL_NULL;       /* Pointer to vertices of surface in maxima
			      calculation.                                */
  double *spar = SISL_NULL;     /* Values of maxima in the parameter area of
			      the second object. Empty in this case.      */
  double *s1,*s2,*sstop;   /* Pointers used to traverse double-arrays.    */
  SISLIntdat *qintdat = SISL_NULL;  /* Pointer to max data structure.*/
  SISLSurf *qs = SISL_NULL;         /* Pointer to 1-dim. surface in maxima-calculation.*/
  SISLObject *qo1 = SISL_NULL;      /* Pointer to object in maxima-calculation.  */
  SISLSurf *qkreg = SISL_NULL;      /* Input surface with k-regularity ensured.  */


  /* Ensure k-regular input surface. */

  if ( ps1 -> cuopen_1 == SISL_SURF_PERIODIC ||
       ps1 -> cuopen_2 == SISL_SURF_PERIODIC )
  {
    /* Cyclic (periodic) surface */

    make_sf_kreg(ps1, &qkreg, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg = ps1;


  /* Check dimension.  */

  if ( qkreg -> idim != idim )  goto err106;

  /* Describe surface with local variables.  */

  kn1 = qkreg -> in1;
  kn2 = qkreg -> in2;
  kk1 = qkreg -> ik1;
  kk2 = qkreg -> ik2;
  st1 = qkreg -> et1;
  st2 = qkreg -> et2;
  ikind = qkreg -> ikind;

  if ( ikind == 2 || ikind == 4 )
  {
    scoef = qkreg -> rcoef;
    /* Allocate space for coeffecients of new surface.  */

    if ( (sc = newarray(2*kn1*kn2, DOUBLE)) == SISL_NULL )  goto err101;

    /* Compute scalar-product of surface-vertices and direction vector. */
    /* Copy over weights. */

    for ( s1=scoef, s2=sc, sstop=s2+2*kn1*kn2;  s2 < sstop;  s1+=idim+1, s2+=2 )
    {
      *s2 = s6scpr(s1, edir, idim);
      *(s2+1) = *(s1+idim);
    }
  }
  else
  {
    scoef = qkreg -> ecoef;
    /* Allocate space for coeffecients of new surface.  */

    if ( (sc = newarray(kn1*kn2, DOUBLE)) == SISL_NULL )  goto err101;

    /* Compute scalar-product of surface-vertices and direction vector. */

    for ( s1=scoef, s2=sc, sstop=s2+kn1*kn2;  s2 < sstop;  s1+=idim, s2++ )
      *s2 = s6scpr(s1, edir, idim);
  }


  /* Create new surface.  */

  qs = newSurf(kn1, kn2, kk1, kk2, st1, st2, sc, qkreg->ikind, 1, 1);
  if ( qs == SISL_NULL )  goto err101;

  /* Create new object and connect surface to object.  */

  qo1 = newObject(SISLSURFACE);
  if ( qo1 == SISL_NULL )  goto err101;
  qo1 -> s1 = qs;

  /* Find maxima.  */

  /* Find maxima. */
  tmax = -(double)HUGE;

  s1161(qo1, &tmax, aepsge, &qintdat, &kstat);
  if ( kstat < 0 )  goto error;

  if (qintdat)
  {

    /* Express maximal points/intervals on output format.  */
    s1880(2, 0, &qintdat->ipoint, qintdat->vpoint,
	  &qintdat->ilist, qintdat->vlist,
	  jpt, gpar, &spar, jcrv, wcurve, &kstat);
    if ( kstat < 0 )  goto error;

    /* Handle periodicity (remove extraneous points) */

    if ( *jpt > 1  &&  idim > 1 && (ps1 -> cuopen_1 == SISL_SURF_PERIODIC ||
				    ps1 -> cuopen_2 == SISL_SURF_PERIODIC ) )
    {
      for ( ki=0; ki < (*jpt); ki++ )
      {
	if ( (ps1 -> cuopen_1 == SISL_SURF_PERIODIC &&
	      (*gpar)[2*ki]   == ps1 -> et1[ps1->in1]) ||
	     (ps1 -> cuopen_2 == SISL_SURF_PERIODIC &&
	      (*gpar)[2*ki+1] == ps1 -> et2[ps1->in2]) )
	{
	  (*jpt)--;
	  (*gpar)[2*ki]   = (*gpar)[2*(*jpt)];
	  (*gpar)[2*ki+1] = (*gpar)[2*(*jpt)+1];
	  ki--;
	}
      }
    }
  }

  /* Extremal points/intervals found.  */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101:
  *jstat = -101;
  s6err("s1921",*jstat,kpos);
  goto out;

  /* Dimensions conflicting.  */

 err106:
  *jstat = -106;
  s6err("s1921",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1921",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space.  */

  if ( qkreg && qkreg != ps1 )  freeSurf(qkreg);
  if (sc) freearray(sc);
  if (spar) freearray(spar);
  if (qo1) freeObject(qo1);
  if (qintdat) freeIntdat(qintdat);

  return;
}
