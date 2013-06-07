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
 * $Id: s1955.c,v 1.3 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1955

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1955(SISLCurve *pc1,SISLCurve *pc2,double aepsco,double aepsge,
	   int *jpt,double **gpar1,double **gpar2,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1955(pc1,pc2,aepsco,aepsge,jpt,gpar1,gpar2,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     SISLCurve    *pc2;
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
* PURPOSE    : Find the closest point between the two curve pc1 and pc2.
*
*
*
* INPUT      : pc1    - Pointer to the first curve in the closest point
*                       problem.
*              pc2    - Pointer to the second curve in the closest point
*                       problem.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single closest points.
*              gpar1  - Array containing the parameter values of the
*                       single closest points in the parameter intervals
*                       of the first curve. The points lie continuous.
*                       Closest curves are stored in wcurve.
*              gpar2  - Array containing the parameter values of the
*                       single closest points in the parameter intervals
*                       of the second curve. The points lie continuous.
*              jcrv   - Number of closest curves.
*              wcurve - Array containing descriptions of the closest
*                       curves. The curves are only described by points
*                       in the parameter area. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
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
* METHOD     : Generate the difference surface between the two curves,
*              and find those points of this surface closest to origo.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1954    - Find closest points between a surface and a point.
*              s1956    - Create difference surface between two curves.
*              freeSurf - Free space occupied by a surface.
*              newIntcurve - Create new closest curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input curves with periodic basis correctly.
*
*********************************************************************
*/
{
  int kstat = 0;               /* Local status variable.                 */
  int kpos = 0;                /* Position of error.                     */
  int ki,kj;                   /* Counters.                              */
  int kdim;                    /* Dimension of the space in which the
				  curves lie.                            */
  double *sorigo = SISL_NULL;       /* Array representing origo.              */
  double *spar1 = SISL_NULL;        /* Parameter-values in the parameter
				  interval of the first curve of the
				  closest curve.                         */
  double *spar2 = SISL_NULL;        /* Parameter-values in the parameter
				  interval of the second curve of the
				  closest curve.                         */
  double *ssing = SISL_NULL;        /* Single closest points between
				  difference surface and origo.          */
  SISLSurf *qsdiff = SISL_NULL;         /* Difference surface.                    */
  SISLIntcurve *qcurve;            /* Pointer to closest curve.              */
  SISLCurve *qkreg1 = SISL_NULL;    /* First incrv with ensured k-regular basis. */
  SISLCurve *qkreg2 = SISL_NULL;    /* Second incrv with ensured k-regular basis. */
  int k1=1,k2=2,k4=4;          /* Constants                              */


  /* Ensure k-regular basis for input curves. */

  if ( pc1 -> cuopen == SISL_CRV_PERIODIC )
  {
    /* Cyclic (i.e. periodic) basis. */

    make_cv_kreg(pc1, &qkreg1, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg1 = pc1;

  if ( pc2 -> cuopen == SISL_CRV_PERIODIC )
  {
    /* Cyclic (i.e. periodic) basis. */

    make_cv_kreg(pc2, &qkreg2, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg2 = pc2;


  /* Test input.  */

  kdim = qkreg1 -> idim;
  if ( kdim != qkreg2 -> idim )  goto err106;

  /* Allocate space for the point origo.  */

  if ( (sorigo = new0array(kdim, DOUBLE)) == SISL_NULL )  goto err101;

  /* Generate the difference surface between the curves.  */

  s1956(qkreg1, qkreg2, &qsdiff, &kstat);
  if ( kstat < 0 )  goto error;

  if ( kstat > 0 )
  {

    /* The curves define a closest interval. Allocate space
       for points defining the interval.                      */

    if ( (spar1 = newarray(k2, DOUBLE)) == SISL_NULL )  goto err101;
    if ( (spar2 = newarray(k2, DOUBLE)) == SISL_NULL )  goto err101;

    if ( kstat == 1 )
    {

      /* The curves are equal.  */

      spar1[0] = *(qkreg1->et + qkreg1->ik - 1);
      spar1[1] = *(qkreg1->et + qkreg1->in);
      spar2[0] = *(qkreg2->et + qkreg2->ik - 1);
      spar2[1] = *(qkreg2->et + qkreg2->in);

    }
    else if ( kstat == 2 )
    {

      /* The curves are equal, but oppositely directed.  */

      spar1[0] = *(qkreg1->et + qkreg1->ik - 1);
      spar1[1] = *(qkreg1->et + qkreg1->in);
      spar2[0] = *(qkreg2->et + qkreg2->in);
      spar2[1] = *(qkreg2->et + qkreg2->ik - 1);
    }

    *wcurve = SISL_NULL;
    if ( (*wcurve = newarray(k1, SISLIntcurve*)) == SISL_NULL )  goto err101;
    *jcrv = 1;

    **wcurve = SISL_NULL;
    if ( (**wcurve = newIntcurve(k2,k1,k1,spar1,spar2,k4)) == SISL_NULL )  goto err101;

  }
  else
  {

    /* Find the closest points between the difference function and origo.*/

    s1954(qsdiff, sorigo, kdim, aepsco, aepsge, jpt, &ssing, jcrv, wcurve, &kstat);
    if ( kstat < 0 )  goto error;

    /* Handle periodicity (remove extraneous points) */

    if ( (*jpt) > 1  &&  kdim > 1  && (pc1 -> cuopen == SISL_CRV_PERIODIC ||
				       pc2 -> cuopen == SISL_CRV_PERIODIC) )
    {
      for ( ki=0;  ki < (*jpt);  ki++ )
      {
	if ( (pc1 -> cuopen == SISL_CRV_PERIODIC &&
	      ssing[2*ki]   == pc1 -> et[pc1->in]) ||
	     (pc2 -> cuopen == SISL_CRV_PERIODIC &&
	      ssing[2*ki+1] == pc2 -> et[pc2->in]) )
	{
	  (*jpt)--;
	  ssing[2*ki]   = ssing[2*(*jpt)];
	  ssing[2*ki+1] = ssing[2*(*jpt)+1];
	  ki--;
	}
      }
    }

    if ( (*jpt) > 0 )
    {

      /* Allocate space for output of single closest points.  */

      *gpar1 = *gpar2 = SISL_NULL;
      if ( (*gpar1 = newarray(*jpt, DOUBLE)) == SISL_NULL )  goto err101;
      if ( (*gpar2 = newarray(*jpt, DOUBLE)) == SISL_NULL )  goto err101;

      /* Copy single closest points to output arrays.  */

      for ( ki=0;  ki < (*jpt);  ki++ )
      {
	*(*gpar1 + ki) = ssing[2*ki];
	*(*gpar2 + ki) = ssing[2*ki+1];
      }
    }

    for ( kj=0;  kj < (*jcrv);  kj++ )
    {

      /* Correct parameter-arrays of closest curve.  */

      qcurve = *(*wcurve + kj);
      qcurve -> ipar1 = qcurve -> ipar2 = 1;

      /* Allocate space for new parameter arrays.  */

      spar1 = spar2 = SISL_NULL;
      if ( (spar1 = newarray(qcurve->ipoint, DOUBLE)) == SISL_NULL )  goto err101;
      if ( (spar2 = newarray(qcurve->ipoint, DOUBLE)) == SISL_NULL )  goto err101;

      /* Copy parameter-values to new arrays.  */

      for ( ki=0;  ki < qcurve->ipoint;  ki++ )
      {
	spar1[ki] = *(qcurve->epar1+2*ki);
	spar2[ki] = *(qcurve->epar1+2*ki+1);
      }

      freearray(qcurve->epar1);
      qcurve -> epar1 = spar1;
      qcurve -> epar2 = spar2;
    }
  }

  /* Closest points found.  */

  *jstat = 0;
  goto out;


  /* Error in space allocation.  */

 err101:
  *jstat = -101;
  s6err("s1955",*jstat,kpos);
  goto out;

  /* Error in input. Dimensions conflicing.  */

 err106:
  *jstat = -106;
  s6err("s1955",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1955",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied local variables.  */

  if ( qkreg1  &&  qkreg1 != pc1 )  freeCurve(qkreg1);
  if ( qkreg2  &&  qkreg2 != pc2 )  freeCurve(qkreg2);
  if (ssing)  freearray(ssing);
  if (sorigo) freearray(sorigo);
  if (qsdiff) freeSurf(qsdiff);

  return;
}
