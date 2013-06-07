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
 * $Id: s1920.c,v 1.3 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1920

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1920(SISLCurve *pc1,double edir[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1920(pc1,edir,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
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
* PURPOSE    : Find the extremal points/intervals of the curve pc1 in
*              the direction edir.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
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
*                       interval of the curve. The points lie continuous.
*                       Extremal curves are stored in wcurve.
*              jcrv   - Number of extremal curves.
*              wcurve - Array containing descriptions of the extremal
*                       curves. The curves are only described by points
*                       in the parameter interval. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The scalar-product of the coefficients of the curve with
*              the vector edir is calculated giving a curve with dimension
*              1. Then the maxima of this curve is calculated.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1880 - Put extremal points/intervals into output format.
*              s1161 - Find maxima of one-dimensional curve.
*              s6scpr - Scalar-product between two vectors.
*              newCurve   - Create new curve.
*              newObject  - Create new object.
*              newIntpt   - Create new extremal point.
*              freeObject - Free the space occupied by a given object.
*              freeIntdat - Free space occupied by an extremal list.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-10.
* REVISED BY : Mike Floater, SI, 1991-04 for a rational curve.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input curve with periodic basis correctly.
*********************************************************************
*/
{
  int ikind;               /* Type of curve pc1 is                     */
  int kstat = 0;           /* Local status variable.                     */
  int kpos = 0;            /* Position of error.                         */
  int ki;                  /* Counter.                                   */
  int kn;                  /* Number of vertices of curve.               */
  int kk;                  /* Order of curve.                            */
  double tmax;             /* Estimate of maximal value of 1-dim. curve. */
  double *st;              /* Pointer to knotvector of curve.            */
  double *scoef;           /* Pointer to vertices of curve.              */
  double *sc = SISL_NULL;       /* Pointer to vertices of curve in maxima
			      calculation.                               */
  double *spar = SISL_NULL;     /* Values of extremal points in the parameter
			      area of the second object. Empty in this case.*/
  double *s1,*s2,*sstop;   /* Pointers used to traverse double-arrays.   */
  SISLIntdat *qintdat = SISL_NULL;  /* Maximum results.                     */
  SISLCurve *qc = SISL_NULL;        /* Pointer to curve in maxima calculation.    */
  SISLObject *qo1 = SISL_NULL;      /* Pointer to object in maxima calculation. */
  SISLCurve *qkreg = SISL_NULL;     /* Input curve with k-regularity ensured. */


  *jpt = 0;
  *jcrv = 0;


  /* Ensure k-regular basis. */

  if ( pc1 -> cuopen == SISL_CRV_PERIODIC )
  {
    /* Cyclic (periodic) basis. */

    make_cv_kreg(pc1, &qkreg, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg = pc1;


  /* Check dimension.  */

  if ( qkreg -> idim != idim )  goto err106;

  /* Describe curve with local variables.  */

  kn = qkreg -> in;
  kk = qkreg -> ik;
  st = qkreg -> et;
  ikind = qkreg -> ikind;

  if ( ikind == 2 || ikind == 4 )
  {
      scoef = qkreg -> rcoef;
      /* Allocate space for coeffecients of new curve.  */

      if ( (sc = newarray(2*kn, DOUBLE)) == SISL_NULL )  goto err101;

      /* Compute scalar-product of curve-vertices and direction vector. */
      /* Copy over weights. */

      for ( s1=scoef, s2=sc, sstop=s2+2*kn;  s2<sstop;  s1+=idim+1, s2+=2 )
      {
          *s2 = s6scpr(s1, edir, idim);
          *(s2+1) = *(s1+idim);
      }
  }
  else
  {
      scoef = qkreg -> ecoef;
      /* Allocate space for coeffecients of new curve.  */

      if ( (sc = newarray(kn, DOUBLE)) == SISL_NULL )  goto err101;

      /* Compute scalar-product of curve-vertices and direction vector. */

      for ( s1=scoef, s2=sc, sstop=s2+kn;  s2<sstop;  s1+=idim, s2++ )
      {
          *s2 = s6scpr(s1, edir, idim);
      }
  }


  /* Create new curve.  */

  qc = newCurve(kn, kk, st, sc, qkreg->ikind, 1, 1);
  if ( qc == SISL_NULL )  goto err101;

  /* Create new object and connect curve to object.  */

  qo1 = newObject(SISLCURVE);
  if ( qo1 == SISL_NULL )  goto err101;
  qo1 -> c1 = qc;

  tmax = -(double)HUGE;

  /* Find maxima. */
  s1161(qo1, &tmax, aepsge, &qintdat, &kstat);
  if ( kstat < 0 )  goto error;

  if (qintdat)
  {

    /* Express maximal points/intervals on output format.  */
    s1880(1, 0, &qintdat->ipoint, qintdat->vpoint,
	  &qintdat->ilist, qintdat->vlist,
	  jpt, gpar, &spar, jcrv, wcurve, &kstat);
    if ( kstat < 0 )  goto error;

    /* Handle periodicity (remove points at periodic curve ends). */
    if ( *jpt > 1 && idim > 1 && pc1->cuopen == SISL_CRV_PERIODIC )
    {
      for ( ki=0; ki < (*jpt); ki++ )
      {
	if ( (*gpar)[ki] == pc1->et[pc1->in] )
	{
	  (*jpt)--;
	  (*gpar)[ki] = (*gpar)[*jpt];
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
  s6err("s1920", *jstat, kpos);
  goto out;

  /* Dimensions conflicting.  */

 err106:
  *jstat = -106;
  s6err("s1920", *jstat, kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1920", *jstat, kpos);
  goto out;

 out:

  /* Free allocated space.  */

  if ( qkreg && qkreg != pc1 ) freeCurve(qkreg);
  if (sc) freearray(sc);
  if (spar) freearray(spar);
  if (qintdat) freeIntdat(qintdat);
  if (qo1) freeObject(qo1);


  return;
}
