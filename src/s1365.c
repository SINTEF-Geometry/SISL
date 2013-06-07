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
 * $Id: s1365.c,v 1.5 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1365

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1365(SISLSurf *ps, double aoffset, double aepsge, double amax,
	   int idim, SISLSurf **rs, int *jstat)
#else
void s1365(ps,aoffset,aepsge,amax,idim,rs,jstat)
     SISLSurf   *ps;
     double aoffset;
     double aepsge;
     double amax;
     int    idim;
     SISLSurf   **rs;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To create a B-spline approximating the offset surface of
*              a NURBS surface.
*
*
*
* INPUT      : ps     - The input NURBS surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              aepsge - Maximal deviation allowed between true offset surface
*                       and the approximated offset surface.
*              amax   - Maximal stepping length. Is negleceted if amax<=aepsge
*                       If amax==0 then a maximal step length of the longest
*                       box side is used.
*              idim   - The dimension of the space (2 or 3).
*
* OUTPUT     :
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer the approximated offset surface
*
* METHOD     : The 4 edge curves of the surface are extracted. Offset curves
*              of these 4 edge curves are approximated and a common
*              basis for the two pairs of opposite offset curves is calculated.
*              Vertices are recomputed.
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*-
* CALLS      : s1435     - Pick a given edge-curve of a B-spline surface.
*              s1360     - Approximate the offset curve of a NURBS curve.
*              s1933     - Put a set of curves on a common basis.
*              s1366     - Create a B-spline surface approximating the offset
*                          surface of a B-spline surface
*
* WRITTEN BY : Per Evensen,  SI, 89-3.
* REWISED BY : Per Evensen,  SI, 90-9; Corrected start/end parameter values of
*                                      common curves.
*
*********************************************************************
*/
{
  SISLCurve *pc[4];
  SISLCurve *rc[4],*rc13[2],*rc24[2];
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int kdim;          /* Dimension of the space in which the surface lies.*/
  int knbcrv = 2;    /* Number of curves in set. */
  /*  int kopen = 1;      Flag telling that the resulting surface should
			be open in both parameter directions. */
  int kn13,kord13;   /* Number of vertices and order of edge curves along
			1. parameter direction. */
  int kn24,kord24;   /* Number of vertices and order of edge curves along
			2. parameter direction. */
  int nder[4];       /* Number of edges along each surface edge.         */
  double sp[4];      /* Parameter value of edge in constatnt direction.  */
  double toffset = (double)0.0; /* Local offset value for extraction of edge-
				   curves. */
  double snorm[3];   /* Local normal vector for extraction of edge-curves. */
  double tstart1,tend1; /* Endpoints of parameter interval in first
			   direction.                                     */
  double tstart2,tend2; /* Endpoints of parameter interval in second
			   direction.                                     */
  double *sknot13=SISL_NULL;/* Pointer to common knot-vector of edge curves along
			  1. parameter direction. */
  double *sknot24=SISL_NULL;/* Pointer to common knot-vector of edge curves along
			  2. parameter direction. */
  int  kk;              /* Loop controller. */

  /* Initialization of variables */
  kdim = ps -> idim;
  for (kk=0; kk<4; kk++)
  {
     nder[kk] = 1;
     pc[kk] = SISL_NULL;
     rc[kk] = SISL_NULL;
  }
  for (kk=0; kk<3; kk++) snorm[kk] = DZERO;

  /* Fetch the 4 edge-curves of surface */

  for (kk=0; kk<4; kk++)
    {
      s1435(ps,kk,&pc[kk],&sp[kk],&kstat);
      if (kstat < 0) goto error;
    }

  /* Create a B-spline curve approximating the offset curve of the 4 edges */

  for (kk=0; kk<4; kk++)
    {
      s1360(pc[kk],toffset,aepsge,snorm,amax,kdim,&rc[kk],&kstat);
      if (kstat<0) goto error;
    }

  /* Rearrange the pointers to the 4 edge curves. */

  rc13[0] = rc[0];
  rc13[1] = rc[2];
  rc24[0] = rc[1];
  rc24[1] = rc[3];

  /* Fetch endpoints of parameter intervals.  */

  tstart1 = *(ps->et1 + ps->ik1 - 1);
  tend1 = *(ps->et1 + ps->in1);
  tstart2 = *(ps->et2 + ps->ik2 - 1);
  tend2 = *(ps->et2 + ps->in2);

  /* Put the edge curves along 1. parameter direction into common basis. */

  s1933(knbcrv,rc13,tstart1,tend1,&sknot13,&kn13,&kord13,&kstat);
  if (kstat<0) goto error;

  /* Put the edge curves along 2. parameter direction into common basis. */

  s1933(knbcrv,rc24,tstart2,tend2,&sknot24,&kn24,&kord24,&kstat);
  if (kstat<0) goto error;

  /* Create a B-spline surface approximating the offset surface of a B-spline
     surface. */

  s1366(ps,aoffset,aepsge,amax,idim,sknot13,kn13,kord13,
	sknot24,kn24,kord24,rs,&kstat);
  if (kstat<0) goto error;

  /* Surface approximated. */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1365",*jstat,kpos);
  goto out;

  out:
     for (kk=0; kk<4; kk++)
     {
	if (pc[kk] != SISL_NULL) freeCurve(pc[kk]);
	if (rc[kk] != SISL_NULL) freeCurve(rc[kk]);
     }
     if (sknot13 != SISL_NULL) freearray(sknot13);
     if (sknot24 != SISL_NULL) freearray(sknot24);

     return;
}
