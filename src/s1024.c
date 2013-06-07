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
 * $Id: s1024.c,v 1.3 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1024

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1024(double center[], double axis[], double equator[], double minor_radius,
          int start_minor, int end_minor, int numb_major, SISLSurf **torus,
	  int *stat)
#else
void s1024(center, axis, equator, minor_radius, start_minor, end_minor,
              numb_major, torus, stat)
     double center[];
     double axis[];
     double equator[];
     double minor_radius;
     int start_minor;
     int end_minor;
     int numb_major;
     SISLSurf  **torus;
     int    *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe octants of a torus as a NURBS. This can also
*              be used to describe the complete torus.
*
*
* INPUT      : center        - Center point of the torus
*              axis          - Normal to the torus plane
*              equator       - Vector from center to start point on the major
*                              circle
*              minor_radius  - Radius of the minor circle
*              start_minor   - Start quadrant on the minor circle (1,2,3 or 4).
*                              This is counted clockwise from the extremum in
*                              the direction of axis
*              end_minor     - End quadrant on the minor circle (1,2,3 or 4)
*                              This is counted clockwise from the extremum in
*                              the direction of axis
*              numb_major    - Number of quadrants on the major circle
*                              (1,2,3 or 4). This is counted counter-clockwise
*                              from equator.
*
*
* OUTPUT     :
*              stat          - status messages
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              torus         - Pointer to the torus produced
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93
* REVISED BY : Christophe Rene Birkeland, SINTEF, Oslo, May 1993
*              (Testing on SISL_NULL-pointers atfer allocations)
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov, 1994. Added check on
*              input parameters to avoid NaN results.
*
*********************************************************************
*/
{
  int kstat;            /* Status variable.                               */
  int kpos=0;           /* Position of error.                             */
  int ki, kj, kl;       /* Indexes in for loops.                          */
  int in1;              /* Number of vertices around the minor circle.    */
  int in2;              /* Number of vertices around the major circle.    */
  int ik1 = 3;          /* Order around the minor circle.                 */
  int ik2 = 3;          /* Order around the major circle.                 */
  double *et1 = SISL_NULL;   /* Knot vector around the minor circle.           */
  double *et2 = SISL_NULL;   /* Knot vector around the major circle.           */
  double *rcoef = SISL_NULL; /* Coefficients of the torus.                     */
  int kind = 2;         /* Rational Bspline surface.                      */
  double weight;        /* Rational weight.                               */
  double major_radius;  /* Radius of the major circle.                    */
  double norm;          /* Length of vectors.                             */
  double x_axis[3];     /* Normalized equator.                            */
  double y_axis[3];     /* Normal to equator in torus plane, normalized.  */
  double z_axis[3];     /* axis normalized.                               */
  double w1, w2;        /* Rational weights in both directions.           */
  double x_comp;        /* Component in local x direction.                */
  double y_comp;        /* Component in local y direction.                */
  double z_comp;        /* Component in local z direction.                */
  int start_coef;       /* Number of starting vertice in minor circle.    */
  int end_coef;         /* Number of second last vertice in minor circle. */


  /* Do necessary initiation and allocation. */

  *torus = SISL_NULL;

  if ( start_minor < 1  ||  start_minor > 4  ||
       end_minor < 1  ||  end_minor > 4  ||
       numb_major < 1  || numb_major > 4 )  goto err151;

  start_coef = 2*(start_minor - 1);
  end_coef   = 2*end_minor + 1;
  weight = (DOUBLE)1.0/sqrt(2.0);
  in1 = 1 + 2*(end_minor + 1 - start_minor);
  in2 = 1 + 2*numb_major;

  major_radius = s6length(equator, 3, &kstat);
  if ( kstat < 0 ) goto error;
  else if ( kstat == 0 ) goto err151;

  for ( ki=0;  ki < 3;  ki++ )
    x_axis[ki] = equator[ki]/major_radius;
  norm = s6length(axis, 3, &kstat);
  if ( kstat < 0 ) goto error;
  else if ( kstat == 0 )  goto err151;

  for ( ki=0;  ki < 3;  ki++ )
    z_axis[ki] = axis[ki]/norm;
  s6crss(z_axis, x_axis, y_axis);

  if ( (et1 = newarray(in1 + ik1, DOUBLE)) == SISL_NULL ) goto err101;
  if ( (et2 = newarray(in2 + ik2, DOUBLE)) == SISL_NULL ) goto err101;
  if ( (rcoef = newarray(4*in1*in2, DOUBLE)) == SISL_NULL ) goto err101;

  /* Initiate the knot vectors. */

  for ( ki=0;  ki < ik1;  ki++ )
    et1[ki] = (DOUBLE)0.;
  for ( ki=0;  ki < (end_minor + 1 - start_minor);  ki++ )
  {
    et1[ik1 + 2*ki]     = (ki + 1)*PIHALF;
    et1[ik1 + 2*ki + 1] = (ki + 1)*PIHALF;
  }
  et1[in1 + ik1 - 1] = (end_minor + 1 - start_minor)*PIHALF;

  for ( ki=0;  ki < ik2;  ki++ )
    et2[ki] = (DOUBLE)0.;
  for ( ki=0;  ki < numb_major;  ki++ )
  {
    et2[ik2 + 2*ki]     = (ki + 1)*PIHALF;
    et2[ik2 + 2*ki + 1] = (ki + 1)*PIHALF;
  }
  et2[in2 + ik2 - 1] = numb_major*PIHALF;

  /* Initiate the coefficient vector. */

  for ( ki=0;  ki < in2;  ki++ )
  {
    if ( ki == 1  ||  ki == 3  ||  ki == 5  ||  ki == 7 )
      w2 = weight;
    else
      w2 = (DOUBLE)1.;

    for ( kj=start_coef;  kj < end_coef;  kj++ )
    {
      if ( kj == 1  ||  kj == 3  ||  kj == 5  ||  kj == 7 )
	w1 = w2*weight;
      else
	w1 = w2;

      if ( ki == 0  ||  ki == 1  ||  ki == 7  ||  ki == 8 )
      {
	if ( kj == 1  ||  kj == 2  ||  kj == 3 )
	  x_comp = major_radius + minor_radius;
	else if ( kj == 5  ||  kj == 6  ||  kj == 7 )
	  x_comp = major_radius - minor_radius;
	else
	  x_comp = major_radius;
      }
      else if ( ki == 3  ||  ki == 4  ||  ki == 5 )
      {
	if ( kj == 1  ||  kj == 2  ||  kj == 3 )
	  x_comp = - major_radius - minor_radius;
	else if ( kj == 5  ||  kj == 6  ||  kj == 7 )
	  x_comp = - major_radius + minor_radius;
	else
	  x_comp = - major_radius;
      }
      else
	x_comp = (DOUBLE)0.;

      if ( ki == 1  ||  ki == 2  ||  ki == 3 )
      {
	if ( kj == 1  ||  kj == 2  ||  kj == 3 )
	  y_comp = major_radius + minor_radius;
	else if ( kj == 5  ||  kj == 6  ||  kj == 7 )
	  y_comp = major_radius - minor_radius;
	else
	  y_comp = major_radius;
      }
      else if ( ki == 5  ||  ki == 6  ||  ki == 7 )
      {
	if ( kj == 1  ||  kj == 2  ||  kj == 3 )
	  y_comp = - major_radius - minor_radius;
	else if ( kj == 5  ||  kj == 6  ||  kj == 7 )
	  y_comp = - major_radius + minor_radius;
	else
	  y_comp = - major_radius;
      }
      else
	y_comp = (DOUBLE)0.;

      if ( kj == 0  ||  kj == 1  ||  kj == 7  ||  kj == 8 )
	z_comp = minor_radius;
      else if ( kj == 3  ||  kj == 4  ||  kj == 5 )
	z_comp = - minor_radius;
      else
	z_comp = (DOUBLE)0.;

      for ( kl=0;  kl < 3;  kl++ )
	rcoef[4*(ki*in1 + kj - start_coef) + kl] =
	  w1*(center[kl] + x_comp*x_axis[kl] +
	      y_comp*y_axis[kl] + z_comp*z_axis[kl]);
      rcoef[4*(ki*in1 + kj - start_coef) + 3] = w1;

    }
  }

  (*torus) = newSurf(in1, in2, ik1, ik2, et1, et2, rcoef, kind, 3, 1);
  if ( (*torus) == SISL_NULL ) goto err101;

  *stat = 0;
  goto out;


  /* Error in curve allocation.  */

err101:
  *stat = -101;
  s6err("s1024", *stat, kpos);
  goto out;


  /* Error in input parameters - ellipse_ratio is zero, bottom_axis or
     axis_dir is of zero length.  */

err151:
  *stat = -151;
  s6err("s1024", *stat, kpos);
  goto out;


  /* Error in lower level routine. */

error:
  *stat = kstat;
  s6err("s1024", *stat, kpos);
  goto out;

out:

  if ( et1 != SISL_NULL ) freearray(et1);
  if ( et2 != SISL_NULL ) freearray(et2);
  if ( rcoef != SISL_NULL ) freearray(rcoef);

  return;
}
