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
 * $Id: s1022.c,v 1.3 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1022

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1022(double bottom_pos[], double bottom_axis[], double ellipse_ratio,
	   double axis_dir[], double cone_angle, double height,
	   SISLSurf **cone, int *stat)
#else
void s1022(bottom_pos, bottom_axis, ellipse_ratio, axis_dir,
	       cone_angle, height, cone, stat)
     double bottom_pos[];
     double bottom_axis[];
     double ellipse_ratio;
     double axis_dir[];
     double cone_angle;
     double height;
     SISLSurf  **cone;
     int    *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe a truncated cone as a NURBS. The cone can be
*              elliptic.
*
*
* INPUT      : bottom_pos    - Center point of the bottom
*              bottom_axis   - One of the bottom axis (major or minor)
*              ellipse_ratio - Ratio betwwen the other axis and bottom_axis
*              axis_dir      - Direction of the cone axis
*              cone_angle    - Angle between axis_dir and the cone at the end
*                              of bottom_axis, positive if the cone is sloping
*                              inwards.
*              height        - Height of the cone, can be negative.
*
*
* OUTPUT     :
*              stat          - status messages
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              cone          - Pointer to the cone produced
*
* METHOD     : The conic is made as a rules surface between two NURBS ellipses.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Added check
*              on input parameters to avoid NaN.
*
*********************************************************************
*/
{
  int kstat;           /* Status variable.                    */
  int kpos=0;            /* Position of error.                  */
  int ki;              /* Index in for loops.                 */
  int in1 = 9;         /* Number of vertices around the cone. */
  int in2 = 2;         /* Number of vertices along the cone.  */
  int ik1 = 3;         /* Order around the cone.              */
  int ik2 = 2;         /* Order along the cone.               */
  double et1[12];      /* Knot vector around the cone.        */
  double et2[4];       /* Knot vector along the cone.         */
  double rcoef[72];    /* Coefficients of the cone.           */
  int kind = 2;        /* Rational Bspline surface.           */
  double b2_axis[3];   /* The second bottom axis.             */
  double top_pos[3];   /* Center point of the top.            */
  double t1_axis[3];   /* The first top axis.                 */
  double t2_axis[3];   /* The second top axis.                */
  double norm_axis[3]; /* axis_dir normalized.                */
  double norm;         /* Vector norm.                        */
  double weight;       /* Rational weight.                    */


  /* Check input parameters. */

  if ( ellipse_ratio == 0.0 )
    goto err151;

  (void) s6length(bottom_axis, 3, &kstat);
  if (kstat == 0) goto err151;

  (void) s6length(axis_dir, 3, &kstat);
  if (kstat == 0) goto err151;


  /* Make the knot vectors. */

  for (ki = 0; ki < 12; ki++)
  {
    if (ki == 0 || ki == 1 || ki == 2)
      et1[ki] = (double)0.0;
    else if (ki == 3 || ki == 4)
      et1[ki] = PIHALF;
    else if (ki == 5 || ki == 6)
      et1[ki] = PI;
    else if (ki == 7 || ki == 8)
      et1[ki] = THREEPIHALF;
    else if (ki == 9 || ki == 10 || ki == 11)
      et1[ki] = TWOPI;
  }
  for (ki = 0; ki < 4; ki++)
  {
    if (ki == 0 || ki == 1)
      et2[ki] = (double)0.0;
    else if (ki == 2 || ki == 3)
      et2[ki] = fabs(height);
  }

  /* Make necessary vectors for the coefficients. */

  norm = s6norm(axis_dir, 3, norm_axis, &kstat);
  if (kstat < 0) goto error;
  s6crss(norm_axis, bottom_axis, b2_axis);
  for (ki = 0; ki < 3; ki++)
    b2_axis[ki] *= ellipse_ratio;

  for (ki = 0; ki < 3; ki++)
    top_pos[ki] = bottom_pos[ki] + height*norm_axis[ki];

  norm = s6length(bottom_axis, 3, &kstat);
  if (kstat < 0) goto error;
  norm = (double)1.0 - height*tan(cone_angle)/norm;
  for (ki = 0; ki < 3; ki++)
    t1_axis[ki] = norm*bottom_axis[ki];

  s6crss(norm_axis, t1_axis, t2_axis);
  for (ki = 0; ki < 3; ki++)
    t2_axis[ki] *= ellipse_ratio;

  /* Make the coefficients. */

  weight = (double)1.0/sqrt(2.0);
  for (ki = 0; ki < 3; ki++)
  {
    rcoef[ki] = bottom_pos[ki] + bottom_axis[ki];
    rcoef[4 + ki] = weight*(bottom_pos[ki] + bottom_axis[ki] + b2_axis[ki]);
    rcoef[8 + ki] = bottom_pos[ki] + b2_axis[ki];
    rcoef[12 + ki] = weight*(bottom_pos[ki] - bottom_axis[ki] + b2_axis[ki]);
    rcoef[16 + ki] = bottom_pos[ki] - bottom_axis[ki];
    rcoef[20 + ki] = weight*(bottom_pos[ki] - bottom_axis[ki] - b2_axis[ki]);
    rcoef[24 + ki] = bottom_pos[ki] - b2_axis[ki];
    rcoef[28 + ki] = weight*(bottom_pos[ki] + bottom_axis[ki] - b2_axis[ki]);
    rcoef[32 + ki] = rcoef[ki];

    rcoef[36 + ki] = top_pos[ki] + t1_axis[ki];
    rcoef[40 + ki] = weight*(top_pos[ki] + t1_axis[ki] + t2_axis[ki]);
    rcoef[44 + ki] = top_pos[ki] + t2_axis[ki];
    rcoef[48 + ki] = weight*(top_pos[ki] - t1_axis[ki] + t2_axis[ki]);
    rcoef[52 + ki] = top_pos[ki] - t1_axis[ki];
    rcoef[56 + ki] = weight*(top_pos[ki] - t1_axis[ki] - t2_axis[ki]);
    rcoef[60 + ki] = top_pos[ki] - t2_axis[ki];
    rcoef[64 + ki] = weight*(top_pos[ki] + t1_axis[ki] - t2_axis[ki]);
    rcoef[68 + ki] = rcoef[36 + ki];
  }
  for (ki = 3; ki < 72; ki += 4)
  {
    if (ki == 3 || ki == 11 || ki == 19 || ki == 27 || ki == 35
	|| ki == 39 || ki == 47 || ki == 55 || ki == 63 || ki == 71)
      rcoef[ki] = (double)1.0;
    else
      rcoef[ki] = weight;
  }

  (*cone) = newSurf(in1, in2, ik1, ik2, et1, et2, rcoef, kind, 3, 1);
  if ((*cone) == SISL_NULL) goto err101;

  *stat = 0;
  goto out;


  /* Error in curve allocation.  */

err101:
  *stat = -101;
  s6err("s1022",*stat,kpos);
  goto out;


  /* Error in input parameters - ellipse_ratio is zero, bottom_axis or
     axis_dir is of zero length.  */

err151:
  *stat = -151;
  s6err("s1022", *stat, kpos);
  goto out;


  /* Error in lower level routine. */

error:
  *stat = kstat;
  s6err("s1022", *stat, kpos);
  goto out;

out:
  return;
}
