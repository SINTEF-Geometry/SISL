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
 * $Id: s1011.c,v 1.3 2001-03-19 15:58:40 afr Exp $
 *
 */


#define S1011

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1011(double start_pos[], double top_pos[], double end_pos[],
 	   double shape, int dim, SISLCurve **arc_seg, int *stat)
#else
void s1011(start_pos, top_pos, end_pos, shape, dim, arc_seg, stat)
     double start_pos[];
     double top_pos[];
     double end_pos[];
     double shape;
     int    dim;
     SISLCurve  **arc_seg;
     int    *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe a conic arc as a NURBS. The arc is given by
*              position at start, shoulder point and end, and a shape factor.
*
*
* INPUT      : start_pos - Start point of segment
*              top_pos   - Shoulder point of segment. This is the intersection
*                          point of the tangents in start_pos and end_pos
*              end_pos   - End point of segment
*              shape     - Shape factor >= 0
*                           shape < 0.5 a ellipse
*                           shape = 0.5 a parabola
*                           shape > 0.5 a hyperbola
*                           shape >= 1  The start and end point lies on
*                                       different branches of the hyperbola.
*                                       We want a single arc segment, therefore
*                                       if shape>=1, shape is put to 0.999999.
*              dim       - The dimension of the curve to be produced
*
*
*
*
* OUTPUT     :
*              stat      - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              arc_seg   - Pointer to the curve produced
*
* METHOD     : The conic is made as a rational B-spline curve according
*              to the following formula:
*
*                     p0 (1-t)(1-t) + 2 s/(1-s) t(1-t) pt + tt p1
*              p(t) = -------------------------------------------
*                         (1-t)(1-t) + 2 s/(1-s) t(1-t) + tt
*
*              where p0 is the start point, pt is the shoulder point, p1 is
*              the end point and s is the shape factor.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93 (Based on s1385,
*              written by Tor Dokken)
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994. Give error if
*              shape < zero.
*********************************************************************
*/
{
  int kpos = 0;       /* Error position.         */
  int ki;             /* Index in for loop.      */
  int in = 3;         /* Number of coefficients. */
  int ik = 3;         /* Order of the curve.     */
  int rdim = dim + 1; /* Rational dimension.     */
  double et[6];       /* Knot vector.            */
  double rcoef[12];   /* Rational coefficients.  */
  int kind = 4;       /* Rational Bezier curve.  */
  double weight;      /* Rational weight.        */


  /* Make sure we get a single arc segment and positive weights. */

  if ( shape >= (DOUBLE)1.0 ) shape = (DOUBLE)0.9999999;
  else if ( shape < (DOUBLE)0.0 )  goto err151;


  /* Make the data needed for curve generation. */

  for ( ki=0;  ki < ik;  ki++ )
  {
    et[ki]      = DZERO;
    et[ik + ki] = (DOUBLE)1.0;
  }

  weight = shape/((DOUBLE)1.0 - shape);

  for ( ki=0;  ki < dim;  ki++ )
  {
    rcoef[ki]          = start_pos[ki];
    rcoef[rdim + ki]   = weight*top_pos[ki];
    rcoef[2*rdim + ki] = end_pos[ki];
  }
  rcoef[dim] =          (DOUBLE)1.0;
  rcoef[dim + rdim] =   weight;
  rcoef[dim + 2*rdim] = (DOUBLE)1.0;

  (*arc_seg) = newCurve(in, ik, et, rcoef, kind, dim, 1);
  if ((*arc_seg) == SISL_NULL) goto err101;

  *stat = 0;
  goto out;


  /* Error in curve allocation.  */

err101:
  *stat = -101;
  s6err("s1011", *stat, kpos);
  goto out;


  /* Error in input parameters - shape is negative.  */

err151:
  *stat = -151;
  s6err("s1011", *stat, kpos);
  goto out;

out:
  return;
}
