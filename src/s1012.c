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
 * $Id: s1012.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define S1012

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1012(double start_pos[], double axis_pos[], double axis_dir[],
 	   double frequency, int numb_quad, int counter_clock, 
	   SISLCurve **helix, int *stat)
#else
void s1012(start_pos, axis_pos, axis_dir, frequency, numb_quad, 
		counter_clock, helix, stat) 
     double start_pos[];
     double axis_pos[];
     double axis_dir[];
     double frequency;
     int numb_quad;
     int counter_clock;
     SISLCurve  **helix;
     int    *stat;      
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a truncated helix  as a NURBS.
*             
*
* INPUT      : start_pos     - Start position on the helix
*              axis_pos      - Point on the helix axis
*              axis_dir      - Direction of the helix axis
*              frequency     - The length along the helix axis for one period 
*                              of revolution
*              numb_quad     - Number of quadrants in the helix
*              counter_clock - Flag for direction of revolution:
*                              = 0 : clockwise
*                              = 1 : counter_clockwise
*
*
*
* OUTPUT     : 
*              stat          - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              helix         - Pointer to the helix produced
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
*
*********************************************************************
*/
{
   int kstat;            /* Status variable.                 */
   int kpos = 0;         /* Error position.                  */
   int ki, kj;           /* Index in for loop.               */
   int mod_count;        /* Modulo index.                    */
   int in;               /* Number of vertices in the helix. */
   int ik = 3;           /* Order of the helix.              */
   int kind = 2;         /* Rational B-spline curve.         */
   int dim = 3;          /* Dimension of geometry.           */
   double *et = SISL_NULL;    /* Knot vector.                     */
   double *rcoef = SISL_NULL; /* Vertices.                        */
   double norm;          /* Norm of vector.                  */
   double radius;        /* Radius of the helix.             */
   double origo[3];      /* Local origo.                     */
   double x_axis[3];     /* Normalized local x_axis.         */
   double y_axis[3];     /* Normalized local y_axis.         */
   double z_axis[3];     /* Normalized local z_axis.         */
   double x_comp;        /* Component on x_axis.             */
   double y_comp;        /* Component on y_axis.             */
   double z_comp;        /* Component on z_axis.             */
   double weight;        /* Rational weight.                 */
   double w1;            /* Rational weight.                 */

 
   weight = (double)1.0/sqrt(2.0);  
   
   /* Allocate space. */
   
   in = 1 + 2*numb_quad;
   et = newarray(in + ik, DOUBLE);
   rcoef = newarray(4*in, DOUBLE);
   
   /* Initiate knot vector. */
   
   for (ki = 0; ki < ik; ki++)
      et[ki] = (double)0.;
   for (ki = 0; ki < numb_quad; ki++)
   {
      et[3 + 2*ki]     = (1 + ki)*PIHALF;
      et[3 + 2*ki + 1] = (1 + ki)*PIHALF;
   }
   et[ik + in - 1] = numb_quad*PIHALF;
   
   /* Initiate coefficient vector. */
   
   norm = s6norm(axis_dir, 3, z_axis, &kstat);
   if (kstat < 0) goto error;
		  
   s6diff(start_pos, axis_pos, 3, x_axis);
   norm = s6scpr(x_axis, z_axis, 3);
   for (ki = 0; ki < 3; ki++)
   {
      origo[ki] = axis_pos[ki] + norm*z_axis[ki];
      x_axis[ki] = start_pos[ki] - origo[ki];
   }
   radius = s6norm(x_axis, 3, x_axis, &kstat);
   
   if (counter_clock == 0)
      s6crss(x_axis, z_axis, y_axis);
   else
      s6crss(z_axis, x_axis, y_axis);
      
   mod_count = 0;
   for (ki = 0; ki < in; ki++)
   {
      if (mod_count == 1 || mod_count == 3 
	  || mod_count == 5 || mod_count == 7)
	 w1 = weight;
      else
         w1 = (double)1.;
      
      if (mod_count == 0 || mod_count == 1 || mod_count == 7)
	 x_comp = radius;
      else if (mod_count == 3 || mod_count == 4 || mod_count == 5)
	 x_comp = - radius;
      else
         x_comp = (double)0.;

      if (mod_count == 1 || mod_count == 2 || mod_count == 3)
	 y_comp = radius;
      else if (mod_count == 5 || mod_count == 6 || mod_count == 7)
	 y_comp = - radius;
      else
         y_comp = (double)0.;
	 
      z_comp = ki*frequency/8;
       
      for (kj = 0; kj < 3; kj++)
	rcoef[4*ki + kj] = w1*(origo[kj] + x_comp*x_axis[kj]
			+ y_comp*y_axis[kj] + z_comp*z_axis[kj]);
      rcoef[4*ki + 3] = w1;
      
      mod_count++;
      if (mod_count == 8) mod_count = 0;
   }
   
   (*helix) = newCurve(in, ik, et, rcoef, kind, dim, 1);
   if (et != SISL_NULL) freearray(et);
   if (rcoef != SISL_NULL) freearray(rcoef);
   if ((*helix) == SISL_NULL) goto err101;
  
   *stat = 0;
   goto out;
  
   /* Error in curve allocation.  */
  
   err101: 
      *stat = -101;
      s6err("s1012",*stat,kpos);
      goto out;
      
   /* Error in lower level routine. */
      
   error:
      *stat = kstat;
      s6err("s1012", *stat, kpos);
      goto out;
    
   out:  
      return;
}
    
