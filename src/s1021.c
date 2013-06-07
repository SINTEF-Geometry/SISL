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
 * $Id: s1021.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1021

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1021(double bottom_pos[], double bottom_axis[], double ellipse_ratio, 
	 double axis_dir[], double height, SISLSurf **cyl, int *stat)
#else
void s1021(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	      height, cyl, stat) 
     double bottom_pos[];
     double bottom_axis[];
     double ellipse_ratio;
     double axis_dir[];
     double height;
     SISLSurf  **cyl;
     int    *stat;      
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a truncated cylinder as a NURBS. The cylinder 
*              can be elliptic.
*             
*
* INPUT      : bottom_pos    - Center point of the bottom
*              bottom_axis   - One of the bottom axis (major or minor)
*              ellipse_ratio - Ratio betwwen the other axis and bottom_axis
*              axis_dir      - Direction of the cylinder axis
*              height        - Height of the cone, can be negative.
*
*
* OUTPUT     : 
*              stat          - status messages  
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              cyl           - Pointer to the cylinder produced
*
* METHOD     : The cylinder is made as a rules surface between two NURBS ellipses.
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
   int kstat;                      /* Status variable.   */
   int kpos=0;                     /* Position of error. */
   double cone_angle = (double)0.; /* No conical slope.  */

   s1022(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	     cone_angle, height, cyl, &kstat);
   if (kstat < 0) goto error;
  
   *stat = 0;
   goto out;
  
   /* Error in lower level routine. */
      
   error:
      *stat = kstat;
      s6err("s1021", *stat, kpos);
      goto out;
  
   out:  
      return;
}
    
