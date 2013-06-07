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
 * $Id: s1245.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1245

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1245(double coef[], int ik, int dim, double point[], 
      double local_tol, int depth, double weight[], double *area, 
      double *moment, int *stat)
#else
void s1245(coef, ik, dim, point, local_tol, depth, weight, area, moment, stat)
     double     coef[];
     int        ik;
     int        dim;
     double     point[];
     double     local_tol;
     int        depth;
     double     weight[];
     double     *area;
     double     *moment;
     int        *stat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the weight point and rotational momentum of
*              an area between a 2D Bezier segment and a 2D point. The area 
*              is also calculated.
*              When the curve is rotating counter-clockwise around the
*              point, the area contribution is positive.
*              When the curve is rotating clockwise around the point,
*              the area contribution is negative.
*
* INPUT      : coef   - Coefficients of the Bezier segment.
*              ik     - Order of the segment in question.
*              dim    - Dimension of geometry (must be 2).
*              point  - The reference point.
*              local_tol- The current tolerance.
*              depth  - Depth of recursion.
*
*
* OUTPUT     : weight - Weight point.
*              area   - Area.
*              moment - Rotational momentum.
*              stat   - Status messages  
*                       = 0 : OK.
*                       < 0 : Error.
*                       > 0 : Warning.
*
*                                  
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Oslo, Norway, 12-94.
*
*********************************************************************
*/
{
   int ki, kj;                /* Index in for loop.         */
   int pos = 0;               /* Position of error.         */
   int index;                 /* Array index.               */
   double add_area;           /* Addition to area.          */
   double add_moment;         /* Addition to area.          */
   double left_area;          /* Left area in recursion.    */
   double right_area;         /* Right area in recursion.   */
   double left_moment;        /* Left moment in recursion.  */
   double right_moment;       /* Right moment in recursion. */
   double vec1[2];            /* Utility vector.            */
   double vec2[2];            /* Utility vector.            */
   double vec3[2];            /* Utility vector.            */
   double vec4[2];            /* Utility vector.            */
   double left_weight[2];     /* Left weight in recursion.  */
   double right_weight[2];    /* right weight in recursion. */
   double* left_coef = SISL_NULL;  /* Left coefficients.         */
   double* right_coef = SISL_NULL; /* Left coefficients.         */


   /* Check input. */

   if (dim != 2)
      goto err106;
   
   /* Check if this is a straight line. */

   if (ik < 3)
   {
      
      /* Straight line. */
      
      for (ki = 0; ki < 2; ki++)
      {
	 vec1[ki] = coef[ki] - point[ki];
         vec2[ki] = coef[2*(ik - 1) + ki] - point[ki];
	 vec3[ki] = (coef[ki] + coef[2*(ik - 1) + ki] + point[ki])/3.;
	 vec4[ki] = coef[2*(ik - 1) + ki] - coef[ki]; 
      }

      *area = (vec1[0]*vec2[1] - vec1[1]*vec2[0])/2.;
      *moment = ((vec1[0]*vec1[0] + vec1[1]*vec1[1])/4. +
		 (vec1[0]*vec4[0] + vec1[1]*vec4[1])/4. +
		 (vec4[0]*vec4[0] + vec4[1]*vec4[1])/12.)/
	         fabs(vec1[0]*vec4[1] - vec1[1]*vec4[0]);
      if (*area < 0)
	 *moment = - (*moment);
      
      weight[0] = (*area)*vec3[0];
      weight[1] = (*area)*vec3[1];
   }
   else
   {
      
      /* Not straight line. */
      
      *area = 0.0;
      *moment = 0.0;
      weight[0] = 0.0;
      weight[1] = 0.0;

      /* Do a calculation on all the coefficients. */

      for (ki = 1; ki < ik; ki++)
      {

	 for (kj = 0; kj < 2; kj++)
         {
	    vec1[kj] = coef[2*(ki - 1) + kj] - point[kj];
            vec2[kj] = coef[2*ki + kj] - point[kj];
    	    vec3[kj] = (coef[2*(ki - 1) + kj] + coef[2*ki + kj] +
			point[kj])/3.;
	    vec4[kj] = coef[2*ki + kj] - coef[2*(ki - 1) + kj]; 
         }

         add_area = (vec1[0]*vec2[1] - vec1[1]*vec2[0]);
         add_moment = ((vec1[0]*vec1[0] + vec1[1]*vec1[1])/4. +
		 (vec1[0]*vec4[0] + vec1[1]*vec4[1])/4. +
		 (vec4[0]*vec4[0] + vec4[1]*vec4[1])/12.)/
	         fabs(vec1[0]*vec4[1] - vec1[1]*vec4[0]);
         if (add_area < 0)
	    add_moment = - add_moment;
      
         weight[0] += add_area*vec3[0];
         weight[1] += add_area*vec3[1];	
	 *area += add_area;
	 *moment += add_moment;
      }
      
      /* Do a calculation on the first and last coefficient. */
      
      for (kj = 0; kj < 2; kj++)
      {
	 vec1[kj] = coef[kj] - point[kj];
         vec2[kj] = coef[2*(ik - 1) + kj] - point[kj];
    	 vec3[kj] = (coef[kj] + coef[2*(ik - 1) + kj] + point[kj])/3.;
	 vec4[kj] = coef[2*(ik - 1) + kj] - coef[kj]; 
      }

      add_area = (vec1[0]*vec2[1] - vec1[1]*vec2[0]);
      add_moment = ((vec1[0]*vec1[0] + vec1[1]*vec1[1])/4. +
		 (vec1[0]*vec4[0] + vec1[1]*vec4[1])/4. +
		 (vec4[0]*vec4[0] + vec4[1]*vec4[1])/12.)/
	         fabs(vec1[0]*vec4[1] - vec1[1]*vec4[0]);
      if (add_area < 0)
	 add_moment = - add_moment;
      
      weight[0] += add_area*vec3[0];
      weight[1] += add_area*vec3[1];	
      *area += add_area;
      *moment += add_moment;  
      
      /* Check the deviation between them. */
      
      if (fabs(*area) < REL_COMP_RES)
      {
	 
	 /* No contribution. */
	 
	 weight[0] = 0.0;
	 weight[1] = 0.0;
	 *area = 0.0;
	 *moment = 0.0;
      }
      
      else if (fabs(2*add_area - *area)/fabs(*area) < local_tol || depth > 20)
      {
	 
	 /* Good enough. */
	 
	 weight[0] /= 4.;
	 weight[1] /= 4.;
	 *area /= 4.;
	 *moment /= 2.;
      }
      
      else
      {
	 
	 /* Not good enough, we have to subdivide. */
	 
	 left_coef = newarray(2*ik, double);
	 right_coef = newarray(2*ik, double);
	 
	 for (ki = 0; ki < 2*ik; ki++)
	 {
	    left_coef[ki] = coef[ki]; 
	    right_coef[ki] = coef[ki];
	 }
	 
	 for (ki = 1; ki < ik; ki++)
	 {
	    for (kj = ki; kj < ik; kj++)
	    {
	       index = 2*(ik - kj + ki - 1);
	       left_coef[index] = (left_coef[index] + left_coef[index - 2])/2.;
	       index++;
	       left_coef[index] = (left_coef[index] + left_coef[index - 2])/2.;
	    }
	 }
	 
	 for (ki = 1; ki < ik; ki++)
	 {
	    for (kj = 0 ; kj < (ik - ki); kj++)
	    {
	       index = 2*kj;
	       right_coef[index] = (right_coef[index] + 
				    right_coef[index + 2])/2.;
	       index++;
	       right_coef[index] = (right_coef[index] + 
				    right_coef[index + 2])/2.;
	    }
	 }
	 
	 /* Make a recursion. */
	 
	 s1245(left_coef, ik, dim, point, local_tol, (depth + 1), left_weight, 
	       &left_area, &left_moment, stat);
	 if (*stat < 0) goto error;
	 
	 s1245(right_coef, ik, dim, point, local_tol, (depth + 1), right_weight, 
	       &right_area, &right_moment, stat);
	 if (*stat < 0) goto error;
	 
	 weight[0] = left_weight[0] + right_weight[0];
	 weight[1] = left_weight[1] + right_weight[1];
	 *area = left_area + right_area;
	 *moment = (left_moment + right_moment)/4.;
	 
	 if (left_coef != SISL_NULL) freearray(left_coef);
	 if (right_coef != SISL_NULL) freearray(right_coef);
      }
      
   }

   goto out;
     
   /* Error in input. */
  
  err106: 
   *stat = -106;
   s6err("s1245",*stat,pos);
   goto out;     
     
   /* Error in lower level function */
     
  error:
   s6err("s1245", *stat, pos);
   goto out;
     
  out:
   return;
}



