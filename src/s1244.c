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
 * $Id: s1244.c,v 1.1 1995-01-03 09:49:21 pfu Exp $
 *
 */


#define S1244

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1244(double knots[], int knot_reg, int first_order, int second_order, 
      int in, int first_index, int second_index, double *integral, int *stat)
#else
void s1244(knots, knot_reg, first_order, second_order, in, first_index, 
	   second_index, integral, stat)
     double knots[];
     int knot_reg;
     int first_order; 
     int second_order;
     int in;
     int first_index;
     int second_index; 
     double *integral; 
     int *stat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To integrate the product of two B-spline basis functions.
*              The order of the functions may be different, but they are
*              defined on the same knot vector.
*              This routine will only perform if
*              (first_order + second_order) < 12 (can easily be extended).
*
* INPUT      : knots    - The common knot vector.
*              knot_reg - k-regularity of the knot vector.
*              first_order  - Order of the first basis.
*              second_order - Order of the second basis.
*              in       - Dimension of the spline space.
*              first_index  - Start knot index of the first basis.
*              second_index - Start knot index of the second basis.
*
*
* OUTPUT     : integral - The resulting integral.
*              stat     - Status messages  
*                         = 0 : OK.
*                         < 0 : Error.
*                         > 0 : Warning.
*
*                                  
* METHOD     : The integration is performed in sequence on each knot 
*              subinterval, by use of Gauss quadrature.
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
   int ki, kj;            /* Index in for loops.                */
   int pos = 0;           /* Position of error.                 */
   int start_knot;        /* Start knot for non-zero integrand. */ 
   int end_knot;          /* End knot for non-zero integrand.   */ 
   int degree;            /* Polynomial degree of the product.  */
   int numb_pnt;          /* Number of evaluation points.       */
   int first_sup;         /* Superfluous knot regularity.       */
   int second_sup;        /* Superfluous knot regularity.       */
   int left;              /* Pointer into knot vector.          */
   double par;            /* Parameter value.                   */
   double sub_integral;   /* Integral on a knot subinterval.    */
   double scale;          /* Scaling factor.                    */
   double nodes[5];       /* Nodes in Gauss integration.        */
   double weights[5];     /* Weights in Gauss integration.      */
   double first_der[12];  /* Evaluation of basis functions.     */
   double second_der[12]; /* Evaluation of basis functions.     */


   /* Initiation. */
   
   first_sup  = knot_reg - first_order;
   second_sup = knot_reg - second_order;
   if (first_sup < 0 || second_sup < 0)
      goto err106;
   
   *integral = 0.0;
   
   start_knot = max(first_index, second_index);
   end_knot = min(first_index + first_order, second_index + second_order);
   if (start_knot >= end_knot)
      goto out;
   
   degree = first_order + second_order - 2;
   numb_pnt = (int) ceil((degree + 1.)/2.);
   numb_pnt = max(numb_pnt, 2);
   if (numb_pnt > 5)
      goto err106;
   
   /* Make a table of Gauss nodes and weights. */
   
   if (numb_pnt == 2)
   {
      nodes[0] = - 0.5773502691;
      nodes[1] = 0.5773502691;
      
      weights[0] = 1.0;
      weights[1] = 1.0;
   }
   else if (numb_pnt == 3)
   {
      nodes[0] = - 0.7745966692;
      nodes[1] = 0.0;
      nodes[2] = 0.7745966692;
      
      weights[0] = 0.5555555555;
      weights[1] = 0.8888888888;
      weights[2] = 0.5555555555;
   } 
   else if (numb_pnt == 4)
   {
      nodes[0] = - 0.8611363115;
      nodes[1] = - 0.3399810435;
      nodes[2] = 0.3399810435;
      nodes[3] = 0.8611363115;
      
      weights[0] = 0.3478548451;
      weights[1] = 0.6521451548;
      weights[2] = 0.6521451548;
      weights[3] = 0.3478548451;
   }
   else
   {
      nodes[0] = - 0.9061798459;
      nodes[1] = - 0.5384693101;
      nodes[2] = 0.0;
      nodes[3] = 0.5384693101;
      nodes[4] = 0.9061798459;
      
      weights[0] = 0.2369268850;
      weights[1] = 0.4786286704;
      weights[2] = 0.5688888888;
      weights[3] = 0.4786286704;
      weights[4] = 0.2369268850;
   }
   
   /* Go through each of the knot subintervals. */

   for (ki = start_knot; ki < end_knot; ki++)
   {
      if ((knots[ki + 1] - knots[ki]) < REL_COMP_RES)
	 continue;
      
      sub_integral = 0.0;
      scale = (knots[ki + 1] - knots[ki])/2.;
      for (kj = 0; kj < numb_pnt; kj++)
      {
	 par = knots[ki] + (nodes[kj] + 1)*scale;
	 
	 left = ki - first_sup;
	 s1220(&knots[first_sup], first_order, in - first_sup,
	       &left, par, 0, first_der, stat);
	 if (*stat < 0) goto error;
	 left = ki - second_sup;
	 s1220(&knots[second_sup], second_order, in - second_sup,
	       &left, par, 0, second_der, stat);
	 if (*stat < 0) goto error;
	 
	 sub_integral += weights[kj]*
	    first_der[first_order - (ki + 1 - first_index)]*
	    second_der[second_order - (ki + 1 - second_index)];
      }
      
      /* Add the contribution from this interval. */
      
      *integral += sub_integral*scale;
   }
   
   goto out;
     
   /* Error, too high order in input. */
  
  err106: 
   *stat = -106;
   s6err("s1244", *stat, pos);
   goto out;     
     
   /* Error in lower level function */
     
  error:
   s6err("s1244", *stat, pos);
   goto out;
     
  out:
     
   return;
}
