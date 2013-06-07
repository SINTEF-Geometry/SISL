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
 * $Id: s2545.c,v 1.8 2001-06-12 11:07:34 jbt Exp $
 *
 */

#define S2545

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
s2545(SISLSurf *surf, int curvature_type, int export_par_val, int pick_subpart,
      double boundary[], int n_u, int n_v, double scale,
      double **garr, int *stat)
#else
void
s2545(surf, curvature_type, export_par_val, pick_subpart, boundary, n_u, n_v,
      scale, garr, stat)
     SISLSurf *surf;
     int curvature_type;
     int export_par_val;
     int pick_subpart;
     double boundary[];
     int n_u;
     int n_v;
     double scale;
     double **garr;
     int *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
* PURPOSE : To compute a set of focal values on an uniform grid
*           in a selected subset of the parameter domain for a
*           NURBS surface. The focal value is a surface position offseted
*           with the surface curvature.
*
*
* INPUT   : surf  	   - The surface to evaluate.
*	    curvature_type - The type of curvature:
*                            0 - Gaussian.
*                            1 - Mean.
*                            2 - Absolute.
*                            3 - Total.
*                            4 - second order Mehlum (curvature).
*                            5 - third order Mehlum (variation of curvature).
*	    export_par_val - Flag telling if the parameter values for each grid
*                            point is to be exported:
*                            0 - False, do not export parameter values,
*                            1 - True, do export parameter values.
*           pick_subpart   - Flag teliing if the grid is to be calculated on a
*                            subpart of the surface:
*                            0 - False, calculate grid on the complete surface,
*                            1 - True, calculate grid on a part of the surface.
*           boundary       - A rectangular subset of the parameter domain.
*              	             [0] - Min value 1. parameter direction.
*                            [1] - Min value 2. parameter direction.
*                            [2] - Max value 1. parameter direction.
*                            [3] - Max value 2. parameter direction.
*                            ONLY USED WHEN pick_subpart = 1. If pick_subpart
*                            = 0, the parameter area of surf is given out here.
*           n_u            - Number of segments in 1. parameter direction.
*           n_v            - Number of segments in 2. parameter direction.
*           scale          - Scaling factor.
*
*
*
* OUTPUT  : garr  	   - Array containing the computed values on the grid.
*		             The allocation is done internally and the dimension
*			     is  (dim+2)*(n_u+1)*(n_v+1) if export_par_val is
*                            true, and dim*(n_u+1)*(n_v+1) if export_par_val is
*                            false.
*                            Each gridpoint consist of dim+2 values
*                            (Ui,Vj,x(Ui,Vj),...) or only the
*                            focal points (x(Ui,Vj),....).
*			     The sequence is running first in the
*                            1. parameter direction.
**
*           stat           - Status message.
*                               > 0      : Warning.
*                               = 0      : Ok.
*                               < 0      : Error.
*
*
* METHOD  :
*
*
* CALLS   :s2540(),s1421(),s6norm()
*
* WRITTEN BY :  Johannes Kaasa,  SINTEF, Oslo, Norway, Sep. 1995.
*
*********************************************************************
*/
{
   int ki, kj, kl;        /* Indices in for loops.            */
   int idx1, idx2;        /* Array indices.                   */
   int incr;              /* Increment in the 3D array.       */
   int leftknot1;         /* Pointer into knot array.         */
   int leftknot2;         /* Pointer into knot array.         */
   double length;         /* Length of normal.                */
   double par[2];         /* Parameter values.                */
   double derive[9];      /* Surface evaluation.              */
   double normal[3];      /* Surface normal.                  */
   double Nnormal[3];     /* Normalized normal.               */
   double *offset = SISL_NULL; /* Offset distances.                */


   /* Generate the scalar valued curvature grid with parameter values. */

   s2540(surf, curvature_type, 1, pick_subpart, boundary, n_u, n_v,
	 &offset, stat);
   if (*stat < 0) goto error;

   /* Allocate the output. */

   incr = (export_par_val? (surf->idim + 2) : surf->idim);

   if (((*garr) = newarray(incr*(n_u + 1)*(n_v + 1), double)) == SISL_NULL)
      goto err101;

   /* Generate the focal points. */

   idx1 = 0;
   idx2 = 0;
   for (ki = 0; ki < (n_v + 1); ki++)
   {
      par[1] = offset[idx1 + 1];

      for (kj = 0; kj < (n_u + 1); kj++)
      {
	 par[0] = offset[idx1];

	 if (export_par_val)
	 {
	    (*garr)[idx2] = par[0];
	    idx2++;
	    (*garr)[idx2] = par[1];
	    idx2++;
	 }

	 /* Calculate the normal. */

	 s1421(surf, 1, par, &leftknot1, &leftknot2, derive, normal, stat);
	 if (*stat < 0) goto error;

	 /* Calculate the point. */

	 if (surf->idim == 1)
	 {
	    (*garr)[idx2] = derive[0] + scale*offset[idx1 + 2];
	    idx2 += 1;
	 }
	 else if (surf->idim == 2)
	 {
	    (*garr)[idx2] = scale*offset[idx1 + 2];
	    idx2 += 1;
	 }
	 else if (surf->idim == 3)
	 {
	    length = s6norm(normal, 3, Nnormal, stat);
	    if (*stat < 0) goto error;
	    for (kl = 0; kl < 3; kl++)
	       (*garr)[idx2 + kl] = derive[kl]
		  + scale*offset[idx1 + 2]*Nnormal[kl];
	    idx2 += 3;
	 }

	 idx1 += 3;
      }
   }


  *stat = 0;
  goto out;

  /* ___________________________________________________________________ */
  /*                         ERROR EXITS                                 */
  /* ___________________________________________________________________ */

  /* Error in space allocation. */
err101:
  *stat = -101;
  s6err("s2545", *stat, 0);
  goto out;

  /* Error in lower level routine. */
error:
  s6err("s2545", *stat, 0);
  goto out;

  /* ___________________________________________________________________ */
  /*                         THE ONE AND ONLY EXIT                       */
  /* ___________________________________________________________________ */

out:
   if (offset) freearray(offset);
  return;

}
