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
 * $Id: crvlintang.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define CRV_LIN_TANG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
   crv_lin_tang(SISLCurve *pc1, double point[], double normal[],
	   double ang_tol, double guess_par, double *iter_par,
	   int *jstat)
#else
     void crv_lin_tang(pc1, point, normal, ang_tol, guess_par,
		       iter_par, jstat)
     SISLCurve   *pc1;
     double point[];
     double normal[];
     double ang_tol;
     double guess_par;
     double *iter_par;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration to find a tangent between the curve
*              pc1 and a line.
*
*
* INPUT      : pc1       - Pointer to the first curve.
*              point     - Original point on the line.
*              normal    - Normal to the line.
*              ang_tol   - Angular tolerance (in radians).
*              guess_par - Guess parameter values in pc1.
*
*
*
* OUTPUT     : iter_par  - Tangential parameter values in pc1.
*              jstat   - status messages  
*                                < 0   : error.
*
*
* METHOD     : Use po_crv_tang.c to find the tangent from point to pc1,
*              and check the tangent direction against normal.
*              Guess_par and iter_par must not be separated by a tangential
*              discontinuity.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SI, March 1992.
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.           */
  int kpos = 0;             /* Position of error.               */
  
  int kder = 0;             /* Evaluates only position.         */
  int kleft = 0;            /* Knot interval pointer.           */
  double iter_pnt[2];       /* Resulting point on the curve.    */ 
  double diffvec[2];        /* Vector between point and curve.  */
  double tangent[2];        /* Tangent along the line.          */
  int kdim = 2;             /* 2 dimensional.                   */
  double iter_ang;          /* Angular deviation.               */
  
  /* Test input.  */
  
  if (pc1->idim != 2) goto err106;
  
  /* Find the tangent from point to pc1. */
  
  po_crv_tang(pc1, point, ang_tol, guess_par, iter_par, &kstat);
  if (kstat < 0) goto error;
  
  /* Check the result. */

  s1221(pc1, kder, *iter_par, &kleft, iter_pnt, &kstat);
  if (kstat < 0) goto error;
  diffvec[0] = iter_pnt[0] - point[0];
  diffvec[1] = iter_pnt[1] - point[1];
  tangent[0] = -normal[1];
  tangent[1] = normal[0];
  iter_ang = s6ang(diffvec, tangent, kdim);
  if (iter_ang < ang_tol)
    *jstat = 1;
  else
    *jstat = 2;

  goto out;
  
  /* Error in input. Conflicting dimensions.  */
  
 err106: *jstat = -106;
  s6err("crv_lin_tang",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("crv_lin_tang",*jstat,kpos);
  goto out;                  
  
 out: return;
}


