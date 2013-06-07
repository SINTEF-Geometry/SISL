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
 * $Id: s1774.c,v 1.2 2001-03-19 15:58:53 afr Exp $
 *
 */
#define S1774

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1774(SISLCurve *crv, double point[], int dim, double epsge,
	   double start, double end, double guess, double *clpar, int *stat)
#else
void s1774(crv, point, dim, epsge, start, end, guess, clpar, stat)
     SISLCurve  *crv;
     double point[];
     int dim;
     double epsge;
     double start;
     double end;
     double guess;
     double *clpar;
     int    *stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameter is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : crv     - The curve in the closest point problem.
*              point   - The point in the closest point problem.
*              dim     - Dimension of the geometry.
*              epsge   - Geometrical resolution.
*              start   - Curve parameter giving the start of the search
*                        interval.
*              end     - Curve parameter giving the end of the search
*                        interval.
*              guess   - Curve guess parameter for the closest point
*                        iteration.
*
*
*
* OUTPUT     : clpar   - Resulting curve parameter from the iteration.
*              stat    - status messages
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Aug 1995
*
*********************************************************************
*/
{
   int kpos = 0;             /* Error indicator. */
   SISLPoint* ppoint = SISL_NULL; /* SISL point.      */
   
   /* Generate a SISL point. */
   
   ppoint = newPoint(point, dim, 0);
   
   /* Call s1771. */
   
   s1771(ppoint, crv, epsge, start, end, guess, clpar, stat);
   if (*stat < 0) goto error;

   goto out;

   /* Error in lower level routine.  */

   error :
   s6err("s1774", *stat, kpos);
   goto out;

 out:    
    if (ppoint != SISL_NULL) freePoint(ppoint);
}

