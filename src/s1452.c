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

#define S1452

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1452(SISLSurf *ps, double aepsge, double aoffset, SISLSurf **rs,
	   int *jstat)
#else
void s1452(ps,aepsge,aoffset,rs,jstat)
     SISLSurf   *ps;
     double aepsge;
     double aoffset;
     SISLSurf   **rs;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To create a B-spline approximating the offset surface of
*              a B-spline surface.
*
*
*
* INPUT      : ps     - The input B-spline surface.
*              aepsge - Maximal deviation allowed between true offset surface
*                       and the approximated offset surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*
* OUTPUT     :
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer the approximated offset surface
*
* METHOD     : The 4 edge curves of the surface are extracted. Offset curves
*              of these 4 edge curves are approximated and a common
*              basis for the two pairs of opposite offset curves is calculated.
*              Vertices are recomputed.
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*-
* CALLS      : s1365     - This routine is performing the actual computation
*                          The existence of the current routine is for
*                          historical reasons.
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 0394.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int kdim = ps->idim;  /* Dimension of geometry space.                  */
  double tmax = (double)0;  /* Dummy.                                    */

  /* Compute offset surface. */

   s1365(ps, aoffset, aepsge, tmax, kdim, rs, &kstat);
   if (kstat < 0) goto error;

  /* Surface approximated. */

  *jstat = kstat;
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1452",*jstat,kpos);
  goto out;

  out: return;
}
