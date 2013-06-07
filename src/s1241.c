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
 * $Id: s1241.c,v 1.4 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1241

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1241(SISLCurve *pcurve, double point[], int dim, double epsge,
      double *area, int *stat)
#else
void s1241(pcurve, point, dim, epsge, area, stat)
     SISLCurve  *pcurve;
     double     point[];
     int        dim;
     double     epsge;
     double     *area;
     int        *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate the area between a 2D curve and a 2D point.
*              When the curve is rotating counter-clockwise around the
*              point, the area contribution is positive.
*              When the curve is rotating clockwise around the point,
*              the area contribution is negative.
*              If the curve is closed or periodic, the area calculated
*              is independent of where the point is situated.
*              The area is calculated exactly for B-spline curves, for
*              NURBS the result is an approximation.
*              This routine will only perform if the order of the curve is
*              less than 7 (can easily be extended).
*
* INPUT      : pcurve - The 2D curve.
*              point  - The reference point.
*              dim    - Dimension of geometry (must be 2).
*              epsge  - Absolute geometrical tolerance.
*
*
* OUTPUT     : area   - Calculated length.
*              stat   - Status messages
*                       = 0 : OK.
*                       < 0 : Error.
*                       > 0 : Warning.
*
*
* METHOD     : The area calculation can be obtained by a sum of integrations:
*              area = Sum(i)Sum(j)(c_x(i)*c_y(j)*(Integral of)(N(i)*dN(j)/dt))
*                     - (c_x(in-1)*c_y(in-1) - c_x(0)*c_y(0))/2,
*              where c_x and c_y are the x and y components of the curve
*              coefficients with regard to the input point, and N(i) is the
*              i'th basis function.
*              The integration is performed by use of Gauss quadrature.
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
   int ki, kj;                    /* Index in for loops.              */
   int pos = 0;                   /* Position of error.               */
   int start;                     /* Start index.                     */
   int stop;                      /* Stop index.                      */
   double max;                    /* Maximum step length.             */
   double cord_lenght;            /* Cord lenght.                     */
   double offset_tol;             /* Offset tolerance.                */
   double knot_span;              /* Span of knot subinterval.        */
   double first_intgr;            /* First integration.               */
   double second_intgr;           /* Second integration.              */
   double dummy[3];               /* Dummy array.                     */
   double* x_comp = SISL_NULL;         /* x-component of the coordinates.  */
   double* y_comp = SISL_NULL;         /* y-component of the coordinates.  */
   SISLCurve* non_rat_crv = SISL_NULL; /* Local non-rational curve.        */
   SISLCurve* non_per_crv = SISL_NULL; /* Local non-periodic curve.        */
   SISLCurve* local_crv = SISL_NULL;   /* Local curve (do not deallocate). */



   /* Check input. */

   if (pcurve->idim != 2 || dim != 2)
      goto err106;

   /* Make local coordinates. Make sure the curve is k-regular, and
      non-rational.                                                   */

   if (pcurve->ikind == 2 || pcurve->ikind == 4)
   {
      /* The offset tolerance is made a bit ad hoc. */

      cord_lenght = 0.0;
      for (ki = 1, kj = 2; ki < pcurve->in; ki++, kj += 2)
	 cord_lenght += sqrt((pcurve->ecoef[kj] - pcurve->ecoef[kj - 2])*
			     (pcurve->ecoef[kj] - pcurve->ecoef[kj - 2]) +
			     (pcurve->ecoef[kj + 1] - pcurve->ecoef[kj - 1])*
			     (pcurve->ecoef[kj + 1] - pcurve->ecoef[kj - 1]));
      if (cord_lenght < REL_COMP_RES)
	 goto err106;
      offset_tol = epsge/cord_lenght;

      max = 0.0;
      s1360(pcurve, 0.0, offset_tol, dummy, max, dim, &non_rat_crv, stat);
      if (*stat < 0) goto error;
      local_crv = non_rat_crv;
   }
   else
      local_crv = pcurve;

   if (local_crv->cuopen == SISL_CRV_PERIODIC )
   {
      s1712 (local_crv, local_crv->et[local_crv->ik - 1],
	     local_crv->et[local_crv->in], &non_per_crv, stat);
      if (*stat < 0) goto error;

      local_crv = non_per_crv;
   }

   x_comp = newarray(local_crv->in, double);
   y_comp = newarray(local_crv->in, double);

   for (ki = 0, kj = 0; ki < local_crv->in; ki++, kj += 2)
   {
      x_comp[ki] = local_crv->ecoef[kj]     - point[0];
      y_comp[ki] = local_crv->ecoef[kj + 1] - point[1];
   }

   /* Add up the contributions */

   *area = 0.;

   /* The outer loop runs along x. */

   for (ki = 0; ki < local_crv->in; ki++)
   {
      start = max(ki - (local_crv->ik - 1), 0);
      stop  = min(ki + (local_crv->ik - 1) + 1, local_crv->in);

      /* The inner loop runs along y. */

      for (kj = start; kj < stop; kj++)
      {

	 /* Do the Gauss quadrature. */

	 knot_span = local_crv->et[kj + local_crv->ik - 1] -
	    local_crv->et[kj];
	 if (kj > 0 && knot_span > REL_COMP_RES)
	 {
	    s1244(local_crv->et, local_crv->ik, local_crv->ik,
		  local_crv->ik - 1, local_crv->in, ki, kj,
		  &first_intgr, stat);
	    if (*stat < 0) goto error;

	    first_intgr *= (local_crv->ik - 1)/knot_span;
	 }
	 else
	    first_intgr = 0.0;

	 knot_span = local_crv->et[kj + local_crv->ik] -
	    local_crv->et[kj + 1];
	 if (kj < (local_crv->in - 1) && knot_span > REL_COMP_RES)
	 {
	    s1244(local_crv->et, local_crv->ik, local_crv->ik,
		  local_crv->ik - 1, local_crv->in, ki, kj + 1,
		  &second_intgr, stat);
	    if (*stat < 0) goto error;

	    second_intgr *= (local_crv->ik - 1)/knot_span;
	 }
	 else
	    second_intgr = 0.0;

	 *area += x_comp[ki]*y_comp[kj]*(first_intgr - second_intgr);
      }
   }

   /* Add the residue. */

   *area += (x_comp[0]*y_comp[0] - x_comp[local_crv->in - 1]*
	     y_comp[local_crv->in - 1])/2.;

   goto out;

   /* Error in input. */

  err106:
   *stat = -106;
   s6err("s1241",*stat,pos);
   goto out;

   /* Error in lower level function */

  error:
   s6err("s1241", *stat, pos);
   goto out;

  out:
   if (non_rat_crv != SISL_NULL) freeCurve(non_rat_crv);
   if (non_per_crv != SISL_NULL) freeCurve(non_per_crv);
   if (x_comp != SISL_NULL) freearray(x_comp);
   if (y_comp != SISL_NULL) freearray(y_comp);

   return;
}
