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
 * $Id: s1243.c,v 1.3 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1243

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1243(SISLCurve *pcurve, double point[], int dim, double epsge,
      double weight[], double *area, double *moment, int *stat)
#else
void s1243(pcurve, point, dim, epsge, weight, area, moment, stat)
     SISLCurve  *pcurve;
     double     point[];
     int        dim;
     double     epsge;
     double     weight[];
     double     *area;
     double     *moment;
     int        *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate the weight point and rotational momentum of
*              an area between a 2D curve and a 2D point. The area is
*              also calculated.
*              When the curve is rotating counter-clockwise around the
*              point, the area contribution is positive.
*              When the curve is rotating clockwise around the point,
*              the area contribution is negative.
*              OBSERVE: FOR CALCULATION OF AREA ONLY, USE s1241().
*
* INPUT      : pcurve - The 2D curve.
*              point  - The reference point.
*              dim    - Dimension of geometry (must be 2).
*              epsge  - Absolute geometrical tolerance.
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
   int ki, kj;                    /* Index in for loops.              */
   int pos = 0;                   /* Position of error.               */
   int numb_seg;                  /* Number of segments.              */
   int depth;                     /* Depth of recursion.              */
   double cord_lenght;            /* Cord lenght.                     */
   double offset_tol;             /* Offset tolerance.                */
   double max;                    /* Maximum step length.             */
   double local_tol;              /* Local tolerance.                 */
   double prev_area;              /* Previous area.                   */
   double bez_area;               /* Bezier area.                     */
   double bez_moment;             /* Bezier moment.                   */
   double dummy[3];               /* Dummy array.                     */
   double bez_weight[2];          /* Bezier weight point.             */
   SISLCurve* non_rat_crv = SISL_NULL; /* Local non-rational curve.        */
   SISLCurve* local_crv = SISL_NULL;   /* Local curve (do not deallocate). */
   SISLCurve* non_per_crv = SISL_NULL; /* Local non-periodic curve.        */
   SISLCurve* bez_crv = SISL_NULL;     /* Bezier convertion of curve.      */


   /* Check input. */

   if (pcurve->idim != 2 || dim != 2)
      goto err106;
   if (pcurve->ik < 1)
      goto err106;
   if ((epsge - 0.0) < REL_COMP_RES)
      goto err106;

   /* Make sure the curve is k-regular, and non-rational. */

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

   /* Convert the curve to Bezier segments. */

   s1730(local_crv, &bez_crv, stat);
   if (*stat < 0) goto error;

   /* Find number of segments. */

   numb_seg = bez_crv->in/bez_crv->ik;

   /* Make a loop with regard to tolerance. */

   local_tol = max(0.1, 10.1*epsge);
   prev_area = 0.0;
   *area = -1.0;

   while (fabs(prev_area - *area) > epsge && local_tol > epsge)
   {

      local_tol *= 0.1;
      prev_area = *area;

      /* Do calculations for each segment. */

      weight[0] = 0.0;
      weight[1] = 0.0;
      *area     = 0.0;
      *moment   = 0.0;

      for (ki = 0; ki < numb_seg; ki++)
      {
	 depth = 1;
	 s1245(&(bez_crv->ecoef[ki*bez_crv->ik*bez_crv->idim]), bez_crv->ik,
	       bez_crv->idim, point, local_tol, depth, bez_weight, &bez_area,
	       &bez_moment, stat);
	 if (*stat < 0) goto error;

	 weight[0] += bez_weight[0];
	 weight[1] += bez_weight[1];
	 *area += bez_area;
	 *moment += bez_moment;
      }

      if (fabs(*area) > REL_COMP_RES)
      {
         weight[0] /= *area;
         weight[1] /= *area;
      }
   }

   goto out;

   /* Error in input. */

  err106:
   *stat = -106;
   s6err("s1243",*stat,pos);
   goto out;

   /* Error in lower level function */

  error:
   s6err("s1243", *stat, pos);
   goto out;

  out:
   if (non_rat_crv != SISL_NULL) freeCurve(non_rat_crv);
   if (non_per_crv != SISL_NULL) freeCurve(non_per_crv);
   if (bez_crv != SISL_NULL) freeCurve(bez_crv);

   return;
}
