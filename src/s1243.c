/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
