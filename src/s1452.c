/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
