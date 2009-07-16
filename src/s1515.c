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
 * $Id: s1515.c,v 1.6 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1515

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1515 (SISLSurf * ps1, double qpoint[], double bvec[], int idim, double aepsco, double aepsge,
       double amax, SISLIntcurve * pintcr, int icur, int igraph, int *jstat)
#else
void
s1515 (ps1, qpoint, bvec, idim, aepsco, aepsge, amax, pintcr, icur, igraph, jstat)
     SISLSurf *ps1;
     double qpoint[];
     double bvec[];
     int idim;
     double aepsco;
     double aepsge;
     double amax;
     SISLIntcurve *pintcr;
     int icur;
     int igraph;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To march the circular silhouette curve described by an intersection
*              curve object, a B-spline surface, point Q and direction B
*              i.e. solution of  f(u,v) = N(u,v) x (P(u,v) - Q) . B
*
*
* INPUT      : ps1    - Pointer to surface.
*              qpoint  - Point Q for circular silhouette
*              bvec  - Direction B for circular silhouette
*              idim   - Dimension of the space in which Q lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              amax   - Maximal allowed step length. If amax <=aepsge
*                       amax is neglected.
*              icur   - Indicator telling if a 3-D curve is to be made
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through
*                            s6move and s6line.
*
*
*
* INPUT/OUTPUT:pintcr - The intersection curve. When coming in as input
*                       only parameter values in the parameter plane
*                       exist. When coming as output the 3-D geometry
*                       and possibly the curve in the parameter plane
*                       of the surface is added.
*
* OUTPUT:      jstat  - status messages
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.
*                         = 0      : ok
*                         < 0      : error
*                         = -185   : No points produced on intersection curve.
*
*
* METHOD     : An implicit description of the problem is made and then
*              a routine for intersecting implicit represented geometry
*              by a B-spline surface is used.
*
* REFERENCES :
*
*-
* CALLS      : s6err, s1313
* WRITTEN BY : Mike Floater, SI, Oslo, Norway, 31 . January 1991
*                a modification of s1319
*
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error                                  */
  int kdeg = 1005;		/* The degree of the implicit equation                */
  /* of the circular silhouette                         */
  int kstat;			/* Local status variable                              */
  double simpli[6];		/* Array containing the implicit description          */
  /* of the circular silhouette                         */
  double snorm[3];		/* Normalized version of direction vector bvec        */



  if (idim != 3)
    goto err104;

  /* Normalize bvec direction vector */

  (void) s6norm (bvec, idim, snorm, &kstat);


  simpli[0] = qpoint[0];
  simpli[1] = qpoint[1];
  simpli[2] = qpoint[2];

  simpli[3] = snorm[0];
  simpli[4] = snorm[1];
  simpli[5] = snorm[2];

  /* Make intersection of implicit surface and B-spline surface */

  s1313 (ps1, simpli, kdeg, aepsco, aepsge, amax, pintcr, icur, igraph, &kstat);
  if (kstat == -185)
    goto err185;
  if (kstat < 0)
    goto error;

  *jstat = kstat;
  goto out;

  /* Dimension not 3 */

err104:
  *jstat = -104;
  s6err ("s1515", *jstat, kpos);
  goto out;

  /* Couldn't march */
  /* Only degenerate or singular guide points */
err185:
  *jstat = -185;
  goto out;

  /* Error in lower level routine.  */

error:
  *jstat = kstat;
  s6err ("s1515", *jstat, kpos);
  goto out;

out:
  return;
}
