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
 * $Id: s1319.c,v 1.3 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1319

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1319(SISLSurf *ps1,double *eview,int idim,double aepsco,double aepsge,
	   double amax,SISLIntcurve *pintcr,int icur,int igraph,int *jstat)
#else
void s1319(ps1,eview,idim,aepsco,aepsge,amax,pintcr,icur,igraph,jstat)
     SISLSurf     *ps1;
     double   *eview;
     int      idim;
     double   aepsco;
     double   aepsge;
     double   amax;
     SISLIntcurve *pintcr;
     int      icur;
     int      igraph;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To march the silhouette curve described by an intersection
*              curve object, a B-spline surface and a view direction.
*
*
* INPUT      : ps1    - Pointer to surface.
*              eview  - View direction
*              idim   - Dimension of the space in which the view direction
*                       lies.
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
* INPUT/OUTPUT:pintcr - The intersection curve. When comming as input
*                       only parameter values it the parameter plane
*                       exist. When comming as output the 3-D geometry
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
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 19 . November 1988
*
*********************************************************************
*/
{
  int kpos=0;         /* Position of error                                  */
  int kdeg=1003;      /* The degree of the implicit equation of the plane   */
  int kstat;          /* Local status variable                              */
  double simpli[4];   /* Array containing the implicit description of plane */
  double snorm[3];    /* Normalized version of normal vector                */



  if (idim != 3) goto err104;

  /* Normalize normal vector */

  (void)s6norm(eview,idim,snorm,&kstat);

  simpli[0] = snorm[0];
  simpli[1] = snorm[1];
  simpli[2] = snorm[2];

  /* Make intersection of implicit surface and B-spline surface */

  s1313(ps1,simpli,kdeg,aepsco,aepsge,amax,pintcr,icur,igraph,&kstat);
  if (kstat == -185) goto err185;
  if (kstat < 0) goto error;

  *jstat = kstat;
  goto out;

  /* Dimension not 3 */

 err104:
  *jstat = -104;
  s6err("s1319",*jstat,kpos);
  goto out;

  /* Couldn't march */

 err185:
  *jstat = -185;
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1319",*jstat,kpos);
  goto out;

 out:
  return;
}
