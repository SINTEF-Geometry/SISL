/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1995 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s2507.c,v 1.10 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2507

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2507(SISLSurf *surf, int ider, double derive[], double normal[],
      double *totalCurvature, int *jstat)
#else
 void s2507(surf, ider, derive, normal, totalCurvature, jstat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double *totalCurvature;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the total curvature T(u,v) of a Surface
*                  for given values (u,v). This is a lower level routine,
*                  used for evaluation of many T(u,v)'s.
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Not used.
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  INPUT/OUTPUT :
*
*  OUTPUT       :
*    totalCurvature - Total curvature of the surface in (u,v).
*        jstat      - Staus messages
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD        :  The total curvature is given by
*
*                      T(x,y) = k1^2 + k2^2,
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      T(u,v) = = k1^2 + k2^2,
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes now sense to speak about T(u,v).
*               (ii) If the surface is closed to degenerate, T(u,v)
*                    can be numerical unstable.
*              (iii) The surface should be C2, since the total c. is calculated
*                    from the second derivatives.
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.  The routine return jstat < 0.
*
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.            Date: 1995-1
* CORRECTED BY :  Johannes Kaasa, SINTEF, Oslo, Norway.          Date: 1995-06
*                 Used absolute valute in square for principal curvature.
* CORRECTED BY :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-8
*                 Calculated the fundamental form coefficients by
*                 calls to s2513.
*****************************************************************************
*/
{
   double fundform[6]; /* The coefficients of the fundamental forms.
			  The sequence is: E, F, G, e, f, g.              */
   double gc;          /* Gaussian curvature.                             */
   double mc;          /* Mean curvature.                                 */
   double k1,k2;       /* Max. and min. principal curvature.              */



   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 2, 0, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      gc = (fundform[3]*fundform[5]-fundform[4]*fundform[4])
	 /((fundform[0]*fundform[2] - fundform[1]*fundform[1])*
	   (fundform[0]*fundform[2] - fundform[1]*fundform[1]));
      mc = 0.5*(fundform[3]*fundform[2] - 2*fundform[4]*fundform[1] 
			    + fundform[5]*fundform[0])
	 /((fundform[0]*fundform[2] - fundform[1]*fundform[1])
	   *sqrt(fundform[0]*fundform[2] - fundform[1]*fundform[1]));

     k1 = mc + sqrt(fabs(mc*mc - gc));
     k2 = mc - sqrt(fabs(mc*mc - gc));

     *totalCurvature = k1*k1 + k2*k2;
   }
  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => T(u,v) = 0 */

    *totalCurvature = 0.0;
  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }

  /* Successful computations  */

  *jstat = 0;
  goto out;


   /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *jstat = -105;
  s6err("s2507", *jstat, 0);
  goto out;
  
error:
  s6err("s2507",*jstat,0);
  goto out;

out:

  return;

}
