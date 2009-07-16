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
 * $Id: s2503.c,v 1.10 2001-06-12 11:07:33 jbt Exp $
 *
 */


#define S2503

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s2503(SISLSurf *surf, int ider, double derive[], double normal[],
	   double *meancurvature, int *jstat)
#else
 void s2503(surf, ider, derive, normal, meancurvature, jstat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double *meancurvature;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the mean curvature H(u,v) of a surface for given
*                  values (u,v). This is a lower level routine, used
*                  for evaluation of many H(u,v)'s.
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Only implemented for ider=0 (derivative order).
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*  OUTPUT       :
*    meancurvature  - Mean curvature of the surface in (u,v) =
*        jstat      - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The mean curvature is given by
*
*                      H(x,y) = 0.5((1+hx^2)hyy-2hxhyhxy+(1+hy^2)hxx)/
*                                      ((1+hx^2+hy^2)^3/2),
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      H(u,v) = 0.5(eG-2fF+gE)/(EG-F*F),
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The routine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*  CALLS        :
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the mean curvature
*                    H(u,v).
*               (ii) If the surface is closed to degenerate, the mean curvature
*                    H(u,v) can be numerical unstable.
*              (iii) If the surface is Cr the curvature calculated is C(r-2).
*                    To get the correct behavior use the sided evaluator s1422
*		     instead of s1421.
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3, if not, jstat = -105 is returned.
*
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.             Date: 1995-1
* CORRECTED BY :  Ulf J Krystad, SINTEF, Oslo, Norway.            Date: 1995-1
*                 Removed knot navigators + some clean up.
* CORRECTED BY :  Johannes Kaasa, SINTEF, Oslo, Norway.           Date: 1995-8
*                 Calculated the fundamental form coefficients by
*                 calls to s2513.
******************************************************************************
*/
{
   double fundform[6]; /* The coefficients of the fundamental forms.
			  The sequence is: E, F, G, e, f, g.         */

   if (ider != 0) goto err178;

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 2, 0, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      *meancurvature = 0.5*(fundform[3]*fundform[2] - 2*fundform[4]*fundform[1] 
			    + fundform[5]*fundform[0])
	 /((fundform[0]*fundform[2] - fundform[1]*fundform[1])
	   *sqrt(fundform[0]*fundform[2] - fundform[1]*fundform[1]));
   }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => H(u,v) = 0 */

    *meancurvature = 0.0;
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
  s6err("s2503", *jstat, 0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2503",*jstat,0);
  goto out;
  
error:
  s6err("s2503",*jstat,0);
  goto out;

out:

  return;
}
