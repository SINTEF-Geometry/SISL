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
 * $Id: s2509.c,v 1.5 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2509

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2509(SISLSurf *surf, int ider, double derive[], double normal[],
      double *mehlum, int *jstat)
#else
 void s2509(surf, ider, derive, normal, mehlum, jstat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double *mehlum;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the second Mehlum curvature M(u,v) of a surface 
*                  for given values (u,v). This is a lower level routine, used
*                  for evaluation of many M(u,v)'s.
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Only implemented for ider=0 (derivative order).
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  OUTPUT       :
*       mehlum      - The second order Mehlum curvature of the surface.
*        jstat      - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The second order Mehlum curvature is given by
*
*                      M(u,v) = 3(eG-2fF+gE)^2/8(EG-F*F)^3 - 
*                               (eg-f*f)/2(EG-F*F)^2.
*
*                  The variables E,F,G,e,f and g
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
*                    (u,v), it makes no sense to speak about the Mehlum
*                    M(u,v).
*               (ii) If the surface is closed to degenerate, the Mehlum
*                    M(u,v) can be numerical unstable.
*              (iii) The dimension of the space in which the surface lies must
*                    be 1,2 or 3, if not, jstat = -105 is returned.
*
*
* WRITTEN BY   :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-8
*****************************************************************************
*/
{
   double fundform[6]; /* The coefficients of the fundamental forms.
			  The sequence is: E, F, G, e, f, g.         */
   double numerator;   /* Value of a numerator.                      */
   double denominator; /* Value of a denominator.                    */

   if (ider != 0) goto err178;

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 2, 0, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      numerator = fundform[3]*fundform[2] - 2*fundform[4]*fundform[1]
	 + fundform[5]*fundform[0];
      denominator = fundform[0]*fundform[2] - fundform[1]*fundform[1];
      
      *mehlum = (3*numerator*numerator)/(8*denominator*denominator*denominator)
	 - (fundform[3]*fundform[5]-fundform[4]*fundform[4])/
	 (2*denominator*denominator);
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => K(u,v) = 0 */

    *mehlum = 0.0;
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
  s6err("s2509",*jstat,0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2509",*jstat,0);
  goto out;
  
error:
  s6err("s2509",*jstat,0);
  goto out;

out:

  return;

}
