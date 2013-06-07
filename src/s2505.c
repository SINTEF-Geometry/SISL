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
 * $Id: s2505.c,v 1.10 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2505

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2505(SISLSurf *surf, int der, double derive[], double normal[],
      double *absCurvature, int *jstat)
#else
 void s2505(surf, der, derive, normal, absCurvature, jstat)
      SISLSurf *surf;
      int der;
      double derive[];
      double normal[];
      double *absCurvature;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the absolute curvature K(u,v) of a Surface
*                  for given values (u,v). This is a lower level routine,
*                  used for evaluation of many K(u,v)'s.
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          der      - Not used.
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  INPUT/OUTPUT :
*
*  OUTPUT       :
*    absCurvature   - Absolute curvature of the surface in (u,v) =
*        jstat      - Staus messages
*                         = 2 : Surface is degenrate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is closed to degenrate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : error.
*
*  METHOD        :  The absolute curvature is given by
*
*                      A(x,y) = |k1| + |k2|,
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      A(u,v) = |k1| + |k2|,
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficents of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The rutine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1421() and s1424().
*
*  LIMITATIONS  :
*                (i) If the surface is degenrated (not regular) at the point
*                    (u,v), it makes now sence to speak about the absolute c.
*                    A(u,v). The routine return jstat == 2.
*               (ii) If the surface is closed to degenrate, the absolute c.
*                    A(u,v) can be numerical unstable. The routine return
*                    jstat == 1.
*              (iii) The surface should be C2, since the absolute c. is calculated
*                    from the second derivativs. But since the rutine is using
*                    right derivatives, the absolute c. will be correct (provided
*                    that the surface is not degenerate).
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
			  The sequence is: E, F, G, e, f, g.         */
   double gc;          /* Gaussian curvature.                        */
   double mc;          /* Mean curvature.                            */



   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {   
      s2513(surf, der, 2, 0, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      gc = (fundform[3]*fundform[5]-fundform[4]*fundform[4])
	 /((fundform[0]*fundform[2] - fundform[1]*fundform[1])*
	   (fundform[0]*fundform[2] - fundform[1]*fundform[1]));
      mc = 0.5*(fundform[3]*fundform[2] - 2*fundform[4]*fundform[1] 
			    + fundform[5]*fundform[0])
	 /((fundform[0]*fundform[2] - fundform[1]*fundform[1])
	   *sqrt(fundform[0]*fundform[2] - fundform[1]*fundform[1]));

      *absCurvature = fabs(mc + sqrt(fabs(mc*mc - gc))) +
       fabs(mc - sqrt(fabs(mc*mc - gc)));
   }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => A(u,v) = 0 */

    *absCurvature = 0.0;
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
  s6err("s2505",*jstat,0);
  goto out;
  
error:
  s6err("s2505",*jstat,0);
  goto out;

out:

  return;

}
