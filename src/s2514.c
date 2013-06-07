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
 * $Id: s2514.c,v 1.3 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2514

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2514(SISLSurf *surf, int ider, double derive[], double normal[],
      double gaussian[], int *stat)
#else
 void s2514(surf, ider, derive, normal, gaussian, stat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double gaussian[];
      int *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the numerator and denominator of the Gaussian 
*                  K(u,v) of a Surface for given
*                  values (u,v). This is a lower level routine, used
*                  for evaluation of many K(u,v)'s.
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Only implemented for ider=0 (derivative order).
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  OUTPUT       :
*     gaussian      - The nominator and denominator of the Gaussian for the 
*                     surface in (parvalue[0],parvalue[1]).
*                     Size = 2.
*        stat       - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The Gaussian is given by
*
*                      K(x,y) = (hxx*hyy-hxy^2)/((1+hx^2+hy^2)^2),
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      K(u,v) = (eg-f*f)/(EG-F*F),
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
*                    (u,v), it makes no sense to speak about the Gaussian
*                    K(u,v).
*               (ii) If the surface is closed to degenerate, the Gaussian
*                    K(u,v) can be numerical unstable.
*              (iii) The surface is Cr the Gaussian calculated is C(r-2).
*                    To get the correct behavior use the sided evaluator s1422
*		     instead of s1421.
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3, if not, stat = -105 is returned.
*
*
* WRITTEN BY   :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-8
*****************************************************************************
*/
{
   double fundform[6]; /* The coefficients of the fundamental forms.
			  The sequence is: E, F, G, e, f, g.         */

   if (ider != 0) goto err178;

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 2, 0, derive, normal, fundform, stat);
      if (*stat < 0) goto error;
      
      gaussian[0] = fundform[3]*fundform[5]-fundform[4]*fundform[4];
      gaussian[1] = (fundform[0]*fundform[2] - fundform[1]*fundform[1])*
	   (fundform[0]*fundform[2] - fundform[1]*fundform[1]);
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => K(u,v) = 0 */

    gaussian[0] = 0.0;
    gaussian[1] = 1.0;
  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }

  /* Successful computations  */

  *stat = 0;
  goto out;


   /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *stat = -105;
  s6err("s2514",*stat,0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *stat = -178;
  s6err("s2514",*stat,0);
  goto out;
  
error:
  s6err("s2514",*stat,0);
  goto out;

out:

  return;

}
