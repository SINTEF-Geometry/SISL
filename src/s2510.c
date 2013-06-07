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
 * $Id: s2510.c,v 1.4 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S2510

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2510(SISLSurf *surf, int ider, int iside1, int iside2, double parvalue[],
      int *leftknot1,int *leftknot2, double *mehlum, int *jstat)
#else
 void s2510(surf, ider, iside1, iside2, parvalue, leftknot1, leftknot2, mehlum, jstat)
      SISLSurf *surf;
      int    ider;
      int    iside1;
      int    iside2;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double *mehlum;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the third order Mehlum curvature M(u,v) of a 
*                  surface for given values (u,v) = (parvalue[0],parvalue[1]), 
*                  where:
*
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*          iside1   - Indicator telling if the derivatives in the first
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*          iside2   - Indicator telling if the derivatives in the second
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*      parvalue     - Parameter-value at which to evaluate. Dimension of
*                     parvalue is 2.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     first parameter direction where parvalue[0] is found,
*                     that is:
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the routine.
*
*     leftknot2     - Pointer to the interval in the knot vector in the
*                     second parameter direction where parvalue[1] is found,
*                     that is:
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the routine.
*
*  OUTPUT       :
*     mehlum        - Third order Mehlum curvature of the surface in (u,v) =
*                     (parvalue[0],parvalue[1]).
*        jstat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The third order Mehlum curvature is given by
*
*                      M(u,v) = (5G^3P^2 
*                             + (EG + 4F^2)
*                             *(9GQ^2 + 9ES^2 + 6GPS + 6EQT) 
*                             + 5E^3T^2 
*                             - 2F(3EG + 2F^2)(PT + 9QS)
*                             - 30F(G^2PQ + E^2ST))
*                             /(16(EG - F^2)^3).
*
*                  The variables E,F,G,P,Q,S and T
*                  are the coefficients of the first and third fundamental form.
*                  They are given by: 
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. 
*                  P = <N,Xuuu> + 3(a*alpha + b*beta), 
*                  Q = <N,Xuuv> + c*alpha + d*beta + 2a*gamma + 2b*delta, 
*                  S = <N,Xuvv> + 2c*gamma + 2d*delta + a*epsilon + b*mu,
*                  T = <N,Xvvv> + 3(c*epsilon + d*mu), 
*                  where N is normalized, and
*                  a = Ff - Ge,
*                  b = Fe - Ef,
*                  c = Fg - Gf,
*                  d = Ff - Eg,
*                  e, f and g being the second fundamental form coefficients
*                  (e = <N,Xuu>, f = <N,Xuv> and g = <N,Xvv>), and
*                  alpha   = <Xuu,Xu>/||N||^2,
*                  beta    = <Xuu,Xv>/||N||^2,
*                  gamma   = <Xuv,Xu>/||N||^2,
*                  delta   = <Xuv,Xv>/||N||^2,
*                  epsilon = <Xvv,Xu>/||N||^2,
*                  mu      = <Xvv,Xv>/||N||^2.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1422() and s2511().
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Mehlum
*                    M(u,v). The routine returns jstat = 2.
*               (ii) If the surface is closed to degenerate, the Mehlum
*                    M(u,v) can be numerical instable. The routine returns
*                    jstat = 1.
*              (iii) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.
*
*
* WRITTEN BY :    Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-9
******************************************************************************
*/
{
  int kwarn = 0;         /* Local staus variable(warning).                  */
  int kistat = 0;        /* Local staus variable.                           */
  double derive[30];     /* Array containing the computed derivatives.      */
  double normal[3];      /* Array containing the computed normalvektor.     */


  if (ider != 0) goto err178;


  if (surf == SISL_NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

     s1422(surf,3,iside1,iside2,parvalue,leftknot1,leftknot2,derive,normal,
	  &kistat);
    if (kistat > 0) kwarn=kistat;

    if (kistat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (kistat != 2) /* The surface is not degenerate */
    {
      s2511(surf, ider, derive, normal, mehlum, &kistat);

      if (kistat < 0)
	goto error;
    }
    else if (kistat == 2) /* The surface is degenerated. */
    {
      *mehlum = 0.0;
      goto war002;
    }

  }


  /* Successful computations  */

  *jstat = kwarn;
  goto out;


   /* The surface is degenerated at (u,v) */
war002:
  *jstat = 2;
  goto out;

  /* Error. Input (surface) pointer is SISL_NULL. */
err150:
  *jstat = -150;
  s6err("s2510", *jstat, 0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2510",*jstat,0);
  goto out;
  /* Error in lower level routine.  */

error:
  *jstat = kistat;
  s6err("s2510",*jstat,0);
  goto out;


out:

  return;

}
