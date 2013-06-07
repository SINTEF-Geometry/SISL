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
 * $Id: s2543.c,v 1.7 1996-08-02 07:29:05 jka Exp $
 *
 */


#define S2543

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s2543(SISLSurf *surf, int ider, double derive[], double normal[], double *k1,
      double *k2, double d1[], double d2[], int *jstat)
#else
 void s2543(surf, ider, derive, normal, k1, k2, d1, d2, jstat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double *k1;
      double *k2;
      double d1[];
      double d2[];
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute principal curvature (k1,k2) with corresponding 
*                  principal directions (d1,d2) of a surface for 
*                  given values (u,v). This is a lower level routine,
*                  used for evaluation of many T(u,v)'s.
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0.
*                       < 0 : No derivative calculated.
*                       = 0 : Principal curvature calculated.
*                       = 1 : Principal curvature and its first derivative 
*                             calculated.
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  INPUT/OUTPUT :
*
*  OUTPUT       :
*         k1        - Max. principal curvature.
*         k2        - Min. principal curvature.
*         d1        - Max. direction of the principal curvature k1, given 
*                     in local coordiantes (with regard to Xu,Xv).
*                     Dimension = 2.
*         d2        - Min. direction of the principal curvature k2, given 
*                     in local coordiantes (with regard to Xu,Xv).
*                     Dimension = 2.
*        jstat      - Status messages
*                         = 0 : Ok.
*                         < 0 : error.
*
*  METHOD        :  The princpal curvatures -k1 and -k2 are eigenvalues of
*                   dN, thus it turn out that we have to solve a 
*                   eigenvalue/eigenvector problem, see references. 
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*                  Elementary Linear Algebra 5e
*                    (Howard Anton, Wiley, ISBN:0-471-84819-0)
*-
*  CALLS        :  
*
*  LIMITATIONS  :
*                (i) If the surface is degenrated (not regular) at the point
*                    (u,v), it makes now sence to speak about curvature.
*               (ii) If the surface is closed to degenrate, the resuts
*                    can be numerical unstable. 
*               (iv) The dimension of the space in which the surface lies must 
*                    be 1,2 or 3.  The routine return istat < 0.
*  
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.            Date: 1995-1
* REWRITTEN BY : Johannes Kaasa, SINTEF, Oslo, Norway.           Date: 1995-9
*****************************************************************************
*/
{
   double denom;        /* Denominator in a fraction.                 */
   double a, b, c;      /* Coefficients of second degree equation.    */
   double sqrt_arg;     /* Square argument.                           */
   double ratio;        /* Ratio of principal direction.              */
   double length;       /* Parameter length.                          */
   double fundform[6];  /* The coefficients of the fundamental forms.
			   The sequence is: E, F, G, e, f, g.         */
   double transform[4]; /* Transformation matrix.
			   The sequence is a11, a12, a21, a22.        */
   double Su[3];        /* Tangent in first parameter direction.      */
   double Sv[3];        /* Tangent in second parameter direction.     */

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      
      /* Set up the tangents. */
      
      if (surf->idim == 1)
      {
	 Su[0] = 1.;
	 Su[1] = 0.;
	 Su[2] = derive[1];
	 Sv[0] = 0.;
	 Sv[1] = 1.;
	 Sv[2] = derive[2];
      }
      else
      {
	 Su[0] = derive[3];
	 Su[1] = derive[4];
	 Su[2] = derive[5];
	 Sv[0] = derive[6];
	 Sv[1] = derive[7];
	 Sv[2] = derive[8];
      }
      
      /* Calculate the fundamental forms. */
      
      s2513(surf, ider, 2, 1, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      /* Calculate the transformation matrix. */
      
      denom = fundform[0]*fundform[2] - fundform[1]*fundform[1];
      
      transform[0] = (fundform[1]*fundform[4] - fundform[2]*fundform[3])/
	 denom;
      transform[1] = (fundform[1]*fundform[5] - fundform[2]*fundform[4])/
	 denom;
      transform[2] = (fundform[1]*fundform[3] - fundform[0]*fundform[4])/
	 denom;
      transform[3] = (fundform[1]*fundform[4] - fundform[0]*fundform[5])/
	 denom;
      
      /* Calculate the principal curvature. */
      
      a = 1.;
      b = transform[0] + transform[3];
      c = transform[0]*transform[3] - transform[1]*transform[2];
      
      sqrt_arg = b*b - 4.*a*c;
      if (sqrt_arg < REL_PAR_RES)
	 goto war100;
      
      *k1 = (- b + sqrt(sqrt_arg))/(2.*a);
      *k2 = (- b - sqrt(sqrt_arg))/(2.*a);
      
      /* Calculate the principal directions. */
      
      /* Maximal curvature direction. */
      
      if (fabs(transform[0] + *k1) < REL_PAR_RES && 
	  fabs(transform[1]) < REL_PAR_RES)
      {
	 
	 /* Parallel to the u direction. */
	 
	 length = 1./sqrt(Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2]);
	 
	 d1[0] = length;
	 d1[1] = 0.;
      }
      else if (fabs(transform[3] + *k1) < REL_PAR_RES && 
	       fabs(transform[2]) < REL_PAR_RES)
      {
	 
	 /* Parallel to the v direction. */
	 
	 length = 1./sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]);
	 
	 d1[0] = 0.;
	 d1[1] = length;
      }
      else if (fabs(transform[0] + *k1) < fabs(transform[1]))
      {
	 ratio = (transform[0] + *k1)/transform[1];
	 length = 1./sqrt((Su[0] - ratio*Sv[0])*(Su[0] - ratio*Sv[0]) +
			  (Su[1] - ratio*Sv[1])*(Su[1] - ratio*Sv[1]) +
			  (Su[2] - ratio*Sv[2])*(Su[2] - ratio*Sv[2]));
	 
	 d1[0] = length;
	 d1[1] = -ratio*length;
      }
      else
      {
	 ratio = transform[1]/(transform[0] + *k1);
	 length = 1./sqrt((Sv[0] - ratio*Su[0])*(Sv[0] - ratio*Su[0]) +
			  (Sv[1] - ratio*Su[1])*(Sv[1] - ratio*Su[1]) +
			  (Sv[2] - ratio*Su[2])*(Sv[2] - ratio*Su[2]));
	 
	 d1[0] = -ratio*length;
	 d1[1] = length;
      }
      
      /* Minimal curvature direction. */
      
      if (fabs(transform[0] + *k2) < REL_PAR_RES && 
	  fabs(transform[1]) < REL_PAR_RES)
      {
	 
	 /* Parallel to the u direction. */
	 
	 length = 1./sqrt(Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2]);
	 
	 d2[0] = length;
	 d2[1] = 0.;
      }
      else if (fabs(transform[3] + *k2) < REL_PAR_RES && 
	       fabs(transform[2]) < REL_PAR_RES)
      {
	 
	 /* Parallel to the v direction. */
	 
	 length = 1./sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]);
	 
	 d2[0] = 0.;
	 d2[1] = length;	 
      }
      else if (fabs(transform[0] + *k2) < fabs(transform[1]))
      {
	 ratio = (transform[0] + *k2)/transform[1];
	 length = 1./sqrt((Su[0] - ratio*Sv[0])*(Su[0] - ratio*Sv[0]) +
			  (Su[1] - ratio*Sv[1])*(Su[1] - ratio*Sv[1]) +
			  (Su[2] - ratio*Sv[2])*(Su[2] - ratio*Sv[2]));
	 
	 d2[0] = length;
	 d2[1] = -ratio*length;
      }
      else
      {
	 ratio = transform[1]/(transform[0] + *k2);
	 length = 1./sqrt((Sv[0] - ratio*Su[0])*(Sv[0] - ratio*Su[0]) +
			  (Sv[1] - ratio*Su[1])*(Sv[1] - ratio*Su[1]) +
			  (Sv[2] - ratio*Su[2])*(Sv[2] - ratio*Su[2]));
	 
	 d2[0] = -ratio*length;
	 d2[1] = length;
      }
      
   }
   else if (surf->idim == 2) /* 2D surface */
   {
      /* The surface lies in a plane => T(u,v) = 0 */
      
      *k1 = 0.0;
      *k2 = 0.0;
      d1[0] = 1.0;
      d1[1] = 0.0;
      d2[0] = 0.0;
      d2[1] = 1.0;
   }
   else /* When surf->idim != 1,2 or 3 */
   {
      goto err105;
   }

   
   /* Successful computations  */
   
   *jstat = 0;
   goto out;
   
   
   /* The surface does not have principal curvatures. */
   war100:
      if (fabs(sqrt_arg) < REL_PAR_RES)
      {
	 *k1 = - b/(2.*a);
	 *k2 = *k1;
      }
      else
      {
	 *k1 = 0.0;
	 *k2 = 0.0;
      }
   d1[0] = 1.0;
   d1[1] = 0.0;
   d2[0] = 0.0;
   d2[1] = 1.0;
   goto out;
   
   /* Error in input, surf->idim != 1,2 or 3. */
   err105:
      *jstat = -105;
   s6err("s2543",*jstat,0);
   goto out;
   
   error:
      s6err("s2543",*jstat,0);
   goto out;
   
   out:
      
      return;

}
