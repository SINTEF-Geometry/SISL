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
 * $Id: s2513.c,v 1.8 1995-09-22 13:45:15 jka Exp $
 *
 */


#define S2513

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2513(SISLSurf *surf, int ider, int type, int normalized, double derive[], 
	 double normal[], double fundform[], int *stat)
#else
   void s2513(surf, ider, type, normalized, derive, normal, fundform, stat)
      SISLSurf *surf;
      int ider;
      int type;
      int normalized;
      double derive[];
      double normal[];
      double fundform[];
      int *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the coefficients of the first and second
*                  fundamental form. 
*  INPUT        :
*     surf         - Pointer to the surface to evaluate.
*     ider         - Only implemented for ider=0 (derivative order).
*     type         - Type of fundamental form:
*                    = 1 : first fundamental form,
*                    = 2 : first and second fundamental form,
*                    = 3 : first, second and third fundamental form.
*     normalized   - Flag telling if the second fundamental form
*                    coefficients shall be calculated from normalized or
*                    not normalized normals:
*                    = 0 : not normalized,
*                    = 1 : normalized,
*     derive       - Array containing derivatives from routine s1421().
*                    Size = idim*6. The sequence is:
*                    S, Su, Sv, Suu, Suv, Svv. 
*     normal       - Array containing the normal from routine s1421().
*                    Size = 3, only used for dim = 3.
*
*  OUTPUT       :
*     fundform     - Array for the fundamental form. The sequence is
*                    (E,F,G,e,f,g,P,Q,S,T,Eu,...,Ev,...,Euu,...,Euv,...,Evv,...).
*     stat         - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :
*
*  LIMITATIONS  :
*               (i) The dimension of the space in which the surface lies must
*                   be 1 or 3, if not, stat = -105 is returned.
*
*
* WRITTEN BY   :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-6
*****************************************************************************
*/
{
   double norm_scale; /* Scalation of the surface normal. */
      
   if (ider != 0) goto err178;
   if (type < 1 || type > 3) goto err178;
      
   if (surf->idim == 1) 
   {
	 
      /* 1D surface. */
	 
      if (normalized == 0)
         norm_scale = 1.;
      else
	 norm_scale = sqrt(1. + derive[1]*derive[1] + derive[2]*derive[2]);
	 
      fundform[0] = 1. + derive[1]*derive[1];
      fundform[1] = derive[1]*derive[2];
      fundform[2] = 1. + derive[2]*derive[2];
	 
      if (type > 1)
      {
	 fundform[3] = derive[3]/norm_scale;
	 fundform[4] = derive[4]/norm_scale;
	 fundform[5] = derive[5]/norm_scale;
      }
      if (type > 2)
      {
	 fundform[6] = derive[6]/norm_scale;
	 fundform[7] = derive[7]/norm_scale;
	 fundform[8] = derive[8]/norm_scale;
	 fundform[9] = derive[9]/norm_scale;
      }
   }
      
   else if (surf->idim == 3) 
   {
	 
      /* 3D surface */
	 
      if (normalized == 0)
         norm_scale = 1.;
      else
	 norm_scale = sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
			   normal[2]*normal[2]);
	 
      fundform[0] = s6scpr(&derive[3], &derive[3], 3);
      fundform[1] = s6scpr(&derive[3], &derive[6], 3);
      fundform[2] = s6scpr(&derive[6], &derive[6], 3);
	 
      if (type > 1)
      {
	 fundform[3] = (s6scpr(normal, &derive[9], 3))/norm_scale;
	 fundform[4] = (s6scpr(normal, &derive[12], 3))/norm_scale;
	 fundform[5] = (s6scpr(normal, &derive[15], 3))/norm_scale;
      }
      
      if (type > 2)
      {
	 fundform[6] = (s6scpr(normal, &derive[18], 3))/norm_scale;
	 fundform[7] = (s6scpr(normal, &derive[21], 3))/norm_scale;
	 fundform[8] = (s6scpr(normal, &derive[24], 3))/norm_scale;
	 fundform[9] = (s6scpr(normal, &derive[27], 3))/norm_scale;
      }
   }
   else 
   {
      /* Not 1D or 3D surface */
	 
      goto err105;
   }
      
            
   /* Successful computations  */
      
   *stat = 0;
   goto out;
      
      
   /* Error in input, surf->idim != 1 or 3 */
err105:
   *stat = -105;
   s6err("s2513",*stat,0);
   goto out;
      
   /* Illegal derivative requested. */
err178:
   *stat = -178;
   s6err("s2513",*stat,0);
   goto out;
      
out:	 
   return;
      
}
