/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"


#define S6BEZPOWSF

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
s6bezpowsf(double *c1, int order11,int order12, int power, 
		 double *Pascal, double *c1_power)
#else
   void s6bezpowsf(c1, order11, order12, power, Pascal, c1_power)
      double *c1;
      int order11;
      int order12;
      int power;
      double *Pascal;
      double *c1_power;
#endif      
{
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Multiply a given 1D Bezier surface with itself till a
*              specified power.
*
*
*
* INPUT      : c1      - Coefficients of Bezier surface.
*              order11 - Order of surface in first par. dir.
*              order12 - Order of surface in second par. dir.
*              power   - The factor into which the surfaces are
*                        to be multiplied.
*              Pascal  - Factors used in curve multiplication.
*
*
*
* OUTPUT     : c1_power - Coefficients of product surface.
*
*
* REFERENCES : s6multsfs
*
* NOTE       : Factors necessary for the multiplication must be
*              precomputed. If the multiplication routine is to be
*              calles severatl times, it speeds up the computations
*              if this factors are called only once.
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SINTEF, 1993. 
* COMMENTED BY : Vibeke Skytt, SINTEF, 06.94.
*
*********************************************************************
*/
  int p;
  int kgrad11=order11-1;
  int kgrad12=order12-1;
  int kdum1,kdum2;
  double *c1_p;
  double *c1_pm1;
  int pgrad11, pgrad12;

  c1_power[0] = 1;
  for(p=1,pgrad11=0,pgrad12=0,c1_pm1=c1_power,c1_p=c1_power+1;p<=power;
      p++,pgrad11+=kgrad11,pgrad12+=kgrad12,c1_pm1=c1_p,c1_p+=(pgrad11+1)*(pgrad12+1))
    {

      s6multsfs(c1, order11,order12, c1_pm1,pgrad11+1,pgrad12+1,
                 Pascal, c1_p, &kdum1,&kdum2);
            }
  
}
