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


#define S1951

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1951(double etau[], double ecoef[], int in, int ik, int idim, 
	 int ilend, int irend, int incont, double efac[])
#else
void s1951(etau, ecoef, in, ik, idim, ilend, irend, incont, efac)
   double etau[];
   double ecoef[];
   int    in;
   int    ik;
   int    idim;
   int    ilend;
   int    irend;
   int    incont;
   double efac[];
#endif     
/*
*********************************************************************
* 
* PURPOSE    : Multiply the coefficients by dtau(-1/2) and express 
*              the incont last coefficients as a weighted sum
*              of the incont first coeffecients. The weights are given
*              in efac.
* 
* 
* INPUT      : ecoef  - Coefficients of spline curve.
*              in     - Number of coefficients.p in
*              idim   - Dimension of geometry space.
*              incont - Number of continuity conditions, i.e. number of
*                       coefficients at the end to be expressed by
*                       coefficients at the start.
*              efac   - Factors.
*              
*
* 
* OUTPUT     : ecoef  - Coefficients of spline curve.
*             
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt,  SINTEF Oslo, 01.95.
*
*********************************************************************
*/
{
   int ki, kj, kr;   /* Counters.  */
   int kstop;
   double tw;

   /* Multiply the part of ec pointed to by kstart and kstop by the
      corresponding parts of the square matrix dtau(-1/2).   */
   
   for (kstop=in-MAX(incont,irend), ki=ilend; ki<kstop; ki++)
     {
	tw = sqrt((double)ik/(etau[ki+ik] - etau[ki]));
	for (kj=0; kj<idim; kj++)
	  ecoef[ki*idim+kj] *= tw;
     }  
   
   /* Express the incont last coefficients by the incont first ones given
      the factors stored in efac. See s1947. */
   
   for (ki=0; ki<incont; ki++)
   {
      for (kr=0; kr<idim; kr++)
      {
	 ecoef[(in-ki-1)*idim+kr] = DZERO;
	 for (kj=0; kj<=ki; kj++)
	    ecoef[(in-ki-1)*idim+kr] += ecoef[kj*idim+kr]*efac[ki*incont+kj];
      }
   }
   
}
   
