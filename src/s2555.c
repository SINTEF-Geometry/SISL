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
 * $Id:
 *
 */


#define S2555

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2555( double     derive[],
       double    *torsion,
       int       *jstat )
#else
void s2555( derive, torsion, jstat )
     double      derive[];
     double     *torsion;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the torsion of a curve based on the
*              derivatives in "derive". We assume that the
*              curve lies in IR^3.
*
*
* INPUT      :
*            derive   - Double array of dimension [4*3]
*                       containing the position and derivative vectors.
*                       These vectors are stored in the following order:
*                       First the position vector, then the tangent vector,
*                       then the 3 components of the second derivative
*                       vector, and so on.
*
*
* OUTPUT     : torsion - The torsion value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : See formula in the book of Do Carmo
*              (Differential Geometry of Curves and Surfaces).
*
* REFERENCES : Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*
* CALLS       : s6length(), s6scpr()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable       */
  double crpr[3], a;    /* Temp. variable.             */


  /* Evaluate torsion */

  s6crss( derive + 3, derive + 6, crpr );

  a = s6length( crpr, 3, &kstat );

  if ( a != 0.0 )
    {
      *torsion = ( s6scpr( derive + 9, crpr, 3 ) ) / ( a*a ) ;
    }
  else
    {
      *torsion = 0.0;

      goto war002;
    }



 /* Successful computations.  */

  *jstat = 0;
  goto out;


  /* Torsion undefined, since the curvature = 0.0.  */
war002:
  *jstat = 2;
  goto out;


 out:

 return;

}
