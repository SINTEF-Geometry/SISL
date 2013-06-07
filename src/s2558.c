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


#define S2558

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2558( double     derive[],
       int        idim,
       double    *VoC,
       int       *jstat )
#else
void s2558( derive, idim, VoC, jstat )
     double      derive[];
     int         idim;
     double     *VoC;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the Variation of Curvature (VoC) of a curve
*              based on the derivatives in "derive". We assume that
*              the curve lies in R, R^2 or R^3.
*
*
* INPUT      :
*            derive   - Double array of dimension [4*idim]
*                       containing the position and derivative vectors.
*                       These vectors are stored in the following order:
*                       First the position vector, then the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*
*
* OUTPUT     : VoC    - The Variation of Curvature (VoC) value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : See formula in given in the reference.
*
*
* REFERENCES :   Westgaard G., Heimann J.:
*                Analysis and Visualization of Geometric Curve Properties,
*                in M. D{\ae}hlen, T. Lyche and L. L. Schumaker (Eds.)
*                Mathematical Methods in CAGD IV, 1998.
*
*
*
* CALLS       : s6length(), s6scpr(), s6crss()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kstat = 0;                         /* Local status variable       */

  double a, atimes2, atimes4, a1dota2;   /* Temp. variable.             */
  double crpr12[3];          		 /* Temp. variable.             */
  double crpr13[3];   			 /* Temp. variable.             */
  double crpr13dotcrpr12, normcrpr12;    /* Temp. variable.             */
  double D[15];   			 /* Derivatives on canonical
					    R^3 form                    */


  /* Put the derivatives on canonical R^3 form */

  if ( idim == 1 )
  {
      D[0]  = 0.0;   D[1] = derive[0];  D[2]  = 0.0;
      D[3]  = 1.0;   D[4] = derive[1];  D[5]  = 0.0;
      D[6]  = 0.0;   D[7] = derive[2];  D[8]  = 0.0;
      D[9]  = 0.0;  D[10] = derive[3];  D[11] = 0.0;
      D[12] = 0.0;  D[13] = derive[4];  D[14] = 0.0;
  }
  else if ( idim == 2 )
  {
      D[0]  = derive[0];   D[1] = derive[1];   D[2]  = 0.0;
      D[3]  = derive[2];   D[4] = derive[3];   D[5]  = 0.0;
      D[6]  = derive[4];   D[7] = derive[5];   D[8]  = 0.0;
      D[9]  = derive[6];  D[10] = derive[7];   D[11] = 0.0;
      D[12] = derive[8];  D[13] = derive[9];   D[14] = 0.0;
  }
  else
  {
      D[0]  = derive[0];   D[1] = derive[1];   D[2]  = derive[2];
      D[3]  = derive[3];   D[4] = derive[4];   D[5]  = derive[5];
      D[6]  = derive[6];   D[7] = derive[7];   D[8]  = derive[8];
      D[9]  = derive[9];  D[10] = derive[10];  D[11] = derive[11];
      D[12] = derive[12]; D[13] = derive[13];  D[14] = derive[14];
  }



  /* Evaluate VoC */

  a = s6length( D + 3, 3, &kstat );

  atimes2 = a*a;
  atimes4 = atimes2*atimes2;

  a1dota2 = s6scpr( D + 3, D + 6, 3 );

  s6crss( D + 3, D + 6, crpr12 );
  s6crss( D + 3, D + 9, crpr13 );


  crpr13dotcrpr12 = s6scpr( crpr13, crpr12, 3 );

  normcrpr12 = s6length( crpr12, 3, &kstat );


  if ( a != 0.0 && normcrpr12 != 0.0 )
    {
      *VoC = ( ( crpr13dotcrpr12 / normcrpr12 ) -
                 ( 3.0 * a1dota2 * normcrpr12 / atimes2 ) ) / atimes4 ;
    }
  else
    {
      *VoC = 0.0;
    }




 /* Successful computations.  */

  *jstat = 0;
  goto out;


 out:

 return;

}
