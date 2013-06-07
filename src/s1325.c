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
 * $Id: s1325.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1325

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s1325(double aradiu,double angle)
#else
double s1325(aradiu,angle)
     double aradiu;
     double angle;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create the tangent length for interpolating a
*              circular arc with an almost equi-oscillating Hermit qubic
*
* INPUT      : aradiu  - The radius of the circular arc
*              angle   - The opening angle of the circular arc
*
* OUTPUT     : s1325   - The proposed tangent length
*
* METHOD     : A second degree equation giving the tanget length is
*              solved
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 30. June 1988
*                                  
*********************************************************************
*/
{
  double tcos,tsin;          /* Dummy variables                     */
  double ta,tb,tc,tl;        /* Dummy variables                     */
  double tconst = (double)1.85530139760811990992528773586425;
                             /* Constant used in the calculation    */
  
  
  
  tcos = cos(angle);
  tsin = sin(angle);
  
  /*  Calculate length of tangents
   *   tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 */
  
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)
           * tcos - (double)0.4*tconst - (double)1.0;
  tl     = aradiu*(-tb+sqrt(tb*tb-4*ta*tc))/((double)2.0*ta);
  
  return(tl);
}
