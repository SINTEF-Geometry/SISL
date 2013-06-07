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
 * $Id: s6crss.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CRSS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6crss(double e1[],double e2[],double e3[])
#else
void s6crss(e1,e2,e3)
     double e1[];
     double e2[];
     double e3[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make the cross product of two 3-D vectors
*
* INPUT      : e1      - First 3-D vector
*              e2      - Second 3-D vector
*
* OUTPUT     : e3      - The vector containing the cross product e1xe2
*
*
* METHOD     : The cross product is calculated by using its definition.
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*
*********************************************************************
*/
{
  e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
  e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
  e3[2] = e1[0]*e2[1] - e1[1]*e2[0];
}
