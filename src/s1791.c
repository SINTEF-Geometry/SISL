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
 * $Id: s1791.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1791

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s1791(double et[],int ik,int in)
#else
int s1791(et,ik,in)
     double et[];
     int    ik;
     int    in;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if it is possible to insert new internal nots
*              any further.
*
*
*
* INPUT      : et     - The knot vector.
*              ik     - The order of the curve/surface.
*              in     - The number of the basic functions.
*
*
*
* OUTPUT     : s1791  - Result of the test.
*                       = 0 : It is not possible to insert new knots
*                             any further.
*                       = 1 : It is possible to insert new knots.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/                                     
{
  register double tstart= et[ik - 1];
  register double tend  = et[in];
  register double tmid  = (tstart+tend)*(double)0.5;
  
  /* Check if it is possible to divide the parameter interval.  */
  
  if (DEQUAL(tmid,tstart) || DEQUAL(tmid,tend)) 
    return  0;
  else 
    return  1;
}





