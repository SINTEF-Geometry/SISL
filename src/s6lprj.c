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
 * $Id: s6lprj.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */
#define S6LPRJ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
s6lprj(double e1[],double e2[],int idim)
#else
double s6lprj(e1,e2,idim)
     double e1[];
     double e2[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the length of the projection of a vector on 
*              an arbitrary axis.
*
* INPUT      : e1      - The vector to be projected
*              e2      - The arbitrary axis vector
*              dim     - Number of dimensions in the space the vectors lie
*
* OUTPUT     : s6lprj  - The length of the projection vector
*
* METHOD     : The length of the projection vector is calculated as:
*                       __     __
*                       e1 dot e2
*              ||e3|| = ---------
*                       __     __  1/2
*                      (e2 dot e2)
*-
* CALLS      : s6scpr, s6length
*
* WRITTEN BY : Per Evensen, SI, Oslo, Norway. 1991-aug-16
*                                  
*********************************************************************
*/
{
  int kstat;
  double scpr1,scpr2,lproj; 

  scpr1 = s6scpr(e1,e2,idim);
  scpr2 = s6length(e2,idim,&kstat);
  if (kstat == 0) scpr2=0.000001;
  
  lproj = scpr1/scpr2;
  return(lproj);
}

