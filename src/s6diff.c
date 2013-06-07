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
 * $Id: s6diff.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DIFF

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6diff(double e1[],double e2[],int idim,double e3[])
#else
void s6diff(e1,e2,idim,e3)
     double e1[];
     double e2[];
     double e3[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the difference between two vectors
*
* INPUT      : e1      - The first vector
*              e2      - The second vector
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : 
*              e3      - The difference of e1 and e2
*
* METHOD     : The difference is calculated by vetor substraction
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-june-1988
*                                  
*********************************************************************
*/
{
  int ki;
  for (ki=0;ki<idim;ki++)
    e3[ki] = e1[ki] - e2[ki];
  return;
}
