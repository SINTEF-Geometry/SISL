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
 * $Id: s6dist.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6dist(double epoint1[],double epoint2[],int idim)
#else
double s6dist(epoint1,epoint2,idim)
     double epoint1[];
     double epoint2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between the points epoint1 and
*              epoint2.
*
*
*
* INPUT      : epoint1 - First point in distance calculation.
*              epoint2 - Second point in distance calculation.
*              idim    - Dimension of the space in which the points lie.
*
*
*
* OUTPUT     : s6dist  - Distance between the points.
*
*
* METHOD     : Compute lenght of the vector epoint1-epoint2.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/                                     
{
  register double *s1,*s2,*s3; /* Pointers used to travers epoint1 and epoint2
				  arrays.                                      */
  register double tdist=DZERO; /* Distance between the points.                 */
  
  for (s1=epoint1,s2=epoint2,s3=epoint1+idim; s1<s3; s1++,s2++)
    tdist += (*s1 - *s2)*(*s1 - *s2);
  
  return(sqrt(tdist));
}
