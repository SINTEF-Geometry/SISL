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
 * $Id: s6ang.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6ANG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6ang(double evec1[],double evec2[],int idim)
#else
double s6ang(evec1,evec2,idim)
     double evec1[];
     double evec2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : s6ang   - Angle in radians between vectors
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*
*********************************************************************
*/                                     
{
  double tscpr,tang,tlength1,tlength2,tcos;
  int    kstat1,kstat2;
  
  tscpr = s6scpr(evec1,evec2,idim);
  
  tlength1 = s6length(evec1,idim,&kstat1);
  tlength2 = s6length(evec2,idim,&kstat2);
  
  if (!kstat1 || !kstat2)
    tang = DZERO;
  else
    {
      tcos = fabs(tscpr/(tlength1*tlength2));
      tcos = MIN((double)1.0,tcos);
      tang = acos(tcos);
    }
  
  return(tang);
}
