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
 * $Id: s6affdist.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6AFFDIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
     s6affdist(double e1[],double e2[],double emat[],int idim)
#else
double s6affdist(e1,e2,emat,idim)
   double e1[];
   double e2[];
   double emat[];
   int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between two points using an
*              affine metric described by the matrix emat.
*
*
*
* INPUT      : e1     - First point.
*              e2     - Second point.
*              emat   - Matrix of affine metric.
*              idim   - Dimension of geometry space.
*              
*
* OUTPUT     : s6affdist - Distance between two points.
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 91-03.
*
*********************************************************************
*/
{
   int ki,kj;              /* Counters.  */
   double tdist = DZERO;   /* Distance.  */
   
   for (ki=0; ki<idim; ki++)
      for (kj=0; kj<idim; kj++)
	 tdist += emat[ki*idim+kj]*(e1[ki]-e2[ki])*(e1[kj]-e2[kj]);
   
   tdist = sqrt(idim*tdist);
   
   return tdist;
}
