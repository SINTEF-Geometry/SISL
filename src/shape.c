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
 * $Id: shape.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SHAPE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      shape(double emid[],double etang[],int idim,int iedge,int *jstat)
#else	 
void shape(emid,etang,idim,iedge,jstat)
     int idim,iedge,*jstat;
     double emid[],etang[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : This routine gives a possibility for the application
*              to adjust the value and the derivatives in the midpoint
*              of the vertex region, i.e. the point in which the region
*              is divided.
*
*
* INPUT      : idim    - Dimension of geometry space.
*              iedge   - Number of edges of the vertex region.
*
*
* INPUT/OUTPUT  : emid    - The value in the midpoint of the region.
*                           Dimension is idim.
*                 etang   - The tangents of the blending surfaces in
*                           the midpoint of the region, along the 
*                           curves which divides the region into 4-sided
*                           blending surfaces. Dimension is iedge*idim.
*                       
*
* OUTPUT     : jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* USE        : The input of the arrays emid and etang may be changed.
*              This will effect the geometry of the blend strongly. Make
*              sure to return sensible midpoint and tangents. The midpoint
*              should lie close to the real midpoint of the region, pulling
*              it close to an edge, may result in blending surfaces with
*              cusps. Also if the tangents is very long, cusps may occur.
*              The numbering of the tangents must be the same as the numbering
*              of the edges, otherwise the G1-continuity will be lost.
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
  int ki;
  double tfac = 1.0;
  
  /*ü printf("Give factor with which to multiply tangent : ");*/
  /*ü scanf("%lf",&tfac);*/
  
  for (ki=0; ki<iedge*idim; ki++) etang[ki] *= tfac;

  *jstat = 0;

  return;
}
