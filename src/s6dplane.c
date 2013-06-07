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
 * $Id: s6dplane.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DPLANE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
  s6dplane(double eq1[],double eq2[],double eq3[],double epoint[],
	   int idim,int *jstat)
#else
double s6dplane(eq1,eq2,eq3,epoint,idim,jstat)
   double eq1[];
   double eq2[];
   double eq3[];
   double epoint[];
   int    idim;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between a plane given by three
*              points in the plane and a fourth point.
*
*
*
* INPUT      : eq1    - Point in the plane
*              eq2    - Point in the plane
*              eq3    - Point in the plane
*              epoint - Point. Dimension is idim.
*              idim   - Dimension of geometry space.
*
*
*
* OUTPUT     : s6dplane - Distance between plane and point.
*              jstat   - status messages  
*                        > 0      : warning
*                        = 0      : ok
*                        < 0      : error
*
*
* METHOD     : 
*              
*
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   -  Scalar product between two vectors.
*              s6diff   -  Difference vector between two vectors.  
*              s6norm   -  Normalize vector.
*              s6crss   -  Cross product between two vectors.
*              s6dist   -  Distance between points.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
   int kstat = 0;         /* Local status varaible.           */
   double tdist;          /* Distance between point and line. */
   double snorm[3];       /* Normal vector to the plane.      */
   double sdiff1[3];      /* Difference vector between points in the plane. */
   double sdiff2[3];      /* Difference vector between points in the plane. */
   double sdiff3[3];      /* Difference vector.               */
   
   /* Test dimension.     */
   
   if (idim != 3) goto err104;
   
   /* Compute difference vectors.  */
   
   s6diff(eq2,eq1,idim,sdiff1);
   s6diff(eq3,eq1,idim,sdiff2);
   s6diff(epoint,eq1,idim,sdiff3);
   
   /* Compute normalized plane normal.  */
   
   s6crss(sdiff1,sdiff2,snorm);
   (void)s6norm(snorm,idim,snorm,&kstat);
   
   /* Compute distance to closest point in plane. */
   
   if (kstat)
      tdist = fabs(s6scpr(sdiff3,snorm,idim));
   else 
      tdist = s6dist(eq1,epoint,idim);   /* Normal of zero length.  */

   /* Set status.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in input, dimension not equal to 3.  */
   
   err104 : *jstat = -104;
   goto out;
   
   out :
      return tdist;
 }
