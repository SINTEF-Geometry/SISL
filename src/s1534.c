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

#define S1534

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1534(double points[],double der10[], double der01[],
	   double der11[],int im1,int im2,int idim,int ipar,
	   int con1,int con2,int con3,int con4, int order1,
	   int order2,int iopen1,int iopen2,SISLSurf **rsurf,int *jstat)
#else
void s1534(points,der10,der01,der11,im1,im2,idim,ipar,con1,con2,con3,
	   con4,order1, order2,iopen1,iopen2,rsurf,jstat)
     double points[];
     double der10[];
     double der01[];
     double der11[];
     int im1;
     int im2;
     int idim;
     int ipar;
     int con1;
     int con2;
     int con3;
     int con4;
     int order1;
     int order2;
     int iopen1;
     int iopen2;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute a B-spline tensor surface interpolating a set
*          of points.
*
* INPUT:
*          points - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          der10  - Array of dimension idim*im1*im2 containing the first
*                   derivatives in the first parameter direction.
*
*          der01  - Array of dimension idim*im1*im2 containing the first
*                   derivatives in the second parameter direction.
*
*          der11  - Array of dimension idim*im1*im2 containing the cross
*                   derivatives (the twists).
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
*          ipar   - Flag showing the desired parametrization to be used:
*                   = 1: Mean accumulated cord-length parameterization.
*                   = 2: Uniform parametrization.
*
*
*                          ^ Second par. direction 
*                          |     
*                          |    (2.)
*                          |-----------|
*                          |           |
*                     (3.) |           | (4.) 
*                          |           |
*                          |           |
*                          |-----------|-> First par. direction
*                               (1.)
*
*          con1      - Additional condition along edge 1:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con2      - Additional condition along edge 2:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con3      - Additional condition along edge 3:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con4      - Additional condition along edge 4:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          order1    - Order of surface in first parameter direction.
*
*          order2    - Order of surface in second parameter direction.
*
*          iopen1    - Open/closed parameter in first parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed surface.
*                      = -1 : closed, periodic surface.
*
*          iopen2    - Open/closed parameter in second parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed surface.
*                      = -1 : closed, periodic surface.
*
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     First, a suitable parametrization is calculated according
*     to indicator variable ipar.
*     Then the interpolation is accomplished by using a one dimensional
*     routine for spline interpolation called several times. 
*     First, the datapoints
*     are considered to be idim*im1 dimentional and so on...
*
*
* REFERENCES :
*
* CALLS      : s1528, s1535.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, May 1993.
*
*********************************************************************
*/                                                               
{
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *par1=SISL_NULL;    /* Transposed positions (in rpos)              */
  double *par2=SISL_NULL;    /* Transposed derivatives (in rder)            */
  
  
  /* Check input */        
  
  if (ipar < 1 || ipar > 3) goto err102;
  
  /* Generate parametrizations */
  
  s1528(idim, im1, im2, points, ipar, iopen1, iopen2, &par1, &par2, &kstat);
  if(kstat < 0) goto error;

  /* Interpolation */

  s1535(points,der10,der01,der11,im1,im2,idim,par1,par2,
	con1,con2,con3,con4,order1, order2, iopen1, iopen2, rsurf,&kstat);
  if (kstat < 0) goto error;

  /* Success */
  
  *jstat = 0;
  goto out;  
  
  /* Error in input data. */

  err102: *jstat = -102;
    s6err("s1534",*jstat,kpos);
    goto out;
    
  /* Error in lower level routine. */

  error:  *jstat =kstat;
    s6err("s1534",*jstat,kpos);
    goto out;
  
  out:
    if(par1 != SISL_NULL) freearray(par1);
    if(par2 != SISL_NULL) freearray(par2);
    return;
}
