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

#define S1987

#include "sislP.h"                                                 


#if defined(SISLNEEDPROTOTYPES)
void s1987(SISLSurf *ps, double aepsge, int *jgtpi, double **gaxis,
	   double *cang,int *jstat)
#else
void s1987(ps,aepsge,jgtpi,gaxis,cang,jstat)
     SISLSurf *ps;
     double aepsge;
     int   *jgtpi;
     double **gaxis;
     double *cang;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the direction cone of a surface.
*
*
* INPUT      : ps        - Surface to treat.
*              aepsge    - Geometry tolerance.
*
* OUTPUT     : jgtpi     - To mark if the angle of the direction cone is
*                          greater than pi.
*                           0 - The direction cone of the surface
*                               is not greater than pi in any
*                               parameter direction.
*                           1 - The direction cone of the surface
*                               is greater than pi in the first
*                               parameter direction.
*                           2 - The direction cone of the surface is greater
*                               than pi in the second parameter direction.                          
*                          10 - The direction cone of a boundary curve of
*                               the surface is greater than pi in the first
*                               parameter direction.
*                          20 - The direction cone of a boundary curve of
*                               the surface is greater than pi in the second
*                               parameter direction.                      
*              gaxis     - Allocated array containing the coordinates of the
*                          center of the cone. It is only computed if
*                          *jgtpi = 0.
*              cang      - The angle from the center to the boundary of the
*                          cone. It is only computed if *jgtpi = 0.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 9403.
*
*********************************************************************
*/                                     
{
   int kstat = 0;        /* Local status variable.  */
   int kpos = 0;
   int kdim = ps->idim;
   
   /* Allocate scratch for the output array. */
   
   if ((*gaxis = newarray(kdim, DOUBLE)) == SISL_NULL) goto err101;
   
   /* Let s1990 compute the cone. */
   
   s1990(ps, aepsge, &kstat);
   if (kstat < 0) goto error;
   
   /* Copy the resulting cone to output parameters. */
   
   *jgtpi = ps->pdir->igtpi;
   *cang = ps->pdir->aang;
   memcopy(*gaxis, ps->pdir->ecoef, kdim, DOUBLE);
   
  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
  err101: 
     *jstat = -101;
  s6err("s1987",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
  error:
     *jstat = kstat;
  s6err("s1987",*jstat,kpos);
  goto out;
     
  out: 
    return;
}
