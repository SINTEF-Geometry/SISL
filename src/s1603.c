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
 * $Id: s1603.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1603

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1603(SISLSurf *psurf,double *cmin1,double *cmin2,double *cmax1,double *cmax2,int *jstat)
#else
void s1603(psurf,cmin1,cmin2,cmax1,cmax2,jstat)
     SISLSurf   *psurf;
     double *cmin1;
     double *cmin2;
     double *cmax1;
     double *cmax2;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline surface
*
* INPUT      : pc     - The B-spline surface.   
*
* OUTPUT     : cmin1  - Start parameter in first parameter directon. 
*              cmin2  - Start parameter in second parameter directon.
*              cmax1  - End   parameter in first parameter directon. 
*              cmax2    End   parameter in second parameter directon.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error  
*             
* METHOD     : 
*
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6err
*              
*
* WRITTEN BY : Qyvind Hjelle SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kpos=0;              /* Position of error          */
  
  /* Check surf pointer */
  
  if (!psurf) goto err118;
  
  /* Pick parametrization */
  
  *cmin1 = psurf->et1[psurf->ik1-1];
  *cmax1 = psurf->et1[psurf->in1];
  *cmin2 = psurf->et2[psurf->ik2-1];
  *cmax2 = psurf->et2[psurf->in2];
  
  *jstat = 0;
  goto out;
  
  /* Error in input, no B-spline surface given */
  
 err118: 
  *jstat = -118;
  s6err("s1603",*jstat,kpos);
  goto out;
  
 out:
  
  return;
}          
