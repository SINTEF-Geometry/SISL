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
 * $Id: s1363.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1363

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1363(SISLCurve *pc,double *cmin,double *cmax,int *jstat)
#else
void s1363(pc,cmin,cmax,jstat)
     SISLCurve  *pc;
     double *cmin;
     double *cmax;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline curve
*
* INPUT      : pc     - The B-spline curve.   
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              cmin   - Start of parametrization of curve
*              cmax   - End of parametrization of curve
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
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kstat;          /* Local status variable                           */
  int kpos=0;         /* Position of error                               */
  
  /* Check if curve is correct */
  
  s1707(pc,&kstat);
  if (kstat<0) goto error;
  
  /* Pick parametrization */
  
  *cmin = pc->et[pc->ik - 1];
  *cmax = pc->et[pc->in];
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1363",*jstat,kpos);
  goto out;
 out:
  
  return;
}          
