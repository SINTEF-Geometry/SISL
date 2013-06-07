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
 * $Id: s1435.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1435

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1435(SISLSurf *ps1,int iedge,SISLCurve **rcedge,double *cpar,int *jstat)
#else
void s1435(ps1,iedge,rcedge,cpar,jstat)
     SISLSurf   *ps1;
     int    iedge;
     SISLCurve  **rcedge;
     double *cpar;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Pick given edge-curve of B-spline surface.
*
*
*
* INPUT      : ps1    - Pointer to surface.
*              iedge  - Number of surface. See figure below.
*
*                           -----------------
*                           !   iedge=2     !
*                           !               !
*                    iedge=3!               !iedge=1
*                           !               !
*                           !               !
*                           !               !
*                           -----------------
*                               iedge=0
*
*
*
* OUTPUT     : rcedge - SISLEdge curve.
*              cpar   - Parameter value of edge in constant direction.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s1436 - Pick curve with constant second parameter.
*              s1437 - Pick curve with constant first parameter.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status parameter.                        */
  int kpos = 0;         /* Position of error.                             */
  double tstart1,tend1; /* Endpoints of parameter interval in first 
			   direction.                                     */
  double tstart2,tend2; /* Endpoints of parameter interval in second 
			   direction.                                     */
  double tpar;          /* Parameter value of curve in constant parameter
			   direction.                                     */
  
  /* Fetch endpoints of parameter intervals.  */
  
  tstart1 = *(ps1->et1 + ps1->ik1 - 1);
  tend1 = *(ps1->et1 + ps1->in1);
  tstart2 = *(ps1->et2 + ps1->ik2 - 1);
  tend2 = *(ps1->et2 + ps1->in2);
  
  /* Find constant parameter of edge. */
  
  if (iedge == 0) tpar = tstart2;
  else if (iedge == 1) tpar = tend1;
  else if (iedge == 2) tpar = tend2;
  else if (iedge == 3) tpar = tstart1;
  
  if (iedge == 0 || iedge == 2)
    {
      
      /* Pick curve with constant second parameter.  */
      
      s1436(ps1,tpar,rcedge,&kstat);
      if (kstat < 0) goto error;
    }
  else if (iedge == 1 || iedge == 3)
    {
      
      /* Pick curve with constant first parameter.  */
      
      s1437(ps1,tpar,rcedge,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* SISLCurve picked.  */
  
  *cpar = tpar;
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1435",*jstat,kpos);
  goto out;
  
 out: return;
}
