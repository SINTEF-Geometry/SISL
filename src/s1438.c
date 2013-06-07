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
 * $Id: s1438.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1438

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1438(SISLCurve  *pc,int iedge,SISLPoint **rpedge,double *cpar,int *jstat)
#else
void s1438(pc,iedge,rpedge,cpar,jstat)
     SISLCurve  *pc;
     int    iedge;
     SISLPoint  **rpedge;
     double *cpar;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Pick given edge-point of B-spline curve.
*
*
*
* INPUT      : pc     - Pointer to curve.
*              iedge  - Number of point. See figure below.
*
*                           
*                    iedge=0-----------------iedge=1
*
*
*
* OUTPUT     : rpedge - SISLEdge point.
*              cpar   - Parameter value of edge.
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 05.89.
*
*********************************************************************
*/                                     
{
  int kpos = 0;         /* Position of error.                             */
  
  if (iedge == 0)
    {
      *cpar = pc->et[pc->ik - 1];
      
      if (((*rpedge) = newPoint(pc->ecoef,pc->idim,1)) == SISL_NULL)
	goto err101;
    }
  else if (iedge == 1)
    {
      *cpar = pc->et[pc->in];
      
      if (((*rpedge)=newPoint(pc->ecoef+pc->idim*(pc->in-1),pc->idim,1))==SISL_NULL)
	goto err101;
    }
  else goto err141;
  
  /* SISLPoint picked.  */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1438",*jstat,kpos);
  goto out;
  
  /* Error in number of edges.  */
  
  err141 : *jstat = -141;
  s6err("s1438",*jstat,kpos);
  goto out;
  
 out: ;
}
