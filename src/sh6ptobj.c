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
 * $Id: sh6ptobj.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6PTOBJ

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
sh6ptobj(double *point, SISLObject *obj, double aepsge,
	   double start[], double result[], int *jstat)
#else
void sh6ptobj(point, obj, aepsge, start, result, jstat)
     double *point;
     SISLObject *obj;
     double aepsge;
     double start[];
     double result[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a point and an object to find a closest point
*              or an intersection point.
*
*
* INPUT      : point     - Point in the intersection.
*              obj       - Object in the intersection.
*              aepsge    - Geometry resolution.
*              start[]   - Start parameter value of the iteration on
*                          the object.
*
*
*
* OUTPUT     : result[]  - Parameter value of the object in intersection
*                          point.
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Sep 1991
*              UJK, SI, Sep 1991, include point object
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  double pstart[2];
  double pend[2];
  SISLPoint *sislpt = SISL_NULL;
  double loc_start[2];
  
  /* Test input.  */
  
  if (obj == SISL_NULL) goto err106;
		   
  if ( obj->iobj == SISLSURFACE)
  {
     if ((sislpt = newPoint(point, obj->s1->idim, 0)) == SISL_NULL)
        goto error;

     memcopy(loc_start,start,2,double);
     
     pstart[0] = obj->s1->et1[obj->s1->ik1 - 1];
     pstart[1] = obj->s1->et2[obj->s1->ik2 - 1];
     pend[0]   = obj->s1->et1[obj->s1->in1];
     pend[1]   = obj->s1->et2[obj->s1->in2];
     
     s1773(sislpt, obj->s1, aepsge,
	   pstart, pend, loc_start, result, &kstat);
     if (kstat < 0) goto error;
  }
  else if ( obj->iobj == SISLCURVE)
  {
     if ((sislpt = newPoint(point, obj->c1->idim, 0)) == SISL_NULL)
        goto error;
     
     pstart[0] = obj->c1->et[obj->c1->ik - 1];
     pend[0]   = obj->c1->et[obj->c1->in];
  
     loc_start[0] = start[0];
     s1771(sislpt, obj->c1, aepsge,
	   pstart[0], pend[0], loc_start[0], result, &kstat);
     if (kstat < 0) goto error;
  }
  else if ( obj->iobj == SISLPOINT)
  {
     if(s6dist(point,obj->p1->ecoef,obj->p1->idim) < aepsge)
	kstat = 1;
     else
        kstat = 2;
  }
  else goto err106;
  
  *jstat = kstat;
  goto out;
  
  /* Error in input. */
  
 err106: *jstat = -106;
  s6err("sh6ptobj",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("sh6ptobj",*jstat,kpos);
  goto out;                  
	 
 out:    if (sislpt != SISL_NULL) freePoint(sislpt);
}
