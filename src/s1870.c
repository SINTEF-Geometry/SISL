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

#define S1870

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1870(SISLSurf *ps1, double *pt1, int idim, double aepsge,
	 int *jpt,double **gpar1,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1870(ps1,pt1,idim,aepsge,jpt,gpar1,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double    *pt1;
     int	idim;
     double   aepsge;
     int      *jpt;
     double   **gpar1;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find all intersections between a NURBS surface
*              and a point.
*
*
*
* INPUT      : ps1    - Pointer to the surface.
*              pt1    - coordinates of the point.
*	       idim   - number of coordinates in pt1.
*              aepsge - Geometry resolution.
*              *jstat    - Flag
*                          = 202 : Complicated point-surface intersection
*                                  in 3D. Perform extra interception test.
*
*
*
*
* OUTPUT     : jpt    - Number of single intersection points.
*              gpar1  - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the surface. The points lie
*                       continuous. Intersection curves are stored in wcurve.
*              jcrv   - Number of intersection curves.
*              wcurve - Array containing descriptions of the intersection
*                       curves. The curves are only described by points
*                       in the parameter plane. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*                       If the curves given as input are degnenerate an
*                       intersection point can be returned as an intersection
*                       curve. Use s1327 to decide if an intersection curve
*                       is a point on one of the curves.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* REFERENCES :
*
*-
* CALLS      : sh1870      - Perform the actual intersection.
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 9403.
*
*********************************************************************
*/
{
  int kstat = 0;           /* Local status variable.                       */
  int kpos = 0;            /* Position of error.                           */
  int trackflag = 0;
  int jtrack;
  int *pretop=SISL_NULL;
  SISLTrack **wtrack=SISL_NULL;
  double aepsco = REL_COMP_RES;

  kstat = *jstat;
  sh1870(ps1, pt1, idim, aepsco, aepsge, trackflag, &jtrack, &wtrack,
	 jpt, gpar1, &pretop, jcrv, wcurve, &kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  /*
   * Intersections found.
   * --------------------
   */

  *jstat = kstat;
  goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1870",*jstat,kpos);
    goto out;

  out:
    return;
}
