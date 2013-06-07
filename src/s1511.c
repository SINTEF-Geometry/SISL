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
 * $Id: s1511.c,v 1.3 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1511

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1511(SISLSurf *ps,double qpoint [],double bvec [],int idim,
	   double aepsco,double aepsge,int *jpt,double **gpar,int *jcrv,
	   SISLIntcurve ***wcurve,int *jstat)
#else
void s1511(ps,qpoint,bvec,idim,aepsco,aepsge,jpt, gpar, jcrv, wcurve, jstat)
     SISLSurf *ps;
     double qpoint[];
     double bvec[];
     int idim;
     double aepsco;
     double aepsge;
     int *jpt;
     double **gpar;
     int *jcrv;
     SISLIntcurve ***wcurve;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find the circular silhouette curves and points of a surface.
*              In addition to the points and curves found by this routine,
*              break curves and edge-curves might be silhouette curves.
*
*
*
* INPUT      : ps  -      Pointer to the surface.
*              qpoint -   A point on the spin axis.
*              bvec  -    The circular silhouette axis direction.
*              idim   -   Dimension of the space in which axis lies.
*              aepsco -   Computational resolution.
*              aepsge -   Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single silhouette points.
*              gpar   - Array containing the parameter values of the
*                       single silhouette points in the parameter
*                       plane of the surface. The points lie continuous.
*                       Silhouette curves are stored in wcurve.
*              jcrv   - Number of silhouette curves.
*              wcurve - Array containing descriptions of the silhouette
*                       curves. The curves are only described by points
*                       in the parameter plane. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES : Main routine written by Mike Floater, SI, 1990.
*
* CALLS      : sh1511, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
*
*********************************************************************
*/
{
  int kstat = 0;              /* Local status variable.                      */
  int kpos = 0;               /* Position of error.                          */
  int i;
  int trackflag = 0;
  int jtrack;
  SISLTrack **wtrack=SISL_NULL;
  int jsurf;
  SISLIntsurf **wsurf=SISL_NULL;
  int *pretop=SISL_NULL;

  sh1511(ps,qpoint, bvec, idim, aepsco, aepsge, trackflag,&jtrack,
	 &wtrack,jpt,gpar,&pretop,jcrv,wcurve,&jsurf,&wsurf,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  for(i=0; i<jsurf; i++)
    freeIntsurf(wsurf[i]);
  if(wsurf != SISL_NULL) freearray(wsurf);

  if(jsurf > 0)
    *jstat=10;
  else
    *jstat = kstat;
  goto out;

  /* Error in lower level routine.  */

  error:
    *jstat = kstat;
    s6err ("s1511", *jstat, kpos);
    goto out;

  out:
    return;
}
