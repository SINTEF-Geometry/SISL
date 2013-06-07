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
 * $Id: s1314.c,v 1.2 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1314

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1314(SISLSurf *ps1,double *epoint,double *enorm,int idim,double aepsco,
	   double aepsge,double amax,SISLIntcurve *pintcr,int icur,
	   int igraph,int *jstat)
#else
void s1314(ps1,epoint,enorm,idim,aepsco,aepsge,amax,pintcr,icur,igraph,jstat)
     SISLSurf     *ps1;
     double   *epoint;
     double   *enorm;
     int      idim;
     double   aepsco;
     double   aepsge;
     double   amax;
     SISLIntcurve *pintcr;
     int      icur;
     int      igraph;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To march an intersection curve desribed by parameter pairs
*              in an intersection curve object, a B-spline surface and
*              a plane.
*
*
* INPUT      : ps1    - Pointer to surface.
*              epoint - SISLPoint in the plane.
*              enorm  - Normal to the plane.
*              idim   - Dimension of the space in which the plane lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              amax   - Maximal allowed step length. If amax <=aepsge
*                       amax is neglected.
*              icur   - Indicator telling if a 3-D curve is to be made 
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:   
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through
*                            s6move and s6line.
*
*
*
* INPUT/OUTPUT:pintcr - The intersection curve. When comming as input
*                       only parameter values it the parameter plane
*                       exist. When comming as output the 3-D geometry
*                       and possibly the curve in the parameter plane
*                       of the surface is added.
*
* OUTPUT:      jstat  - status messages  
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.*                                         = 3      : Marching not succeded
*                         = 0      : ok
*                         < 0      : error
*
*
* METHOD     : An implicit description of the plane is made and then
*              a routine for intersecting implicit represented geometry
*              by a B-spline surface is used.
*
* REFERENCES :
*
*-
* CALLS      : s6err, s1313
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 2. July 1988
* Revised by : Tor Dokken, SI, Oslo, Norway, Mars 1989
*              Corrected use of none normalized use of normal vectors
*
*********************************************************************
*/
{            
  int kpos=0;         /* Position of error                                  */
  int kdeg=1;         /* The degree of the implicit equation of the plane   */
  int kstat;          /* Local status variable                              */
  double simpli[4];   /* Array containing the implicit description of plane */
  double snorm[3];    /* Normalized version of normal vector                */
  
  
  
  if (idim != 3) goto err104;
  
  /* Normalize normal vector */
  
  (void)s6norm(enorm,idim,snorm,&kstat);
  
  simpli[0] = snorm[0];
  simpli[1] = snorm[1];
  simpli[2] = snorm[2];
  simpli[3] = -s6scpr(epoint,snorm,idim);
  
  
  /* Make intersection of implicit surface and B-spline surface */
  
  s1313(ps1,simpli,kdeg,aepsco,aepsge,amax,pintcr,icur,igraph,&kstat);
  if (kstat == -185) goto err185;
  if (kstat < 0) goto error;
  
  *jstat = kstat;
  goto out;
  
  /* Dimension not 3 */
  
 err104: 
  *jstat = -104;                
  s6err("s1314",*jstat,kpos);
  goto out;
  
  /* Couldn't march */
  
 err185:
  *jstat = -185;
  goto out;
  
  /* Error in lower level routine.  */
  
 error: 
  *jstat = kstat;     
  s6err("s1314",*jstat,kpos);
  goto out;
  
 out:
  return;
}                                               
