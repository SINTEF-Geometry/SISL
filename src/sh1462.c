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
 * $Id: sh1462.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH1462

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
typedef void (*fshapeProc)(double [],double [],int,int,int *);
#else
typedef void (*fshapeProc)();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
  sh1462(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
	 double etwist[],double etang[],double eder[],int *jstat)
#else	 
void sh1462(fshape,vboundc,icurv,etwist,etang,eder,jstat)
     fshapeProc    fshape;
     double        etwist[],etang[],eder[];
     SISLCurve     *vboundc[];
     int           icurv,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a three sided vertex region, evaluate the first
*              blending surface in the corner lying in the middle of the
*              vertex region. Compute the tangent vectors in the middle 
*              vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 6, i.e. 2*icurv.
*              icurv   - Number of sides. icurv = 3.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is 3*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Evaluate the Gregory Charrot function in the midpoint of the
*              vertex region. Compute the wanted derivatives.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : sh1466  - Evaluate the Gregory Charrot function.  
*
* WRITTEN BY : Vibeke Skytt, SI, 03.90.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Status variable.  */
  int kder = 2;          /* Number of derivatives to evaluate.  */
  int ki;                /* Counter.  */
  int kdim = 3;          /* Dimension of geometry space.  */
  double tonethird = (double)1.0/(double)3.0;  /* 1/3     */
  double tonesixth = (double)1.0/(double)6.0;  /* 1/6     */
  double sbar[3];        /* Barycentric coordinates of point to evaluate. */
  double sder[18];       /* Value and derivatives of blending. */

  /* Set up the barycentric coordinates of the midpoint of the region. */

  sbar[0] = sbar[1] = sbar[2] = tonethird;
  
  /* Evaluate the Gregory Charrot function at the midpoint. */

  sh1466(vboundc,etwist,kder,sbar,sder,&kstat);
  if (kstat < 0) goto error;
   
  /* Compute tangent vectors.  */

  for (ki=0; ki<kdim; ki++)
    {
      etang[ki] = -sder[kdim+ki]*tonethird + sder[2*kdim+ki]*tonesixth;
      etang[kdim+ki] = sder[kdim+ki]*tonesixth - sder[2*kdim+ki]*tonethird;
      etang[2*kdim+ki] = sder[kdim+ki]*tonesixth + sder[2*kdim+ki]*tonesixth;
    }
  
  /* Application driven routine to alter the midpoint and tangents in the
     midpoint.  */

  fshape(sder,etang,kdim,icurv,&kstat);
  if (kstat < 0) goto error;
  
  /* Copy value and 1. derivatives of first patch.  */

  memcopy(eder,sder,kdim,DOUBLE);
  memcopy(eder+kdim,etang+2*kdim,kdim,DOUBLE);
  memcopy(eder+2*kdim,etang,kdim,DOUBLE);
  
  /* Compute 2. derivatives.  */

  for (ki=0; ki<kdim; ki++)
    {
      eder[3*kdim+ki] = sder[3*kdim+ki]*tonesixth*tonesixth 
	+ (double)2.0*sder[4*kdim+ki]*tonesixth*tonesixth 
	  + sder[5*kdim+ki]*tonesixth*tonesixth;
      eder[4*kdim+ki] = -sder[3*kdim+ki]*tonesixth*tonethird 
	+ sder[4*kdim+ki]*tonesixth*(tonesixth - tonethird)
	  + sder[5*kdim+ki]*tonesixth*tonesixth;
      eder[5*kdim+ki] = sder[3*kdim+ki]*tonethird*tonethird 
	- (double)2.0*sder[4*kdim+ki]*tonethird*tonesixth 
	  + sder[5*kdim+ki]*tonesixth*tonesixth;
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in a lower level function.  */

 error:
  *jstat = kstat;
  goto out;
  
  out :
    return;
}
