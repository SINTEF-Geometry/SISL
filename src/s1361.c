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
 * $Id: s1361.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */
#define S1361

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1361(double epnt1[],double epnt2[],int idim,
	   double gmidd[],double gmtang[],int *jstat)
#else
void s1361(epnt1,epnt2,idim,gmidd,gmtang,jstat)
     double epnt1[];
     double epnt2[];
     int    idim;
     double gmidd[];
     double gmtang[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To determint if the shape described by the input, points,
*              tangents, curvatures and radius of curvature when used
*              for making an Hermit segment, have a shape close to
*              a circular arc.
*
* INPUT      : epnt1   - Start point with tangent, curvature and radius
*                        of curvature.
*              epnt2   - End point with tangent, curvature and radius
*                        of curvature.
*              idim    - The dimension of the space the point lie in
*              gmidd   - The middle point of the Bezier segement
*              gmtang  - The tangent at the middle of the Bezier segment
*
*
* OUTPUT     : jstat   - Status variable
*                         0 - Shape not acceptabel
*                         1 - Shape acceptable
*
* METHOD     : 
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 18. Oct 1988
*              UJK, October 1990, Included test for negative curvature. 
*********************************************************************
*/                        
{
  double tang1,tang2;    /* Tangent lengths */
  double tscal1,tscal2;  /* The cosine of an angle         */
  double ta1,ta2;        /* An angle                       */
  double tlength;        /* The length of a vector         */
  double tdiff;          /* Difference between two numbers */
  double tdist;          /* Distance between points        */
  double tv2,tv3;        /* Vertex compnents               */
  int    ki;             /* Variable in loop               */
  int    kstat;          /* Local status variable          */
  
  /* Find angle between tangents of epnt1 and epnt2, we assume that
     the tangents are normalized */
  
  tscal1  = s6scpr(epnt1+idim,epnt2+idim,idim);
  
  if (tscal1 >= DZERO)
    tscal1  = MIN((double)1.0,tscal1);
  else
    tscal1  = MAX((double)-1.0,tscal1);
  
  ta1 = acos(tscal1);
  
  if (fabs(ta1) < ANGULAR_TOLERANCE) ta1 = DZERO;
  
  
  /* Make distance between epnt1 and epnt2 */
  
  tdist = s6dist(epnt1,epnt2,idim);
  
  /* Make tangent lengths for start and end points */
  
  if (DNEQUAL(ta1,DZERO))
    {
      /*  Make tangents based on radius of curvature */
      
      tang1 = s1325(epnt1[3*idim],ta1);
      tang2 = s1325(epnt2[3*idim],ta1);
    }
  
  /* Make sure that the tangent does not explode due to numeric errors, and
     make a controlled tangent when the radius is zero or almost zero  */
  
  /* UJK, October 90, must include the case negative curvature */
  if (DEQUAL(ta1,DZERO) || tang1 > tdist || epnt1[3*idim] <= DZERO)
    tang1 = tdist/(double)3.0;
  if (DEQUAL(ta1,DZERO) || tang2 > tdist || epnt2[3*idim] <= DZERO) 
    tang2 = tdist/(double)3.0;
  
  
  /* We now know the Bezier polygon of the Hermit curve. Make angles
     between line 1 and 2 and between line 2 and 3. Make length of line 3
     */
  
  tscal1 = DZERO;
  tscal2 = DZERO;
  tlength = DZERO;
  
  for (ki=0;ki<idim;ki++)
    {
      /*  Make difference between second and third vertex, and accumulte
	  scalar products between polygon lines */
      tv2 = epnt1[ki] + tang1*epnt1[ki+idim]; 
      tv3 = epnt2[ki] - tang2*epnt2[ki+idim];
      tdiff = tv3 - tv2 ;
      tlength += tdiff*tdiff;
      tscal1  += tdiff*epnt1[ki+idim]; 
      tscal2  += tdiff*epnt2[ki+idim];
      
      /*  Make midpoint and tangent at midpoint */
      
      gmidd[ki]  = (epnt1[ki] + (double)3.0*(tv2+tv3) + epnt2[ki])/(double)8.0;
      gmtang[ki] = (epnt2[ki] + tv3 - tv2 - epnt1[ki])/(double)8.0;
      
    }                                                      
  tlength = sqrt(tlength);
  if (DEQUAL(tlength,DZERO)) tlength = (double)1.0;
  
  tscal1 = tscal1/tlength;
  tscal2 = tscal2/tlength;

  if (tscal1 >= DZERO)
    tscal1  = MIN((double)1.0,tscal1);
  else
    tscal1  = MAX((double)-1.0,tscal1);

  if (tscal2 >= DZERO)
    tscal2  = MIN((double)1.0,tscal2);
  else
    tscal2  = MAX((double)-1.0,tscal2);
  
  ta1    = acos(tscal1);
  ta2    = acos(tscal2);
  
  /* Normalize tangent at midpoint */
  
  (void)s6norm(gmtang,idim,gmtang,&kstat);
  
  
  /* Make total angular change of polygon */
  
  ta1 = fabs(ta1) + fabs(ta2);
  
  /* If total angular change is greater than PI/3 or the length is greater
     than 0.45 x the distance don't accept. The last condition make sure
     that the middle span in the polygon is less that what we get when
     we have a 90 circular arc. The first condition makes sure that the
     polygon direction is not oscillating too much*/
  
  if (ta1 > (double)1.0 || tlength > (double)0.45*tdist)
    *jstat = 0;
  else
    *jstat = 1;
}

