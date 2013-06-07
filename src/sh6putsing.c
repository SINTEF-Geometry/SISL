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
 * $Id: sh6putsing.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6PUTSING
#include "sislP.h"                            

#if defined (SISLNEEDPROTOTYPES)
void
      sh6putsing( SISLIntpt *psource, SISLIntpt *pdest)
#else
	 
	 void sh6putsing(psource,pdest)
	    
	    SISLIntpt *psource,*pdest;
#endif

/*
*********************************************************************
*                                                                   
* PURPOSE    : To set some approximative values in a singular point
*              based on symmetry.
*
*
*
* INPUT      : psource - The previous point on the curves.
*              pdest   - The point that is to be fit in.
*                        Posisiton must be ok.
* OUTPUT     : pdest   - Geometric values, ie tangent is changed.

*
* METHOD     : For each curve we mirror the tangent in psource
*              around the difference vector pdest-psource.
*
* REFERENCES :
*
*-
* CALLS      : s6scpr, s6diff, s6norm
*
* WRITTEN BY : UJK
*
*********************************************************************
*/
{
  int kdim,kstat;                
  double alfa;
  double diffv[3];
  double delta[3];
  int ki;
  
  kdim = 3;
  s6diff(pdest->geo_track_3d,psource->geo_track_3d,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_3d+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_3d[kdim+ki] = alfa*delta[ki] - psource->geo_track_3d[kdim+ki];

  pdest->geo_track_3d[9] = -(double) 1.0;
  
  kdim = 2;
  s6diff(pdest->geo_track_2d_1,psource->geo_track_2d_1,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_2d_1+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_2d_1[kdim+ki] = alfa*delta[ki] - psource->geo_track_2d_1[kdim+ki];
  
  pdest->geo_track_2d_1[6] = -(double) 1.0;

  kdim = 2;
  s6diff(pdest->geo_track_2d_2,psource->geo_track_2d_2,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_2d_2+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_2d_2[kdim+ki] = alfa*delta[ki] - psource->geo_track_2d_2[kdim+ki];

  pdest->geo_track_2d_1[6] = -(double) 1.0;

}
