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
 * $Id: sh6ptouchp.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6PUTTOUCH
#include "sislP.h"                            

#if defined (SISLNEEDPROTOTYPES)
void
      sh6puttouch( SISLIntpt *psource, SISLIntpt *pdest, int seq)
#else
	 
	 void sh6puttouch(psource,pdest, seq)
	    
	    SISLIntpt *psource,*pdest;
	    int seq;
#endif

/*
*********************************************************************
*                                                                   
* PURPOSE    : To set right direction in a touch point (Singular, but no branch).
*              The invariant to hold is that the direction on the curve 
*              in the parameter space is parallel to the delta vector between
*              the to points.
*
*
*
* INPUT      : psource - The previous point on the curves.
*              pdest   - The point that is to be fit in.
*                        Posisiton must be ok.
*              seq     - The sequncing of the points
*                        +1 - psource comes before pdest.
*                        -1 - psource comes after pdest.
* OUTPUT     : pdest   - Geometric values, ie tangent is changed.

*
* METHOD     : 
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
  int kdim2=2;
  int kdim3=3;
  double diffv[3];
  int ki;
  double dot;
  /* ______________________ */
  
  if (psource->iinter == SI_ORD)
    sh6putsing(psource, pdest);
  else
    {
       kdim2 = 2;
       s6diff(pdest->geo_track_2d_1,psource->geo_track_2d_1,kdim2,diffv);
       dot = s6scpr(pdest->geo_track_2d_1+kdim2, diffv, kdim2);
       if (dot * seq < 0)
	 {
	    /* Turn direction */
	    for (ki=0;ki<kdim2;ki++) 
	      {
		 pdest->geo_track_2d_1[kdim2+ki] = -pdest->geo_track_2d_1[kdim2+ki];
		 pdest->geo_track_2d_2[kdim2+ki] = -pdest->geo_track_2d_2[kdim2+ki];
	      }  
	    
	    for (ki=0;ki<kdim3;ki++) 
	      pdest->geo_track_3d[kdim3+ki] = -pdest->geo_track_3d[kdim3+ki];
	 }
       
    }
}
