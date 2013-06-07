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
 * $Id: maketracks.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define MAKE_TRACKS
#include "sislP.h"

#if defined (SISLNEEDPROTOTYPES)
void make_tracks (SISLObject * po1, SISLObject * po2, int ideg,
		  double eimpli[], int icrv, SISLIntlist ** vlist,
		  int *jtrack,SISLTrack *** wcrv, double aepsge, int *jstat)
#else
void make_tracks (po1, po2, ideg,eimpli, icrv, vlist,
		  jtrack,wcrv, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     int ideg;
     double eimpli[];
     int icrv;
     SISLIntlist **vlist;
     int *jtrack;
     SISLTrack ***wcrv;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : An empty function, (to ensure similarity with other
*              versions).
*
*
* INPUT      : po1    - Pointer first object.
*              po2    - Pointer second object.
*                       internal format.
*              icrv   - Number of lists in vlist.
*              vlist  - Array representing intersection curves on the
*                       internal format.
*              ideg   - Type of track
*                           = 0, Bspline vs Bspline
*                           = 1, Bspline vs Plane
*                           = 2, Bspline vs Quadric surface
*                           = 1001 Bspline vs Torus surface
*                           = 1003 Bspline silhouette line, parallel projection
*                           = 1004 Bspline silhouette line, perspective projection
*                           = 1005 Bspline silhouette line, circular projection
*
*              eimpli[16]  Description of the implicit surface.
*              aepsge - Geometry tolerance
*
*OUTPUT:       jtrack - No of tracks made.
*              wtrack - Array containing pointers to tracks.
*              jstat - status messages
*                       >0:warning
*                       = 0:ok
*                       <0:error
*
*
*
*METHOD:   Refining datapoints and
*          Hermite interpolation using curvature radius.
*
*REFERENCES:
*
*
*-
*CALLS:     control_and_refine_ss(_si)
*          s1359 - Hermite interpolation of curve using curvature radius.
*
*WRITTEN BY:Ulf J.Krystad, SI, 30.06 .91.
*
*********************************************************************
*/
{

  *jstat  = 0;
  *jtrack = 0;

}

