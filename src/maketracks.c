/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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

