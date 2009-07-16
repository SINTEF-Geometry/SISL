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
 * $Id: s1372.c,v 1.4 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1372

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1372(SISLCurve *pc1,double epoint[],double edirec[],double aradiu,
	   int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1372(pc1,epoint,edirec,aradiu,idim,aepsco,aepsge,
           jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     double   epoint[];
     double   edirec[];
     double   aradiu;
     int      idim;
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find all intersections between a curve and a cylinder if
*              the dimension is equal to three and a curve and a circle
*              if the dimension is two.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              epoint - SISLPoint on cylinder axis
*              edirec - Direction of cylinder axis
*              aradiu - Radius of the circle or sphere
*              idim   - Dimension of the space in which the cylinder/circle
*                       lies. idim should be equal to two or three.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : *jpt   - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the curve. The points lie continuous.
*                       Intersection curves are stored in wcurve.
*              *jcrv  - Number of intersection curves.
*              wcurve  - Array containing descriptions of the intersection
*                       curves. The curves are only described by points
*                       in the parameter interval. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The vertices of the curve is put into the equation of the
*              sphere/circle achieving a curve in the one-dimentional space.
*              The zeroes of this curve is found.
*
*
* REFERENCES : Main routine written by Tor Dokken, SI, 1988.
*
* CALLS      : sh1371,sh1372, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Updated to
*              handle 2D input as specified in header.
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


  if ( idim == 2 )
    sh1371(pc1,epoint,aradiu,idim,aepsco,aepsge,
	   trackflag,&jtrack,&wtrack,jpt,gpar,&pretop,jcrv,wcurve,&kstat);
  else
    sh1372(pc1,epoint,edirec,aradiu,idim,aepsco,aepsge,
	   trackflag,&jtrack,&wtrack,jpt,gpar,&pretop,jcrv,wcurve,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  /*
   * Intersections found.
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1372",kstat,kpos);
    goto out;

  out:
    return;
}
