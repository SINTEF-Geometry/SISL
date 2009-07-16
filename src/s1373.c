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
 * $Id: s1373.c,v 1.3 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1373

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1373(SISLCurve *pc1,double etop[],double eaxis[],double econe[],int idim,
	   double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1373(pc1,etop,eaxis,econe,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     double   etop[];
     double   eaxis[];
     double   econe[];
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
* PURPOSE    : Find all intersections between a curve and a cone.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              etop   - Top point on cone.
*              eaxis  - Point on the cone axis (other than top point).
*              econe  - SISLPoint on cone surface.
*              idim   - Dimension of the space in which the cone
*                       lies. idim should be equal to three.
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
*              cone achieving a curve in the one-dimentional space.
*              The zeroes of this curve is found.
*
*
* REFERENCES : Main routine written by Tor Dokken, SI, 1988.
*
* CALLS      : sh1373, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Updated
*              header to reflect that only 3D input is accepted.
*
*********************************************************************
*/
{
  int kstat = 0;              /* Local status variable.                      */
  int kpos = 0;               /* Position of error.                          */
  int trackflag = 0;
  int jtrack;
  SISLTrack **wtrack=SISL_NULL;
  int *pretop=SISL_NULL;

  sh1373(pc1,etop,eaxis,econe,idim,aepsco,aepsge,
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
    s6err("s1373",*jstat,kpos);
    goto out;

  out:
    return;
}
