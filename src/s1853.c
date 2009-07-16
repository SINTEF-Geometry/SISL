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
 * $Id: s1853.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1853

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1853(SISLSurf *ps1,double epoint[],double edirec[],double aradius,
	   int idim,double aepsco,double aepsge,int *jpt,double **gpar,
	   int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1853(ps1,epoint,edirec,aradius,idim,aepsco,aepsge,
	   jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   epoint[];
     double   edirec[];
     double   aradius;
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
* PURPOSE    : Find all intersections between a tensor-product surface
*              and a cylinder.
*
*
*
* INPUT      : ps1      - Pointer to surface.
*              epoint   - SISLPoint on the axis of the cylinder.
*              edirec   - The direction vector of the axis of the cylinder.
*              aradius  - Radius of the cylinder.
*              idim     - Dimension of the space in which the cylinder lies.
*              aepsco   - Computational resolution.
*              aepsge   - Geometry resolution.
*                                                                        
*
*
* OUTPUT     : *jpt   - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       plane of the surface. The points lie continuous. 
*                       Intersection curves are stored in wcurve.
*              *jcrv  - Number of intersection curves.
*              wcurve  - Array containing descriptions of the intersection
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
* METHOD     : The vertices of the surface are put into the equation of the
*              cylinder achieving a surface in the one-dimentional space.
*              Then the zeroes of this surface is found.
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1853, s6err.
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

  sh1853(ps1,epoint,edirec,aradius,idim,aepsco,aepsge,trackflag,&jtrack,
	 &wtrack,jpt,gpar,&pretop,jcrv,wcurve,&jsurf,&wsurf,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  for(i=0; i<jsurf; i++)
    freeIntsurf(wsurf[i]);
  if(wsurf != SISL_NULL) freearray(wsurf);

  if(jsurf > 0) 
    *jstat=10;
  else 
    *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error : 
    *jstat = kstat;
    s6err("s1853",*jstat,kpos);
    goto out;

  out:
    return;
}                                               
