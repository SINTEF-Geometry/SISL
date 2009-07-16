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
 * $Id: s1856.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1856

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1856(SISLSurf *ps1,double epoint[],double edir[],int idim,
	   double aepsco,double aepsge,int *jpt,double **gpar,
	   int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1856(ps1,epoint,edir,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   epoint[];
     double   edir[];
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
*              and an infinite straight line.
*
*
*
* INPUT      : ps1    - Pointer to surface.
*              epoint - SISLPoint on the line.
*              edir   - Direction vector of the line.
*              idim   - Dimension of the space in which the line lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
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
* METHOD     : The line is described as the intersection between two
*              planes. The vertices of the surface are put into the equation 
*              of this planes achieving a surface in the two-dimentional 
*              space. Then the zeroes of this surface is found.
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1856, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
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

  sh1856(ps1,epoint,edir,idim,aepsco,aepsge,trackflag,&jtrack,
	 &wtrack,jpt,gpar,&pretop,jcrv,wcurve,&kstat);
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
    s6err("s1856",*jstat,kpos);
    goto out;

  out:
    return;
}                                               
