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
 * $Id: s1369.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1369

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1369(SISLSurf *ps,double ecentr[],double enorm[],double abigr,
	   double asmalr,int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1369(ps,ecentr,enorm,abigr,asmalr,idim,aepsco,aepsge,
           jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps;
     double   ecentr[];
     double   enorm[];
     double   abigr;
     double   asmalr;
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
* PURPOSE    : Find all intersections between a surface and a torus. 
*
* INPUT      : ps     - Pointer to the surface.
*              ecentr - The center of the torus (lying in the symmetry plane)
*              enorm  - Normal of symmetry plane
*              abigr  - Distance from ecentr to center circle of torus
*              asmalr - The radius of the torus surface
*              idim   - Dimension of the space in which the torus
*                       lies. idim should be equal to two or three.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the curve. The points lie continuous. 
*                       Intersection curves are stored in wcurve.
*              jcrv   - Number of intersection curves.
*              wcurve - Array containing descriptions of the intersection
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
* METHOD     : The vertices of the surface is put into the equation of the
*              torus achieving a surface in the one-dimensional space.
*              The zeroes of this surface are found.
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1369, s6err.
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

  sh1369(ps,ecentr,enorm,abigr,asmalr,idim,aepsco,aepsge,
	 trackflag,&jtrack,&wtrack,jpt,gpar,&pretop,jcrv,wcurve,
	 &jsurf,&wsurf,&kstat);
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
    s6err("s1369",kstat,kpos);
    goto out;

  out:
    return;
}                                               
                                           
                       
       
