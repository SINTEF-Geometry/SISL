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
 * $Id: s1503.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1503

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1503(SISLSurf *ps1,double base[],double norm[],double axisA[],
	   double alpha,double ratio,int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1503(ps1,base,norm,axisA,alpha,ratio,idim,aepsco,aepsge,
	   jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   base[];
     double   norm[];
     double   axisA[];
     double   alpha;
     double   ratio;
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
*              and a cone.
*
*
*
* INPUT      : ps1    - Pointer to the surface.
*              base   - Base point of cone
*              norm   - Direction of cone axis
*              axisA  - One of the two ellipse axis vectors
*              alpha  - The opening angle of the cone at axisA
*              ratio  - The ratio of axisA to axisB
*              idim     - Dimension of the space in which the cone lies.
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
* METHOD     : The vertices of the surface are put into the equation of 
*              the cone resulting in a surface in one-dimensional space.
*              Then the zeroes of this surface are found.
*
*
* REFERENCES : Main routine written by Mike Floater, SI, 1990.
*
* CALLS      : sh1503, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
*
*
*********************************************************************
*/
{            
  int kstat = 0;           /* Local status variable.                       */
  int kpos = 0;            /* Position of error.                           */
  int i;
  int trackflag = 0;
  int jtrack;
  int *pretop=SISL_NULL;
  int jsurf;
  SISLIntsurf **wsurf=SISL_NULL;
  SISLTrack **wtrack=SISL_NULL;

  sh1503(ps1,base,norm,axisA,alpha,ratio,idim,aepsco,aepsge,
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
    s6err("s1503",*jstat,kpos);
    goto out;

  out:
    return;
}                                               
