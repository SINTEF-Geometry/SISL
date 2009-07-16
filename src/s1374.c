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
 * $Id: s1374.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1374

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1374(SISLCurve *pc1,double earray[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1374(pc1,earray,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     double   earray[];
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
* PURPOSE    : Find all intersections between a curve and a quadric
*              curve if idim=2 and a quadric surface if idim=3.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              earray - Matrix of dimension (idim+1)x(idim+1) describing
*                       the conic surface such for idim=2:
*                                         T
*                        (x,y,1) A (x,y,1)  = 0 is the implicit equation
*                                               for the curve
*
*                       For idim=3:
*                        (x,y,z,1) A (x,y,z,1) = 0 is the implicit equation
*                                                  for the surface
*
*              idim   - Dimension of the space in which the plane/line
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
*              cone achieving a curve in the one-dimentional space.
*              The zeroes of this curve is found.
*
*
* REFERENCES : Main routine written by Tor Dokken, SI, 1988.
*
* CALLS      : sh1374, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
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

  sh1374(pc1,earray,idim,aepsco,aepsge,
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

  error: 
    *jstat = kstat;
    s6err("s1374",*jstat,kpos);
    goto out;

  out:
    return;
}                                               
                                           
                       
