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
 * $Id: s1860.c,v 1.4 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1860

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1860(SISLSurf *ps,double eview[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1860(ps,eview,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps;
     double   eview[];
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
* PURPOSE    : Find the silhouette curves and points of a surface when
*              the surface is viewed from a specific direction. In
*              addition to the points and curves found by this routine,
*              break curves and edge-curves might be silhouette curves.
*
*
*
* INPUT      : ps  - Pointer to the surface.
*              eview  - The direction vector of the viewing.
*              idim   - Dimension of the space in which eview lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single silhouette points.
*              gpar   - Array containing the parameter values of the
*                       single silhouette points in the parameter
*                       plane of the surface. The points lie continuous.
*                       Silhouette curves are stored in wcurve.
*              jcrv   - Number of silhouette curves.
*              wcurve - Array containing descriptions of the silhouette
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
* METHOD     :
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1860, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Jan. 1995.
*              Made deallocation of 'wsurf' safer.
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

  sh1860(ps,eview,idim,aepsco,aepsge,trackflag,&jtrack,&wtrack,
	 jpt,gpar,&pretop,jcrv,wcurve,&jsurf,&wsurf,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  if ( wsurf != SISL_NULL )
  {
    for( i=0;  i < jsurf;  i++ )
      if ( wsurf[i] != SISL_NULL )  freeIntsurf(wsurf[i]);
    freearray(wsurf);
  }

  if(jsurf > 0)
    *jstat=10;
  else
    *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1860",*jstat,kpos);
    goto out;

  out:
    return;
}
