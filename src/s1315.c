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
 * $Id: s1315.c,v 1.2 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1315

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1315(SISLSurf *ps1,double *ecentr,double aradiu,int idim,
	   double aepsco,double aepsge,double amax,
	   SISLIntcurve *pintcr,int icur,int igraph,int *jstat)
#else
void s1315(ps1,ecentr,aradiu,idim,aepsco,aepsge,amax,pintcr,icur,igraph,jstat)
     SISLSurf     *ps1;
     double          *ecentr;
     double          aradiu;
     int             idim;
     double          aepsco;
     double          aepsge;
     double          amax;
     SISLIntcurve *pintcr;
     int             icur;
     int             igraph;
     int             *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To march an intersection curve desribed by parameter pairs
*              in an intersection curve object, a B-spline surface and
*              a sphere.
*
*
* INPUT      : ps1    - Pointer to surface.
*              ecentr - Center of the sphere
*              aradiu - Radius of sphere
*              idim   - Dimension of the space in which the plane lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              amax   - Maximal allowed step length. If amax <=aepsge
*                       amax is neglected.
*              icur   - Indicator telling if a 3-D curve is to be made 
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:   
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through
*                            s6move and s6line.
*
*
*
* INPUT/OUTPUT:pintcr - The intersection curve. When comming as input
*                       only parameter values it the parameter plane
*                       exist. When comming as output the 3-D geometry
*                       and possibly the curve in the parameter plane
*                       of the surface is added.
*
* OUTPUT:      jstat  - status messages  
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.*                                         = 3      : Marching not succeded
*                         = 0      : ok
*                         < 0      : error
*
*
* METHOD     : An implicit description of the sphere is made and then
*              a routine for intersecting implicit represented geometry
*              by a B-spline surface is used.
*
* REFERENCES :
*
*-
* CALLS      : s6err, s1313, s1321
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 2. July 1988
*
*********************************************************************
*/
{            
  int kpos=0;         /* Position of error                                  */
  int kdeg=2;         /* The degree of the implicit equation of the plane   */
  int knumb=1;        /* Number of implicit representations to be made      */
  int kstat;          /* Local status variable                              */
  double simpli[16];  /* Array containing the implicit description of sphere*/
  
  
  if (idim != 3) goto err104;
  
  /* Make description of sphere */
  
  s1321(ecentr,aradiu,idim,knumb,simpli,&kstat);
  if (kstat < 0) goto error;                              
  
  /* Make intersection of implicit surface and B-spline surface */
  
  s1313(ps1,simpli,kdeg,aepsco,aepsge,amax,pintcr,icur,igraph,&kstat);
  if (kstat == -185) goto err185;
  if (kstat < 0) goto error;
  
  *jstat = kstat;
  goto out;
  
  /* Dimension not 3 */
  
 err104: 
  *jstat = -104;                
  s6err("s1315",*jstat,kpos);
  goto out;
  
  /* Couldn't march */
  
 err185:
  *jstat = -185;
  goto out;
  
  /* Error in lower level routine.  */
  
 error: 
  *jstat = kstat;     
  s6err("s1315",*jstat,kpos);
  goto out;
  
 out:
  return;
}                                               
