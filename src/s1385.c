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
 * $Id: s1385.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1385

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1385(double ep0[],double ept[],double ep1[],double as,int idim,
	   double aepsge,SISLCurve **rc,int *jstat)      
#else
void s1385(ep0,ept,ep1,as,idim,aepsge,rc,jstat)
     double ep0[];
     double ept[];
     double ep1[];
     double as;
     int    idim;
     double aepsge;
     SISLCurve  **rc;
     int    *jstat;      
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a conic arc by the start point, the inter-
*              section point of the tangents, the end point and a shape
*              factor. The conic arc is converted to a spline curve.
*             
*
* INPUT      : ep0    - Start point of segment
*              ept    - Intersection point of tangents
*              ep1    - End point of segment
*              as     - Shape factor
*                        as < 0.5 a ellipse
*                        as = 0.5 a parabola
*                        as > 0.5 a hyperbola
*                        as >= 1  The start and end point lies on
*                                 different branches of the parabola.
*                                 If as>=1 then as is put to 0.999999
*              idim   - The dimension of the curve to be produced
*              aepsge - Tolerance to be used for the conversion
*
*
*
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rc     - Pointer to the curve produced
*
* METHOD     : The conic is made as a rational B-spline curve according
*              to the following formula:
*
*                     ep0 (1-t)(1-t) + 2 s/(1-s) t(1-t) ept + tt ep1
*              p(t) = ----------------------------------------------
*                         (1-t)(1-t) + 2 s/(1-s) t(1-t) + tt        
*
*              The the conic arc is converted to a B-spline curve within
*              the tolerance specified.
*
* REFERENCES :
*
*-                                      
* CALLS      : 
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, May 1988
*
*********************************************************************
*/
{
  int kk=3;              /* Polynomial order                       */
  int kn=3;              /* Number of vertices                     */
  int kind=4;            /* Rational Bezier curve                  */
  int kpos=0;            /* Position of error                      */
  int kstat;             /* Local status variable                  */
  int ki;                /* Lopp variable                          */
  double tdum;           /* Homogenous coordinate                  */
  double toffset=DZERO;  /* Approximation is with zero offset      */
  double tmax;           /* Maximal step length                    */
  double st[6];          /* Knots                                  */
  double *scoef=SISL_NULL;    /* Vertices                               */
  double *sdum=SISL_NULL;     /* Dummy normal vector                    */
  SISLCurve *qc=SISL_NULL;        /* Pointer to rational description of conic */
  
  /* Allocate scrath for vertex vector */
  
  scoef = newarray((idim+1)*3,DOUBLE);
  if (scoef==SISL_NULL) goto err101;
  sdum  = new0array(idim,DOUBLE);
  if (sdum == SISL_NULL) goto err101;
  
  if (as>=(double)1.0) as = (double)0.9999999;
  memcopy(scoef,ep0,idim,DOUBLE);
  scoef[idim] = (double)1.0;
  tdum = as/((double)1.0-as);
  scoef[2*idim+1] = tdum;
  for (ki=0;ki<idim;ki++)
    {
      if (DEQUAL(tdum,DZERO))
        {
	  scoef[idim+1+ki] = ept[ki];
        }
      else
        {
	  scoef[idim+1+ki] = tdum*ept[ki];
        }
    }
  
  memcopy(scoef+2*idim+2,ep1,idim,DOUBLE);
  scoef[3*idim+2] = (double)1.0;
  
  st[0] = DZERO;
  st[1] = st[0];                                                     
  st[2] = st[0];
  st[3] = (double)1.0;
  st[4] = st[3];
  st[5] = st[4];
  
  qc = newCurve(kn,kk,st,scoef,kind,idim,1);
  if (qc == SISL_NULL) goto err101;
  
  /* Convert to spline representation */
  
  tmax = s6dist(ep0,ep1,idim);
  
  s1360(qc,toffset,aepsge,sdum,tmax,idim,rc,&kstat);
  if (kstat<0) goto error;
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("S1385",*jstat,kpos);
  goto out;
  
  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("S1385",*jstat,kpos);
  goto out;
  
  
  
  /* Free allocated arrays */
 out:
  
  if (scoef != SISL_NULL) freearray(scoef);
  if (sdum  != SISL_NULL) freearray(sdum);
  if (qc    != SISL_NULL) freeCurve(qc);
  
  return;
}
    
