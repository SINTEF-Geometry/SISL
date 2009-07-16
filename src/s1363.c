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
 * $Id: s1363.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1363

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1363(SISLCurve *pc,double *cmin,double *cmax,int *jstat)
#else
void s1363(pc,cmin,cmax,jstat)
     SISLCurve  *pc;
     double *cmin;
     double *cmax;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline curve
*
* INPUT      : pc     - The B-spline curve.   
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              cmin   - Start of parametrization of curve
*              cmax   - End of parametrization of curve
*
* METHOD     : 
*
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6err
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kstat;          /* Local status variable                           */
  int kpos=0;         /* Position of error                               */
  
  /* Check if curve is correct */
  
  s1707(pc,&kstat);
  if (kstat<0) goto error;
  
  /* Pick parametrization */
  
  *cmin = pc->et[pc->ik - 1];
  *cmax = pc->et[pc->in];
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1363",*jstat,kpos);
  goto out;
 out:
  
  return;
}          
