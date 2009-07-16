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
 * $Id: s1603.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1603

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1603(SISLSurf *psurf,double *cmin1,double *cmin2,double *cmax1,double *cmax2,int *jstat)
#else
void s1603(psurf,cmin1,cmin2,cmax1,cmax2,jstat)
     SISLSurf   *psurf;
     double *cmin1;
     double *cmin2;
     double *cmax1;
     double *cmax2;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline surface
*
* INPUT      : pc     - The B-spline surface.   
*
* OUTPUT     : cmin1  - Start parameter in first parameter directon. 
*              cmin2  - Start parameter in second parameter directon.
*              cmax1  - End   parameter in first parameter directon. 
*              cmax2    End   parameter in second parameter directon.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error  
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
* WRITTEN BY : Qyvind Hjelle SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kpos=0;              /* Position of error          */
  
  /* Check surf pointer */
  
  if (!psurf) goto err118;
  
  /* Pick parametrization */
  
  *cmin1 = psurf->et1[psurf->ik1-1];
  *cmax1 = psurf->et1[psurf->in1];
  *cmin2 = psurf->et2[psurf->ik2-1];
  *cmax2 = psurf->et2[psurf->in2];
  
  *jstat = 0;
  goto out;
  
  /* Error in input, no B-spline surface given */
  
 err118: 
  *jstat = -118;
  s6err("s1603",*jstat,kpos);
  goto out;
  
 out:
  
  return;
}          
