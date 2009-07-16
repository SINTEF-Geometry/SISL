/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1986

#include "sislP.h"                                                 


#if defined(SISLNEEDPROTOTYPES)
void s1986(SISLCurve *pc, double aepsge, int *jgtpi, double **gaxis,
	   double *cang,int *jstat)
#else
void s1986(pc,aepsge,jgtpi,gaxis,cang,jstat)
     SISLCurve *pc;
     double aepsge;
     int   *jgtpi;
     double **gaxis;
     double *cang;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the direction cone of a curve.
*
*
* INPUT      : pc        - Curve to treat.
*              aepsge    - Geometry tolerance.
*
* OUTPUT     : jgtpi     - To mark if the angle of the direction cone is
*                          greater than pi.
*                           0 - The direction cone of the curve
*                               is not greater than pi. 
*                           1 - The direction cone of the curve
*                               is greater than pi.
*              gaxis     - Allocated array containing the coordinates of the
*                          center of the cone. It is only computed if
*                          *jgtpi = 0.
*              cang      - The angle from the center to the boundary of the
*                          cone. It is only computed if *jgtpi = 0.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 9403.
*
*********************************************************************
*/                                     
{
   int kstat = 0;        /* Local status variable.  */
   int kpos = 0;
   int kdim = pc->idim;
   
   /* Allocate scratch for the output array. */
   
   if ((*gaxis = newarray(kdim, DOUBLE)) == SISL_NULL) goto err101;
   
   /* Let s1991 compute the cone. */
   
   s1991(pc, aepsge, &kstat);
   if (kstat < 0) goto error;
   
   /* Copy the resulting cone to output parameters. */
   
   *jgtpi = (pc->pdir->igtpi > 0) ? 1 : 0;
   *cang = pc->pdir->aang;
   memcopy(*gaxis, pc->pdir->ecoef, kdim, DOUBLE);
   
  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
  err101: 
     *jstat = -101;
  s6err("s1986",*jstat,kpos);
  goto out;
    
  /* Error in lower level routine. */
  
  error:
     *jstat = kstat;
  s6err("s1986",*jstat,kpos);
  goto out;
     
  
  out: 
    return;
}
