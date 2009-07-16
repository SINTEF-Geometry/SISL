/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1708

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1708(SISLSurf *ps,int *jstat)
#else
void s1708(ps,jstat)
     SISLSurf *ps;
     int   *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if a B-spline surface is correct.
*
* INPUT      : ps      - SISLSurf to treat.
*
* OUTPUT     : jstat   - status messages
*                        > 0      : warning (1&2 only used when 
*                                            cuopen==SISL_CRV_PERIODIC in at
*                                            least one direction)
*                                      = 1: Cyclic but not full freedom.
*                                      = 2: Not cyclic.
*                                      = 8: Non-positive rational weights.
*                       = 0      : ok
*                       < 0      : error
*
* CALLS      :
*
* WRITTEN BY : Christophe Rene Birkeland, SI-SINTEF, May 1993.
*
**********************************************************************/
{
  int type = 1;
  int stat = 0;              /* Status for lower level routines   */
  int kstat1= 0;
  int kstat2= 0;
  int kpos = 2;              /* Indicator of parameter direction
			      * for error message.                */
  int step = 0;
  register double *s1,*s2;   /* Pointers used in loop.            */
  SISLCurve *curve=SISL_NULL;     /* Local curve                       */

  if(ps->ikind == 3) type = 3;


  /*
   *  Check second parameter direction 
   *  
   */

  curve = newCurve(ps->in2, ps->ik2, ps->et2, ps->ecoef,
		   type, ps->idim * ps->in1, 0);
  if(curve == SISL_NULL) goto err101;

  s1707( curve, &stat);
  if(stat<0) goto error;

  /* Free curve element used in 2. par. direction */

  freeCurve(curve);
  curve = SISL_NULL;
  

  /*
   *  Check first parameter direction
   *
   */

  kpos = 1;

  /*
   * NOT necessary to transpose coefficient array. 
   * Coefficients are not checked
   * in routine s1707 for curve-kind 1 || 3
   *
   */

  curve = newCurve(ps->in1, ps->ik1, ps->et1, ps->ecoef,
		   type, ps->idim * ps->in2, 0);
  if (curve == SISL_NULL) goto err101;

  s1707( curve, &stat);
  if(stat<0) goto error;

  kpos = 0;

  /* Check weights in case of rational surface */
  if(ps->ikind == 2 || ps->ikind == 4)
    {
      step = ps->idim + 1;
      for (s1 = ps->rcoef + ps->idim, s2 = ps->rcoef + ps->in1*ps->in2*step; 
           s1 < s2;
           s1+= step)
        if (*s1 <= 0) goto war08;
    }

  /* Check if surface really is cyclic */
  if(ps->cuopen_1 == SISL_CRV_PERIODIC || ps->cuopen_2 == SISL_CRV_PERIODIC)
    {
      kpos = 1;
      test_cyclic_knots(ps->et1,ps->in1,ps->ik1,&kstat1);
      if (kstat1 < 0) goto error;
      kpos = 2;
      test_cyclic_knots(ps->et2,ps->in2,ps->ik2,&kstat2);
      if (kstat2 < 0) goto error;
      if ((kstat1 == 0 && ps->cuopen_1 == SISL_CRV_PERIODIC)
	  || (kstat2 == 0 && ps->cuopen_2 == SISL_CRV_PERIODIC))
	goto war02;
      if ((kstat1 == 1 && ps->cuopen_1 == SISL_CRV_PERIODIC)
	  || (kstat2 == 1 && ps->cuopen_2 == SISL_CRV_PERIODIC))
	goto war01;
    }
      
  /* SUCCESS. */

  *jstat = 0;
  goto out;

  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector does not give
   * full freedom. */
  
  war01:
    *jstat = 1;
    goto out;
  
  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector not cyclic. */
  
  war02:
    *jstat = 2;
    goto out;
  
  /* Warning: Non-positive rational coefficients. */
  
  war08:
    *jstat = 8;
    goto out;
  
  /* Error in allocations */

  err101:
    *jstat = -101;
    s6err("s1708", *jstat, kpos);
    goto out;

  /* Lower level error */

  error:
    *jstat = stat;
    s6err("s1708",*jstat,kpos);
    goto out;

 out: 
  if (curve != SISL_NULL) freeCurve( curve );
  return;
}
