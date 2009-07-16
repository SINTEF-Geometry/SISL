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
 * $Id: s1707.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1707

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1707(SISLCurve *pc,int *jstat)
#else
void s1707(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if a B-spline curve is correct.
*
* INPUT      : pc     - SISLCurve to treat.
*
* OUTPUT     : jstat     - status messages
*                        > 0      : warning (1&2 only used when 
*                                            cuopen=SISL_CRV_PERIODIC)
*                                      = 1: Cyclic but not full freedom.
*                                      = 2: Not cyclic.
*                                      = 8: Non-positive rational weights.
*                       = 0      : ok
*                       < 0      : error
*
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Christophe Rene Birkeland, SI-SINTEF, May 1993.
*
**********************************************************************/
{

  int kpos=0;              /* Position of error. */
  int kstat=0;
  int step = 0;
  register double *s1,*s2; /* Pointers used in loop. */
  
  if (!pc) goto err150;

  if (pc->ik > pc->in) goto err111;
  
  if (pc->ik <= 0) goto err110;
  
  if (pc->in <= 0) goto err159;
  
  if (pc->idim <= 0) goto err102;
  
  if (pc->et[pc->in+pc->ik-1] <= *pc->et) goto err112;
  
  for (s1=pc->et,s2=pc->et+pc->in+pc->ik-1; s1<s2; s1++)
    if (s1[1] < *s1) goto err112;

  /* Check rational coefficients */
  if(pc->ikind == 2 || pc->ikind == 4)
    {
      step = pc->idim + 1;
      for (s1 = pc->rcoef + pc->idim, s2 = pc->rcoef + pc->in*step; 
	   s1 < s2; 
	   s1+= step)
	if (*s1 <= 0) goto war08;
    }

  /* Check if curve really is cyclic */
  if(pc->cuopen == SISL_CRV_PERIODIC)
    {
      test_cyclic_knots(pc->et,pc->in,pc->ik,&kstat);
      if (kstat < 0) goto error;
      if (kstat == 0) goto war02;
      if (kstat == 1) goto war01;
    }
      

  /* Updating output. No errors ! */
  
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
  
  /* Dimension less than 1. */
  
  err102:
    *jstat = -102;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order less than 1. */
  
  err110:
    *jstat = -110;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order greater than number of vertices. */
   
  err111:
    *jstat = -111;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Error in knotvector. */
  
  err112:
    *jstat = -112;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Null pointer. */
  
  err150:
    *jstat = -150;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Number of vertices less than 1. */
  
  err159:
    *jstat = -159;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine */
      
  error:
    *jstat = kstat;
    s6err("s1707",*jstat,kpos);
    goto out;

  out: 
    return;
}
