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
 * $Id: s1399.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1399

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1399(SISLCurve *pc,double astart,double astop)
#else
void s1399(pc,astart,astop)
     SISLCurve  *pc;
     double astart;
     double astop;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Change the knotvector to go from astart to astop.
*
* INPUT      : pc      - The curve.
*              astart  - Parametervalue at new startpoint.
*              astop   - Parametervalue at new endpoint.
*
*-
* CALLS      :
*
* WRITTEN BY : Morten Daehlen, SI, 88-09.
*
********************************************************************/
{
  int  kk= pc->ik;             /* Order of the input curve.             */
  int  kn= pc->in;             /* Number of vertices in the input curve.*/
  double *st=SISL_NULL;             /* Pointers used in loop.                */ 
  double a,b;
  int ii, kpos=0, kstat=0;
  if (!pc) goto out;
  
  if((st = newarray(kk+kn,DOUBLE)) == SISL_NULL) goto err101;
  
  a = pc -> et[kk-1];
  b = pc -> et[kn];
  
  for (ii=0;ii<kn+kk;ii++)
    st[ii] = astart+(((pc -> et[ii])-a)/(b-a))*(astop-astart);
  for (ii=0;ii<kn+kk;ii++)
    pc -> et[ii] = st[ii];
  goto out;
  
  /* Error in scratch allocation */
  
  err101: 
    kstat = -101;
    s6err("s1399",kstat,kpos);
    goto out;
      
  out:
    if (st != SISL_NULL) freearray(st);
    return;
}
