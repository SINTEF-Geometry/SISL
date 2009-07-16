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
 * $Id: s1232.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1232

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1232(double et1[],int in,int ik,
	   double afak1,double afak2,double et2[],int *jstat)
#else
void s1232(et1,in,ik,afak1,afak2,et2,jstat)
     double et1[];
     int    in;
     int    ik;
     double afak1;
     double afak2;
     double et2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Produce a knot vector that is stretched in the start
*              and/or end compared to an input knot vector.
*
*
*
* INPUT      : et1    - Original knot vector.
*              in     - Number of vertices of the curve which the
*                       knot vector belongs to.
*              ik     - Order of the curve which the knot vector
*                       belongs to.
*              afak1  - How much the new knot vector is to be stretched
*                       at the start compared with the original one. 
*                       afak1 >= 0. The parameter interval is extended 
*                       afak1 times the length of the parameter interval
*                       in the start of the interval.
*              afak2  - How much the new knot vector is to be stretched
*                       at the end compared with the original one. 
*                       afak2 >= 0.
*
*
*
* OUTPUT     : et2    - The stretched knot vector.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Replace the ik first and last knots with knots at the
*              new endpoints of the knot vector. The interior knots are
*              left unchanged.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/                   
{
  int kpos = 0;   /* Position of error.                */
  int ki;         /* Counter.                          */
  double tleng;   /* Length of parameter interval.     */
  double tstart;  /* New start value of knot vector.   */
  double tend;    /* New end value of knot vector.     */
  
  /* Test input.  */
  
  if (ik < 1) goto err110;
  if (in < ik) goto err111;
  
  /* Test if knot vector degenerated.  */
  
  tleng = et1[in] - et1[ik-1];
  if (tleng <= (double)0.0) goto err112;  
  
  /* Copy input knot vector to output knot vector.   */
  
  memcopy(et2,et1,in+ik,double);
  
  if (afak1 > (double)0.0)
    {
      
      /* Extend basis at start of curve.  */
      
      tstart = et1[ik-1] - tleng*afak1;
      for (ki=0; ki<ik; ki++) et2[ki] = tstart;
    }
  
  if (afak2 > (double)0.0)
    {
      
      /* Extend basis at end of curve.  */
      
      tend = et1[in] + tleng*afak2;
      for (ki=in; ki<in+ik; ki++) et2[ki] = tend;
    }
  
  /* New knot vector produced.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Order less than 1.  */
  
 err110: *jstat = -110;
  s6err("s1232",*jstat,kpos);
  goto out;
  
  /* Error in input. Number of vertices less than order.  */
  
 err111: *jstat = -111;
  s6err("s1232",*jstat,kpos);
  goto out;
  
  /* Error in input. Knot vector degenerated.  */
  
 err112: *jstat = -112;
  s6err("s1232",*jstat,kpos);
  goto out;
  
 out: return;
}
