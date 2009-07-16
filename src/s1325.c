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
 * $Id: s1325.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1325

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s1325(double aradiu,double angle)
#else
double s1325(aradiu,angle)
     double aradiu;
     double angle;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create the tangent length for interpolating a
*              circular arc with an almost equi-oscillating Hermit qubic
*
* INPUT      : aradiu  - The radius of the circular arc
*              angle   - The opening angle of the circular arc
*
* OUTPUT     : s1325   - The proposed tangent length
*
* METHOD     : A second degree equation giving the tanget length is
*              solved
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 30. June 1988
*                                  
*********************************************************************
*/
{
  double tcos,tsin;          /* Dummy variables                     */
  double ta,tb,tc,tl;        /* Dummy variables                     */
  double tconst = (double)1.85530139760811990992528773586425;
                             /* Constant used in the calculation    */
  
  
  
  tcos = cos(angle);
  tsin = sin(angle);
  
  /*  Calculate length of tangents
   *   tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 */
  
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)
           * tcos - (double)0.4*tconst - (double)1.0;
  tl     = aradiu*(-tb+sqrt(tb*tb-4*ta*tc))/((double)2.0*ta);
  
  return(tl);
}
