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
 * $Id: s6crss.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CRSS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6crss(double e1[],double e2[],double e3[])
#else
void s6crss(e1,e2,e3)
     double e1[];
     double e2[];
     double e3[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make the cross product of two 3-D vectors
*
* INPUT      : e1      - First 3-D vector
*              e2      - Second 3-D vector
*
* OUTPUT     : e3      - The vector containing the cross product e1xe2
*
*
* METHOD     : The cross product is calculated by using its definition.
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*
*********************************************************************
*/
{
  e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
  e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
  e3[2] = e1[0]*e2[1] - e1[1]*e2[0];
}
