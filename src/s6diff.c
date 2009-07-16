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
 * $Id: s6diff.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DIFF

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6diff(double e1[],double e2[],int idim,double e3[])
#else
void s6diff(e1,e2,idim,e3)
     double e1[];
     double e2[];
     double e3[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the difference between two vectors
*
* INPUT      : e1      - The first vector
*              e2      - The second vector
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : 
*              e3      - The difference of e1 and e2
*
* METHOD     : The difference is calculated by vetor substraction
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-june-1988
*                                  
*********************************************************************
*/
{
  int ki;
  for (ki=0;ki<idim;ki++)
    e3[ki] = e1[ki] - e2[ki];
  return;
}
