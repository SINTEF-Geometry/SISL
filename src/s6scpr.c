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
 * $Id: s6scpr.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6SCPR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6scpr(double e1[],double e2[],int idim)
#else
double s6scpr(e1,e2,idim)
     double e1[];
     double e2[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make the scalar product of two vectors
*
* INPUT      : e1      - The first vector in the scalar product
*              e2      - The second vector in the scalar product
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : s6scpr  - The value of the scalar product
*
* METHOD     : The scalar product is calculated according to the definition
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*                                  
*********************************************************************
*/
{
  register int ki;
  register double tsum=DZERO; 

  for (ki=0;ki<idim;ki++)
    tsum += e1[ki]*e2[ki];

  return(tsum);
}
