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
 * $Id: s6norm.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */
#define S6NORM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6norm(double e1[],int idim,double e2[],int *jstat)
#else
double s6norm(e1,idim,e2,jstat)
     double e1[];
     int    idim;
     double e2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To normalize a vector to length 1
*
* INPUT      : e1      - The vector to be normalized
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : jstat   - Status message
*                         0 - The length of the vector is zero
*                         1 - The length of the vector is one
*              s6norm  - The actual length of the vector
*              e2      - The nomalized version of the vector
*
* METHOD     : The length of the input vector is calulated, and the
*              output is assigned the values of the input vector and 
*              divided by the length.
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*                                  
*********************************************************************
*/
{
  register int ki;      /* Running variable in loop */
  register double tsum=DZERO; /* Dummy variables in summing loop */
  
  /* If the dimension is 1 the length of the vector is the same as the
   *  absolute value of the number */
  
  if (idim == 1)
    tsum = fabs(e1[0]);
  else
    {
      for (ki=0;ki<idim;ki++)
	tsum += (e1[ki]*e1[ki]);
      
      tsum = sqrt(tsum);
    }
  
  if (DNEQUAL(tsum,DZERO))
    for (ki=0;ki<idim;ki++)
      e2[ki] = e1[ki]/tsum;
  else
    {
      for (ki=0;ki<idim;ki++)
        e2[ki] = DZERO;
      goto mes00;
    }
  goto mes01;

  /* Length of vector is zero    */

 mes00: 
  *jstat = 0;
       goto out;

  /* Length of vector different from zero   */

 mes01: 
  *jstat = 1;
  goto out;
out: return(tsum);
}
