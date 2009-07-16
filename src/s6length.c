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
 * $Id: s6length.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */
#define S6LENGTH

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6length(double e1[],int idim,int *jstat)
#else
double s6length(e1,idim,jstat)
     double e1[]; 
     int    idim;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the length of a vector
*
* INPUT      : e1      - The vector 
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : jstat   - Status message
*                         0 - The length of the vector is zero
*                         1 - The length of the vector is not zero
*              s6norm  - The actual length of the vector
*
* METHOD     : The length of the input vector is calulated. 
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-june 1988
*                                  
*********************************************************************
*/
{
  register int ki;            /* Running variable in loop */
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
    goto mes01;

  /* Length of vector is zero    */

  *jstat = 0;
  goto out;

  /* Length of vector different from zero   */

 mes01: *jstat = 1;
        goto out;

 out: return(tsum);
}
