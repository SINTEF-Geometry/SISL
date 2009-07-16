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
 * $Id: s1393.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1393

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1393(int n1,SISLCurve *pc1[],SISLCurve *sc1[],SISLCurve *ec1[],int *jstat)
#else
void s1393(n1,pc1,sc1,ec1,jstat)
     int   n1;
     SISLCurve *pc1[];
     SISLCurve *sc1[];
     SISLCurve *ec1[];
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* Purpose :  Split curve at midpoint, turn second curve, 
*            and normalize parameterinterval.
*
* Input     : pc1       - Pointers to first boundary curves
*             n1        - Number of curves
*
* Output    : sc1       - Pointers to first part of curve.
*             ec1       - Pointers to second part of curve.
*
*             jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                      
*-
* Calls      : s6err - error messages.
*              s1710 - split curve at given parameter value.
*              s1706 - Turn parametrization of curve.
*              s1399 - Normalize parameterinterval.
*
* Written by : Mortend Daehlen, SI, Aug. 88.
*
*********************************************************************
*/                                     
{
  int kpos = 0;
  int ki;
  int kstat = 0;
  double ax,astart,astop;
  SISLCurve *h1,*h2;
  
  astart = DZERO;
  astop  = (double)1.0;
  
  /* For each curve in pc1 split/turn and normalize. */
  for (ki=0;ki<n1;ki++)
    {
      
      /* Split */
      
      ax=(pc1[ki]->et[pc1[ki]->in]-(pc1[ki]->et[(pc1[ki]->ik)-1]))/(double)2.0;
      s1710(pc1[ki],ax,&h1,&h2,&kstat); 
      if (kstat < 0) goto error;
      
      /* Turn */
      
      s1706(h2);
      if (kstat < 0) goto error;
      
      /* Normalize */
      
      s1399(h1,astart,astop);
      if (kstat < 0) goto error;
      s1399(h2,astart,astop);
      if (kstat < 0) goto error;
      sc1[ki] = h1;
      ec1[ki] = h2;
    }
  
  *jstat=0;
  goto out;
  
  /* Error in lower level routine.   */
  
 error: *jstat = kstat;
  s6err("s1393",*jstat,kpos);
  goto out;
  
 out: return;
}
