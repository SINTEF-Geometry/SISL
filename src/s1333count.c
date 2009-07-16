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
 * $Id: s1333count.c,v 1.2 2001-03-19 15:58:45 afr Exp $
 *
 */


#define S1333_COUNT
#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1333_count(int inbcrv,SISLCurve *vpcurv[],int *jcont,int *jstat)
#else
void s1333_count(inbcrv,vpcurv,jcont,jstat)
     int    	inbcrv;
     SISLCurve  *vpcurv[];
     int        *jcont;
     int        *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To count the continuity at the start and end of the curves
*              based on the multiplicity of knots and return the continuity.
*
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcurv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              *jcount - The continuity at the start or end based on
*                        multiplicity of knots
*-
* CALLS      : 
*
* WRITTEN BY : Tor Dokken  SI  Oslo,Norway.  Feb 1992
*
*********************************************************************
*/
{
  int kmult1,kmult2,kmult;   /* Multiplicities */
  int kcont=0;               /* Continuity so far */
  int kpos=0;
  int kleft = 0;
  int kstat;
  int ki;
  SISLCurve *curve=SISL_NULL;     /* Pointer to curve being tested */

  *jcont = -1;
  

  for (ki=0 ; ki<inbcrv ; ki++)
    {
       curve = vpcurv[ki];
       kmult1 = s6knotmult(curve->et,curve->ik,curve->in,&kleft,
                           curve->et[curve->ik-1],&kstat);
       if (kstat<0)goto error;
       
       kmult2 = s6knotmult(curve->et,curve->ik,curve->in,&kleft,
                           curve->et[curve->in],&kstat);
       if (kstat<0)goto error;
       kmult = MAX(kmult1,kmult2);
       kmult = MIN(kmult,curve->ik);
       
       if (ki==0)
	 {
	   
	   kcont = curve->ik - kmult - 1;
	 }
       else
	 {   
           kcont = MIN(kcont, curve->ik - kmult - 1); 
	 }
     }
  

  /* Task done */
  
  *jcont = kcont;
  
  *jstat = 0;
  goto out; 
  
  /* Error in lower level routine.  */

  error : 
    *jstat = kstat;     
  s6err("s1333_count",*jstat,kpos);
  goto out;
 out:
   
  return;
}
