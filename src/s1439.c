/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1439.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1439

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1439(SISLSurf *ps1,double apar,int idirec,SISLCurve **rcurve,int *jstat)
#else
void s1439(ps1,apar,idirec,rcurve,jstat)
     SISLSurf   *ps1;
     double apar;
     int idirec;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make constant parameter curve in the surface. The 
*              constant parameter value used is apar and is in the 
*              idirec parameter direction.
*
*
*
* INPUT      : ps1    - Surface.
*              apar   - Parameter value to use when picking out constant
*                       parameter curve.
*              idirec - Parameter direction in which to pick (must be 1 or 2)
*
*
* OUTPUT     : rcurve - Constant parameter curve.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES :
*
* CALLS      : s1436, s1437 - These two routines do the job, which
*                             one is called depends on what parameter
*                             direction to pick from.
*
* WRITTEN BY : Christophe Rene Birkelan, SINTEF Oslo, July 1993.
*
*********************************************************************
*/                                     
{
  int kpos = 0;      /* Position of error.                            */

  if(idirec == 1)
    {
      s1437(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else if(idirec == 2)
    {
      s1436(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else
    goto err115;

  /* Success !  Curve picked */

  goto out;


  /* Error in input parameter idirec.  */
  
  err115: 
    *jstat = -115;
    s6err("s1439",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine.  */
  
  error: 
    s6err("s1439",*jstat,kpos);
    goto out;
  
  out: 
    return;
}
