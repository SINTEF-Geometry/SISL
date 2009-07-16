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
 * $Id: s1240.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1240

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1240(SISLCurve *pcurve,double aepsge,double *clength,int *jstat)
#else
void s1240(pcurve,aepsge,clength,jstat)
     SISLCurve  *pcurve;
     double aepsge;
     double *clength;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Calculate the length of a B-spline curve. The length
*              calculated will not deviate more than (aepsge/the
*              length calculated) from the real length of the curve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : clength - The length of the curve.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist,s1251,make_cv_kreg.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;  /* Local status variable.                          */
  int kpos = 0;   /* Position of error.                              */
  int ki;         /* Counter.                                        */
  int kdim;       /* Dimension of the space in which the curve lies. */
  int kn;         /* Number of vertices of curve.                    */
  int kcalc;      /* Indicates if correct length of curve is found.  */
  double tlength; /* Length of curve.                                */
  double tprev;   /* Previous length of curve calculated.            */
  double teps;    /* Local tolerance.                                */
  double *s1;     /* Pointer used to traverse real array.            */
  SISLCurve *qc=SISL_NULL;  /* k-regular local curve.                     */
  
  if (pcurve->cuopen == SISL_CRV_PERIODIC)
    {
       /* Make curve k-regular. */
       
       make_cv_kreg(pcurve,&qc,&kstat);
       if (kstat < 0) goto error;
    }
  else qc = pcurve;
       
  /* Copy curve information to local parameters. */
  
  kdim = qc -> idim;
  kn   = qc -> in;
  
  /* Calculate length of control polygon.  */
  
  tlength = 0;
  for (ki=1,s1=qc->ecoef+kdim; ki<kn; ki++,s1+=kdim)
    tlength += s6dist(s1-kdim,s1,kdim);
  
  /* Set up local tolerance.  */
  
  teps = aepsge*100;
  
  kcalc = 0;
  while (kcalc == 0)
    {
      teps = teps/2.0;
      tprev = tlength;
      
      /* Compute length of curve.  */
      
      s1251(qc,teps,&tlength,&kstat);
      if (kstat < 0) goto error;
      
      /* Test if the error is within the tolerance. */
      
      if (fabs(tprev-tlength)/MAX(tprev,tlength) < aepsge) kcalc = 1;
      
    }
  
  /* Length of curve calculated. */
  
  *clength = tlength;
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine. */
  
  error : *jstat = kstat;
  s6err("s1240",*jstat,kpos);
  goto out;
  
 out: 
    
    /* Free local curve.  */
    
    if (qc != SISL_NULL && qc != pcurve) freeCurve(qc);
    
    return;
}

