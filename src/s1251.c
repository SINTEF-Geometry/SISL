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
 * $Id: s1251.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1251

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1251(SISLCurve *pcurve,double aepsco,double *clength,int *jstat)
#else
void s1251(pcurve,aepsco,clength,jstat)
     SISLCurve  *pcurve;
     double aepsco;
     double *clength;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Calculate the length of a B-spline curve. The length
*              calculated will not deviate more than aepsco from the 
*              real length of the curve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*              aepsco - Computer resolution.
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
* CALLS      : s6dist,s1251,s1710,freeCurve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int ki;                   /* Counter.                                    */
  int kn;                   /* Number of vertices of curve.                */
  int kdim;                 /* Dimension of the space in which the curve lies.*/
  double tpol = (double)0.0;/* Length of control polygon.                  */
  double tdist =(double)0.0;/* Distance between first and last vertex.     */
  double tdum;              /* Help variable.                              */
  double tmid;              /* Midpoint of parameter interval of curve.    */
  double tlength1,tlength2; /* Length of sub-curves.                       */
  double *s1;               /* Pointer used to traverse coefficient array. */
  SISLCurve *qc1 = NULL;        /* First sub-curve.                            */
  SISLCurve *qc2 = NULL;        /* Second sub-curve.                           */
  
  /* Copy properties of curve to local parameters. */
  
  kn = pcurve -> in;
  kdim = pcurve -> idim;
  
  /* Calculate length of control polygon. */
  
  for (ki=1,s1=pcurve->ecoef+kdim; ki<kn; ki++,s1+=kdim)
    tpol += s6dist(s1-kdim,s1,kdim);
  
  /* Calculate distance from first to last vertex.  */
  
  tdist = s6dist(pcurve->ecoef,pcurve->ecoef+(kn-1)*kdim,kdim);
  
  /* Test if the length of the curve can be approximated by the
     distance from the first to the last vertex.                */
  
  if (DEQUAL(tpol + tdist,(double)0.0)) tdum = (double)0.0;
  else tdum = (tpol - tdist)/(tpol + tdist);
  
  if (tdum < aepsco)
    
    /* The length of the curve is found.  */
    
    *clength = tdist;
  else
    {
      
      /* Subdivide the curve at the midpoint. */
      
      tmid = ((double)0.5)*(*(pcurve->et+pcurve->ik-1) + *(pcurve->et+kn));
      
      s1710(pcurve,tmid,&qc1,&qc2,&kstat);
      if (kstat < 0) goto error;
      
      /* Compute length of first sub-curve. */
      
      s1251(qc1,aepsco,&tlength1,&kstat);
      if (kstat < 0) goto error;
      
      /* Compute length of second sub-curve.  */
      
      s1251(qc2,aepsco,&tlength2,&kstat);
      if (kstat < 0) goto error;
      
      *clength = tlength1 + tlength2;
    }
  
  /* Length of curve found.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1251",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by sub-curves.  */
  
  if (qc1 != NULL) freeCurve(qc1);
  if (qc2 != NULL) freeCurve(qc2);
  
  return;
}

