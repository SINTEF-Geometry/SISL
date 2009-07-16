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
 * $Id: s1237.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1237

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1237(SISLSurf *psurf,int inmb1,int inmb2,double aepscu,int *jstat)
#else
void s1237(psurf,inmb1,inmb2,aepscu,jstat)
     SISLSurf *psurf;
     int         inmb1;
     int         inmb2;
     double      aepscu;
     int         *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Draw constant parameter lines in a B-spline surface.
*
*
*
* INPUT      : psurf  - Pointer to the surface.
*              inmb1  - Number of constant parameter lines to be drawn
*                       in first parameter direction.
*              inmb2  - Number of constant parameter lines to be drawn
*                       in second parameter direction.
*              aepscu - The maximal distance allowed between the curves
*                       drawn and the surface.
*
*
*
* OUTPUT     : 
*              jstat  - status messages  
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
* CALLS      : s1236  - Find constant parameter values where a curve
*                       will be drawn.
*              s1436  - Pick curve with constant second parameter from
*                       surface.
*              s1437  - Pick curve with constant first parameter from
*                       surface.
*              s1605  - Approximate curve with a sequence of straight lines.
*              s6drawseq - Draw a sequence of straight lines.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;           /* Local status variable.               */
  int kpos = 0;            /* Position of error.                   */
  int ki;                  /* Counter.                             */
  int knbpnt;              /* Number of points in line sequence.   */
  double *spar1 = SISL_NULL;    /* Values of constant parameter curves
			      in first parameter direction.        */
  double *spar2 = SISL_NULL;    /* Values of constant parameter curves
			      in second parameter direction.       */
  double *spoint = SISL_NULL;   /* Sequence of straight lines 
			      approximating a curve.               */
  SISLCurve *qc = SISL_NULL;        /* Constant parameter curve.            */
  
  /* Test dimension of surface.  */
  
  if (psurf -> idim != 3) goto err104;
  
  /* Allocate space for arrays containing constant parameter values. */
  
  if ((spar1 = newarray(inmb1,double)) == SISL_NULL) goto err101;
  if ((spar2 = newarray(inmb2,double)) == SISL_NULL) goto err101;
  
  /* Find parameter values to be used to make curves with constant
     parameter values in second direction.                         */
  
  s1236(psurf->et2,psurf->in2,psurf->ik2,inmb2,spar2,&kstat);
  if (kstat < 0) goto error;
  
  for (ki=0; ki<inmb2; ki++)
    {
      
      /* Pick curve with constant second parameter direction. */
      
      s1436(psurf,spar2[ki],&qc,&kstat);
      if (kstat < 0) goto error;
      
      /* Approximate the curve by a sequence of straight lines. */
      
      s1605(qc,aepscu,&spoint,&knbpnt,&kstat);
      if (kstat < 0) goto error;
      
      /* Draw the curve as a sequence of straight lines.  */
      
      s6drawseq(spoint,knbpnt);
      
      /* Prepare for next curve to draw.  */
      
      if (qc != SISL_NULL) freeCurve(qc);   qc = SISL_NULL;
      if (spoint != SISL_NULL) freearray(spoint);  spoint = SISL_NULL;
    }
  
  /* Find parameter values to be used to make curves with constant 
     parameter values in first direction.                          */
  
  s1236(psurf->et1,psurf->in1,psurf->ik1,inmb1,spar1,&kstat);
  if (kstat < 0) goto error;
  
  for (ki=0; ki<inmb1; ki++)
    {
      
      /* Pick curve with constant first parameter direction.  */
      
      s1437(psurf,spar1[ki],&qc,&kstat);
      if (kstat < 0) goto error;
      
      /* Approximate the curve by a sequence of straight lines. */
      
      s1605(qc,aepscu,&spoint,&knbpnt,&kstat);
      if (kstat < 0) goto error;
      
      /* Draw the curve as a sequence of straight lines.  */
      
      s6drawseq(spoint,knbpnt);
      
      /* Prepare for next curve to draw.  */
      
      if (qc != SISL_NULL) freeCurve(qc);   qc = SISL_NULL;
      if (spoint != SISL_NULL) freearray(spoint);  spoint = SISL_NULL;
    }
  
  /* The surface is drawn.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1237",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */
  
 err104: *jstat = -104;
  s6err("s1237",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
  error : *jstat = kstat;
  s6err("s1237",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local arrays etc.  */
  
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);
  if (spoint != SISL_NULL) freearray(spoint);
  if (qc != SISL_NULL) freeCurve(qc);
  
  return;
}
