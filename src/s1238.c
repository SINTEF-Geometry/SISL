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
 * $Id: s1238.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1238

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1238(SISLSurf *psurf,SISLCurve *pcurve,int inmb1,int inmb2,
	   double aepsco,double aepscu,int *jstat)
#else
void s1238(psurf,pcurve,inmb1,inmb2,aepsco,aepscu,jstat)
     SISLSurf   *psurf;
     SISLCurve  *pcurve;
     int    inmb1;
     int    inmb2;
     double aepsco;
     double aepscu;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Draw constant parameter lines in a B-spline surface
*              that is limited by a closed B-spline curve lying in 
*              the parameter plane.
*
*
*
* INPUT      : psurf  - Pointer to the surface.
*              pcurve - SISLCurve limiting the part of the surface that
*                       is to be drawn.
*              inmb1  - Number of constant parameter lines to be drawn
*                       in first parameter direction.
*              inmb2  - Number of constant parameter lines to be drawn
*                       in second parameter direction.
*              aepsco - Computer resolution.
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
*              s1239  - Pick segments of curve that lies within the
*                       curve pcurve.
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
  int ki,kj;               /* Counters.                            */
  int knbpnt;              /* Number of points in line sequence.   */
  int kmax = 20;           /* Maximum number of curves to draw 
			      along one parameter line.            */
  int kcurve = 0;          /* Number of curves along one parameter
			      line.                                */
  double *spar1 = NULL;    /* Values of constant parameter curves
			      in first parameter direction.        */
  double *spar2 = NULL;    /* Values of constant parameter curves
			      in second parameter direction.       */
  double *spoint = NULL;   /* Sequence of straight lines 
			      approximating a curve.               */
  SISLCurve *qc = NULL; /* Constant parameter curve.            */
  SISLCurve *uc[20];    /* Constant parameter curves.           */
  
  for (ki=0; ki<20; ki++) uc[ki] = NULL;
  
  /* Test dimension of curve and surface.  */
  
  if (pcurve -> idim != 2) goto err108;
  if (psurf -> idim != 3) goto err104;
  
  /* Test if the limiting curve is closed.  */
  
  s1364(pcurve,aepscu,&kstat);
  if (kstat < 0) goto error;
  if (kstat != 1) goto err114;
  
  /* Allocate space for arrays containing constant parameter values. */
  
  if ((spar1 = newarray(inmb1,double)) == NULL) goto err101;
  if ((spar2 = newarray(inmb2,double)) == NULL) goto err101;
  
  /* Find parameter values to be used to make curves with constant
     parameter values in second direction.                         */
  
  s1236(psurf->et2,psurf->in2,psurf->ik2,inmb2,spar2,&kstat);
  if (kstat < 0) goto error;
  
  for (ki=0; ki<inmb2; ki++)
    {
      
      /* Pick curve with constant second parameter direction. */
      
      s1436(psurf,spar2[ki],&qc,&kstat);
      if (kstat < 0) goto error;
      
      /* Pick the part of the curves that lie within the curve pcurve. */
      
      s1239(qc,1,spar2[ki],pcurve,aepsco,aepscu,uc,kmax,&kcurve,&kstat);
      if (kstat < 0) goto error;
      
      for (kj=0; kj<kcurve; kj++)
	{
	  
	  /* Approximate the curve by a sequence of straight lines. */
	  
	  s1605(uc[kj],aepscu,&spoint,&knbpnt,&kstat);
	  if (kstat < 0) goto error;
	  
	  /* Draw the curve as a sequence of straight lines.  */
	  
	  s6drawseq(spoint,knbpnt);
	  
	  /* Prepare for next curve to draw.  */
	  
	  if (uc[kj] != NULL) freeCurve(uc[kj]);   uc[kj] = NULL;
	  if (spoint != NULL) freearray(spoint);  spoint = NULL;
	}
      kcurve = 0;
      if (qc != NULL) freeCurve(qc);   qc = NULL;
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
      
      /* Pick the part of the curves that lie within the curve pcurve. */
      
      s1239(qc,0,spar1[ki],pcurve,aepsco,aepscu,uc,kmax,&kcurve,&kstat);
      if (kstat < 0) goto error;
      
      for (kj=0; kj<kcurve; kj++)
	{
	  
	  /* Approximate the curve by a sequence of straight lines. */
	  
	  s1605(uc[kj],aepscu,&spoint,&knbpnt,&kstat);
	  if (kstat < 0) goto error;
	  
	  /* Draw the curve as a sequence of straight lines.  */
	  
	  s6drawseq(spoint,knbpnt);
	  
	  /* Prepare for next curve to draw.  */
	  
	  if (uc[kj] != NULL) freeCurve(uc[kj]);   uc[kj] = NULL;
	  if (spoint != NULL) freearray(spoint);  spoint = NULL;
	}
      kcurve = 0;
      if (qc != NULL) freeCurve(qc);  qc = NULL;
    }
  
  /* The surface is drawn.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1238",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */
  
 err104: *jstat = -104;
  s6err("s1238",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2.  */
  
 err108: *jstat = -108;
  s6err("s1238",*jstat,kpos);
  goto out;
  
  /* Error in curve-description. Open curve when close expected.  */
  
 err114: *jstat = -114;
  s6err("s1238",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
  error : *jstat = kstat;
  s6err("s1238",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local arrays etc.  */
  
  if (spar1 != NULL) freearray(spar1);
  if (spar2 != NULL) freearray(spar2);
  if (spoint != NULL) freearray(spoint);
  if (qc != NULL) freeCurve(qc);
  
  return;
}
