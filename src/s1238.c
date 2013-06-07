/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: s1238.c,v 1.2 2001-03-19 15:58:42 afr Exp $
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
  double *spar1 = SISL_NULL;    /* Values of constant parameter curves
			      in first parameter direction.        */
  double *spar2 = SISL_NULL;    /* Values of constant parameter curves
			      in second parameter direction.       */
  double *spoint = SISL_NULL;   /* Sequence of straight lines 
			      approximating a curve.               */
  SISLCurve *qc = SISL_NULL; /* Constant parameter curve.            */
  SISLCurve *uc[20];    /* Constant parameter curves.           */
  
  for (ki=0; ki<20; ki++) uc[ki] = SISL_NULL;
  
  /* Test dimension of curve and surface.  */
  
  if (pcurve -> idim != 2) goto err108;
  if (psurf -> idim != 3) goto err104;
  
  /* Test if the limiting curve is closed.  */
  
  s1364(pcurve,aepscu,&kstat);
  if (kstat < 0) goto error;
  if (kstat != 1) goto err114;
  
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
	  
	  if (uc[kj] != SISL_NULL) freeCurve(uc[kj]);   uc[kj] = SISL_NULL;
	  if (spoint != SISL_NULL) freearray(spoint);  spoint = SISL_NULL;
	}
      kcurve = 0;
      if (qc != SISL_NULL) freeCurve(qc);   qc = SISL_NULL;
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
	  
	  if (uc[kj] != SISL_NULL) freeCurve(uc[kj]);   uc[kj] = SISL_NULL;
	  if (spoint != SISL_NULL) freearray(spoint);  spoint = SISL_NULL;
	}
      kcurve = 0;
      if (qc != SISL_NULL) freeCurve(qc);  qc = SISL_NULL;
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
  
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);
  if (spoint != SISL_NULL) freearray(spoint);
  if (qc != SISL_NULL) freeCurve(qc);
  
  return;
}
