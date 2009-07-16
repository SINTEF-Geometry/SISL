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
 * $Id: evalcrvarc.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define EVAL_CRV_ARC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   eval_crv_arc(SISLCurve *pc1,double center[],double radius,int ider,
	   double epar[],int *ilfs,double eder[],int *jstat)
#else
void eval_crv_arc(pc1,center,radius,ider,epar,ilfs,eder,jstat)
     SISLCurve *pc1;
     double    center[];
     double    radius;
     int       ider;
     double    epar[];
     int       *ilfs;
     double    eder[];
     int       *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the surface expressed by ((pc2x - pc1x,pc2y -
*              pc1y)*(-(dpc1/ds)y,(dpc1/ds)x), (pc2x - pc1x,pc2y - pc1y)
*              *(-(dpc2/dt)y,(dpc2/dt)x)), at the parameter values given
*              in epar, where pc2 is a circle. Compute ider derivatives.
*
*
*
* INPUT      : pc1    - Pointer to the first curve.
*              center - Center of the circle.
*              radius - Radius of the circle.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.    
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*              epar   - Parameter-value at which to calculate. Dimension
*                       of epar is 2.
*
*                
*
* INPUT/OUTPUT : ilfs  - Pointer to the interval in the knotvector
*                        in first parameter direction where epar[0] 
*                        is found. The relation
*                          
*                          et1[ilfs] <= epar[0] < et1[ilfs+1]
* 
*                        where et1 is the knotvektor should hold.
*                        ilfs is set equal to zero at the first call
*                        to the routine. 
*             
*
*
*
* OUTPUT     : eder   - Array where the derivative of the curve in
*                       apar is placed. The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative, etc. Dimension of eder is 
*                       idim*(1+2+...+(ider+1)).
*              jstat  - status messages  
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SI, March 1992. 
*
*********************************************************************
*/                                     
{
  int kstat=0;         /* Local status variable.                         */
  int kpos=0;          /* The position of error.                         */
  int kder = ider + 1; /* Number of necessary curve derivatives.         */ 
  double crv1[8];      /* The derivatives of the first curve.            */
  double crv2[8];      /* The derivatives of the second curve.           */
  double diffvec[2];   /* Difference vector between the two curves.      */
  double ang_cos;      /* Cosine of the angle in the circle.             */
  double ang_sin;      /* Sine of the angle in the circle.               */
  
  /* Check the input. */
  
  if (pc1->idim != 2) goto err102;
  if (ider > 2) goto err103;
  
  /* Evaluate the curve. */
  
  s1221(pc1, kder, epar[0], ilfs, crv1, &kstat);
  if (kstat < 0) goto error;
  
  /* Evaluate the circle. */

  ang_cos = cos(epar[1]);
  ang_sin = sin(epar[1]);
  
  crv2[0] = center[0] + radius*ang_cos;
  crv2[1] = center[1] + radius*ang_sin;
  crv2[2] = - radius*ang_sin;
  crv2[3] = radius*ang_cos;
  crv2[4] = - radius*ang_cos;
  crv2[5] = - radius*ang_sin;
  crv2[6] = radius*ang_sin;
  crv2[7] = - radius*ang_cos;
  
  diffvec[0] = crv2[0] - crv1[0];
  diffvec[1] = crv2[1] - crv1[1];
  
  /* Calculate the position on the surface. */
  
  eder[0] = diffvec[1]*crv1[2] - diffvec[0]*crv1[3];
  eder[1] = diffvec[1]*crv2[2] - diffvec[0]*crv2[3];
  
  /* Calculate the first derivatives. */
  
  if (ider > 0)
    {
       eder[2] = diffvec[1]*crv1[4] - diffvec[0]*crv1[5];
       eder[3] = crv1[2]*crv2[3] - crv1[3]*crv2[2];
       eder[4] = crv1[2]*crv2[3] - crv1[3]*crv2[2];
       eder[5] = diffvec[1]*crv2[4] - diffvec[0]*crv2[5];
    }
  
  /* Calculate the second derivatives. */
  
  if (ider > 1)
    {
       eder[6] = crv1[2]*crv1[5] - crv1[3]*crv1[4]
	  + diffvec[1]*crv1[6] - diffvec[0]*crv1[7];
       eder[7] = crv1[4]*crv2[3] - crv1[5]*crv2[2];
       eder[8] = crv1[4]*crv2[3] - crv1[5]*crv2[2];
       eder[9] = crv1[2]*crv2[5] - crv1[3]*crv2[4];
       eder[10] = crv1[2]*crv2[5] - crv1[3]*crv2[4];
       eder[11] = crv2[3]*crv2[4] - crv2[2]*crv2[5]
	  + diffvec[1]*crv2[6] - diffvec[0]*crv2[7];
    }
	  
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("eval_crv_arc",*jstat,kpos);
  goto out;
  
  /* Dimension unequal to 2. */
  
  err102: *jstat = -102;
  s6err("eval_crv_arc",*jstat,kpos);
  goto out;

  /* More than 2 derivatives. */
  
  err103: *jstat = -103;
  s6err("eval_crv_arc",*jstat,kpos);

  out:
  
  return;
}

