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
 * $Id: eval_2_crv.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define EVAL_2_CRV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   eval_2_crv(SISLCurve *pc1,SISLCurve *pc2,int ider,double epar[],int *ilfs,
	   int *ilft,double eder[],int *jstat)
#else
void eval_2_crv(pc1,pc2,ider,epar,ilfs,ilft,eder,jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     int       ider;
     double    epar[];
     int       *ilfs;
     int       *ilft;
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
*              in epar. Compute ider derivatives.
*
*
*
* INPUT      : pc1    - Pointer to the first curve.
*            : pc2    - Pointer to the second curve.
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
*                ilft  - Corresponding to ilfs in the second parameter
*                        direction.
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
  
  /* Check the input. */
  
  if (pc1->idim != 2 || pc2->idim != 2) goto err102;
  if (ider > 2) goto err103;
  
  /* Evaluate the curves. */
  
  s1221(pc1, kder, epar[0], ilfs, crv1, &kstat);
  if (kstat < 0) goto error;

  s1221(pc2, kder, epar[1], ilft, crv2, &kstat);
  if (kstat < 0) goto error;
  
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
  s6err("eval_2_crv",*jstat,kpos);
  goto out;
  
  /* Dimension unequal to 2. */
  
  err102: *jstat = -102;
  s6err("eval_2_crv",*jstat,kpos);

  /* More than 2 derivatives. */
  
  err103: *jstat = -103;
  s6err("eval_2_crv",*jstat,kpos);

  out:
  
  return;
}

