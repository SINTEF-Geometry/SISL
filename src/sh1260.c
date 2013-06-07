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

#define SH1260

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1260(double aconst,SISLCurve *vcurve[],int icurve,int *jstat)
#else	 
void sh1260(aconst,vcurve,icurve,jstat)
     double aconst;
     SISLCurve *vcurve[];
     int icurve;
     int *jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Check length of tangent vectors at the endpoints of a
*              curve compared to the size of the curve. If the vectors
*              are too long, reparametrize the curve and a number of
*              corresponding curves.
*
*
*
* INPUT      : aconst     - Constant used to check when the tangent
*                           vectors at the endpoints of the curve are
*                           too long.
*              icurve     - Number of curves. icurve >= 1.
*              
*
* INPUT/OUTPUT : vcurve   - Array containing a curve set. The curves
*                           are expected to have the same parametrization.
*                           If the tangents in the endpoints of the
*                           first curve is too long, all curves are
*                           reparametrized. Dimension of the array is icurve.
*                       
*
* OUTPUT     : jstat      - status messages  
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
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221    - Evaluate curve.  
*              s6diff   - Difference vector between two vectors.  
*              s6scpr   - Scalar product between two vectors.  
*              s6length - Length of vector.   
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int kstat = 0;      /* Status variable.   */
  int ki;             /* Counter.           */
  int kder = 1;       /* Number of derivatives of curve to evaluate. */
  int kleft = 0;      /* Parameter used in curve evaluation.         */
  int kdim = vcurve[0]->idim;   /* Dimension of first curve in the curve set.     */
  double tpar1 = *(vcurve[0]->et + vcurve[0]->ik - 1);   /* Start of parameter
                                                            interval of 1. curve. */
  double tpar2 = *(vcurve[0]->et + vcurve[0]->in);       /* End of parameter
                                                            interval of 1. curve. */
  double *sder1 = SISL_NULL;  /* Value of 1. curve in start of parameter interval.   */
  double *sder2 = SISL_NULL;  /* Value of 1. curve in end of parameter interval.     */
  double *sdiff = SISL_NULL;  /* Difference vector between endpoints of 1. curve.    */
  double tdiff;          /* Length of sdiff.  */
  double t1,t2;          /* Length of the components of the tangent vectors in the 
                            endpoints along sdiff, compared with the length of sdiff.*/
  double tscal;          /* The factor with which to scale the curves.  */
  double tnewend;        /* New endpoint of parameter interval.         */
  double *s1;            /* Pointer used to traverse knot vector of curve. */  
  double *s2;            /* Pointer used to stop traversing knot vector.   */
  SISLCurve *qcpt;           /* Pointer to curve in curve set. */
  
  /* Test input.  */

  if (icurve < 1) goto err110;
  for (ki=1; ki<icurve; ki++)
    {
      qcpt = vcurve[ki];
      if (qcpt->idim != kdim) goto err106;
      if (*(qcpt->et+qcpt->ik-1) != tpar1) goto err112;
      if (*(qcpt->et+qcpt->in) != tpar2) goto err112;
    }
      
  /* Allocate scratch for value and derivatives of the first curve
     in the endpoints.  */

  if ((sder1 = newarray(2*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((sder2 = newarray(2*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((sdiff = newarray(kdim,DOUBLE)) == SISL_NULL) goto err101;

  /* Evaluate the first curve in the endpoints.  */

  s1221(vcurve[0],kder,tpar1,&kleft,sder1,&kstat);
  if (kstat < 0) goto error;
  
  s1221(vcurve[0],kder,tpar2,&kleft,sder2,&kstat);
  if (kstat < 0) goto error;

  /* Compute difference vector between endpoints of position curve. */

  s6diff(sder2,sder1,kdim,sdiff);
  
  /* Compute length of difference vector.  */

  tdiff = s6length(sdiff,kdim,&kstat);
  
  /* Compute length of the component of the tangent in the first endpoint
     along the difference vector compared to the distance between the
     endpoints. The result lies between 0 and 1.  */

  t1 = s6scpr(sder1+kdim,sdiff,kdim)/(tdiff*tdiff);   

/*  t1 = s6length(sder1+kdim,kdim,&kstat);  */
  
  /* Compute length of the component of the tangent in the second endpoint
     along the difference vector compared to the distance between the
     endpoints. The result lies between 0 and 1.  */

  t2 = s6scpr(sder2+kdim,sdiff,kdim)/(tdiff*tdiff);  
  
/*  t2 = s6length(sder2+kdim,kdim,&kstat);   */
  
  /* Check if any of the tangents are too long.  */

  if (MAX(t1,t2) > aconst)
    {
      /* One of the tangents is too long. Reparametrize to reduce tangent
	 length.   */

      tscal = MAX(t1,t2)/aconst;
      
      /* Find new endpoint of parameter interval. The startpoint is kept. */

      tnewend = tscal*(tpar2-tpar1) + tpar1;

      for (ki=0; ki<3; ki++)
	{
	  qcpt = vcurve[ki];
	  
	  /* Traverse position and u- and v-derivative curves.  */

	  for (s1=qcpt->et,s2=qcpt->et+qcpt->in+qcpt->ik; s1<s2; s1++)
	    *s1 = tpar1 + (*s1 - tpar1)*(tnewend - tpar1)/(tpar2 - tpar1);
	}
      
    }
 
  /* Testing and evt. reparametrization performed.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Number of curve is less than one.  */

  err110 :
    *jstat = -110;
  goto out;
  
  /* Error in input. Conflicting dimensions.  */

  err106 :
    *jstat = -106;
  goto out;
  
  /* Error in input. The curves do not have the same parameter interval.  */

  err112 :
    *jstat = -112;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :
    
    /* Free space occupied by local arrays and curves.  */

    if (sder1 != SISL_NULL) freearray(sder1);
  if (sder2 != SISL_NULL) freearray(sder2);
  if (sdiff != SISL_NULL) freearray(sdiff);
  
  return;
}     
