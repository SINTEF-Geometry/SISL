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

#define S1708

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1708(SISLSurf *ps,int *jstat)
#else
void s1708(ps,jstat)
     SISLSurf *ps;
     int   *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if a B-spline surface is correct.
*
* INPUT      : ps      - SISLSurf to treat.
*
* OUTPUT     : jstat   - status messages
*                        > 0      : warning (1&2 only used when 
*                                            cuopen==SISL_CRV_PERIODIC in at
*                                            least one direction)
*                                      = 1: Cyclic but not full freedom.
*                                      = 2: Not cyclic.
*                                      = 8: Non-positive rational weights.
*                       = 0      : ok
*                       < 0      : error
*
* CALLS      :
*
* WRITTEN BY : Christophe Rene Birkeland, SI-SINTEF, May 1993.
*
**********************************************************************/
{
  int type = 1;
  int stat = 0;              /* Status for lower level routines   */
  int kstat1= 0;
  int kstat2= 0;
  int kpos = 2;              /* Indicator of parameter direction
			      * for error message.                */
  int step = 0;
  register double *s1,*s2;   /* Pointers used in loop.            */
  SISLCurve *curve=SISL_NULL;     /* Local curve                       */

  if(ps->ikind == 3) type = 3;


  /*
   *  Check second parameter direction 
   *  
   */

  curve = newCurve(ps->in2, ps->ik2, ps->et2, ps->ecoef,
		   type, ps->idim * ps->in1, 0);
  if(curve == SISL_NULL) goto err101;

  s1707( curve, &stat);
  if(stat<0) goto error;

  /* Free curve element used in 2. par. direction */

  freeCurve(curve);
  curve = SISL_NULL;
  

  /*
   *  Check first parameter direction
   *
   */

  kpos = 1;

  /*
   * NOT necessary to transpose coefficient array. 
   * Coefficients are not checked
   * in routine s1707 for curve-kind 1 || 3
   *
   */

  curve = newCurve(ps->in1, ps->ik1, ps->et1, ps->ecoef,
		   type, ps->idim * ps->in2, 0);
  if (curve == SISL_NULL) goto err101;

  s1707( curve, &stat);
  if(stat<0) goto error;

  kpos = 0;

  /* Check weights in case of rational surface */
  if(ps->ikind == 2 || ps->ikind == 4)
    {
      step = ps->idim + 1;
      for (s1 = ps->rcoef + ps->idim, s2 = ps->rcoef + ps->in1*ps->in2*step; 
           s1 < s2;
           s1+= step)
        if (*s1 <= 0) goto war08;
    }

  /* Check if surface really is cyclic */
  if(ps->cuopen_1 == SISL_CRV_PERIODIC || ps->cuopen_2 == SISL_CRV_PERIODIC)
    {
      kpos = 1;
      test_cyclic_knots(ps->et1,ps->in1,ps->ik1,&kstat1);
      if (kstat1 < 0) goto error;
      kpos = 2;
      test_cyclic_knots(ps->et2,ps->in2,ps->ik2,&kstat2);
      if (kstat2 < 0) goto error;
      if ((kstat1 == 0 && ps->cuopen_1 == SISL_CRV_PERIODIC)
	  || (kstat2 == 0 && ps->cuopen_2 == SISL_CRV_PERIODIC))
	goto war02;
      if ((kstat1 == 1 && ps->cuopen_1 == SISL_CRV_PERIODIC)
	  || (kstat2 == 1 && ps->cuopen_2 == SISL_CRV_PERIODIC))
	goto war01;
    }
      
  /* SUCCESS. */

  *jstat = 0;
  goto out;

  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector does not give
   * full freedom. */
  
  war01:
    *jstat = 1;
    goto out;
  
  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector not cyclic. */
  
  war02:
    *jstat = 2;
    goto out;
  
  /* Warning: Non-positive rational coefficients. */
  
  war08:
    *jstat = 8;
    goto out;
  
  /* Error in allocations */

  err101:
    *jstat = -101;
    s6err("s1708", *jstat, kpos);
    goto out;

  /* Lower level error */

  error:
    *jstat = stat;
    s6err("s1708",*jstat,kpos);
    goto out;

 out: 
  if (curve != SISL_NULL) freeCurve( curve );
  return;
}
