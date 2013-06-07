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
 * $Id: s1707.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1707

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1707(SISLCurve *pc,int *jstat)
#else
void s1707(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if a B-spline curve is correct.
*
* INPUT      : pc     - SISLCurve to treat.
*
* OUTPUT     : jstat     - status messages
*                        > 0      : warning (1&2 only used when 
*                                            cuopen=SISL_CRV_PERIODIC)
*                                      = 1: Cyclic but not full freedom.
*                                      = 2: Not cyclic.
*                                      = 8: Non-positive rational weights.
*                       = 0      : ok
*                       < 0      : error
*
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Christophe Rene Birkeland, SI-SINTEF, May 1993.
*
**********************************************************************/
{

  int kpos=0;              /* Position of error. */
  int kstat=0;
  int step = 0;
  register double *s1,*s2; /* Pointers used in loop. */
  
  if (!pc) goto err150;

  if (pc->ik > pc->in) goto err111;
  
  if (pc->ik <= 0) goto err110;
  
  if (pc->in <= 0) goto err159;
  
  if (pc->idim <= 0) goto err102;
  
  if (pc->et[pc->in+pc->ik-1] <= *pc->et) goto err112;
  
  for (s1=pc->et,s2=pc->et+pc->in+pc->ik-1; s1<s2; s1++)
    if (s1[1] < *s1) goto err112;

  /* Check rational coefficients */
  if(pc->ikind == 2 || pc->ikind == 4)
    {
      step = pc->idim + 1;
      for (s1 = pc->rcoef + pc->idim, s2 = pc->rcoef + pc->in*step; 
	   s1 < s2; 
	   s1+= step)
	if (*s1 <= 0) goto war08;
    }

  /* Check if curve really is cyclic */
  if(pc->cuopen == SISL_CRV_PERIODIC)
    {
      test_cyclic_knots(pc->et,pc->in,pc->ik,&kstat);
      if (kstat < 0) goto error;
      if (kstat == 0) goto war02;
      if (kstat == 1) goto war01;
    }
      

  /* Updating output. No errors ! */
  
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
  
  /* Dimension less than 1. */
  
  err102:
    *jstat = -102;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order less than 1. */
  
  err110:
    *jstat = -110;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order greater than number of vertices. */
   
  err111:
    *jstat = -111;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Error in knotvector. */
  
  err112:
    *jstat = -112;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Null pointer. */
  
  err150:
    *jstat = -150;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Number of vertices less than 1. */
  
  err159:
    *jstat = -159;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine */
      
  error:
    *jstat = kstat;
    s6err("s1707",*jstat,kpos);
    goto out;

  out: 
    return;
}
