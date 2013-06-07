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
 * $Id: s1236.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1236

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1236(double et[],int in,int ik,int inpar,double epar[],int *jstat)
#else
void s1236(et,in,ik,inpar,epar,jstat)
     double et[];
     int    in;
     int    ik;
     int    inpar;
     double epar[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To find the constant parameter values to be used in
*              surface/curve drawing.
*
*
*
* INPUT      : et     - Knot vector.
*              in     - Number of vertices corresponding to et.
*              ik     - Order corresponding to et.
*              inpar  - Number of constant parameter values wanted.
*
*
*
* OUTPUT     : epar   - Array containing constant parameter values.
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
* CALLS      : s1235 - Calculate break points of the knot vector.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;   /* Local status variable.                             */
  int kpos = 0;    /* Position of error.                                 */
  int ki,kj;       /* Counters.                                          */
  int kpar;        /* Remaining number of parameter values to be set
		      at non-break points.                               */
  int knumb;       /* Number of parameter values between two break points.*/
  int knbreak = 0; /* Number of break points.                            */
  double tend;     /* Value of the end of the parameter interval.        */
  double tlast;    /* Value of last break point.                         */
  double tval;     /* Value of current break point.                      */
  double tdist;    /* Distance between two parameter values.             */
  double *spar;    /* Pointer used to traverse epar.                     */
  double *sbreak = SISL_NULL;  /* Array containing break points.              */
  
  /* Test input.  */
  
  if (ik < 1) goto err110;
  if (in < ik) goto err111;
  
  /* Find break points, including endpoints of the parameter interval. */
  
  s1235(et,in,ik,&knbreak,&sbreak,&kstat);
  if (kstat < 0) goto error;
  
  /* Find number of constant parameter values not at break points. */
  
  kpar = inpar - knbreak;
  
  /* Adjust number of constant parameter values if the number of break point
     is greater than the number of parameter values wanted.  */

  if (kpar < 0)
    {
      sbreak[1] = sbreak[knbreak-1];
      knbreak = 2;
      kpar = inpar - 2;
    }
  
  /* Set first parameter value and end of parameter interval.  */
  
  spar = epar;
  tlast = *spar = sbreak[0];
  tend = sbreak[knbreak-1];
  
  for (spar++,ki=1; ki<knbreak; ki++)
    {
      
      /* Set constant parameter values between two break points.  */
      
      tval = sbreak[ki];
      tdist = tval - tlast;
      knumb = (int)(kpar*tdist/(tend-tlast));
      
      /* Remaining non-break point parameter values to be set.  */
      
      kpar -= knumb;
      
      /* Distance between parameter values. */
      
      tdist /= (knumb + 1);
      
      /* Set parameter values.  */
      
      for (kj=0; kj<knumb; kj++) *(spar++) = tlast + (kj+1)*tdist;
      
      /* Set parameter value at break point. */
      
      *(spar++) = tval;
      tlast = tval;
    }
  
  /* Parameter values found. */
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Order less than one.  */
  
 err110: *jstat = -110;
  s6err("s1236",*jstat,kpos);
  goto out;
  
  /* Error in input. Number of vertices less than order.  */
  
 err111: *jstat = -111;
  s6err("s1236",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1236",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local array.  */
  
  if (sbreak != SISL_NULL) freearray(sbreak);
  
  return;
}
