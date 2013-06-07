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
 * $Id: s1235.c,v 1.3 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1235

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1235(double et[],int in,int ik,int *jnbreak,double **gbreak,int *jstat)
#else
void s1235(et,in,ik,jnbreak,gbreak,jstat)
     double et[];
     int    in;
     int    ik;
     int    *jnbreak;
     double **gbreak;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find break points in a knot vector. The first and last
*              parameter values are break values.
*
*
*
* INPUT      : et     - Knot vector to find break points in.
*              in     - Number of vertices of the curve corresponding
*                       to the knot vector.
*              ik     - Order of the curve corresponding to et.
*
*
*
* OUTPUT     : jnbreak - Number of break points found.
*              gbreak  - Array containing parameter values of break points.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The knot vector has a break point at a knot if the 
*              multiplicity of the knot is ik-1.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Vibeke Skytt, SINTEF, 9801. Correction in loop counter
*                                          for the periodic case.
*
*********************************************************************
*/
{
  int kpos = 0;   /* Position of error.                    */
  int kj;         /* Counter.                              */
  int kbreak;     /* Current number of break points found. */
  int kmult;      /* Multiplisity of current knot.         */
  double tprev;   /* Value of previous knot.               */
  double *sbreak; /* Pointer into break point array.       */
  double *st;     /* Pointer into knot vector.             */
  
  /* Allocate space for an array that is at least as great as the
     number of break points.                                       */
  
  *gbreak = SISL_NULL;
  if ((*gbreak = newarray(in+2,double)) == SISL_NULL) goto err101;
  
  /* Set local pointer to and counter of break points.  */
  
  sbreak = *gbreak;
  kbreak = 0;
  
  /* Find break point in start of parameter interval and internal breaks. */
  
  tprev = et[ik-1];
  kmult = ik - 1;
  for (st=et+ik,kj=ik; kj<in; st++,kj++)
    {
      
      if (*st == tprev) kmult++;
      else
	{
	  if (kmult >= ik-1)
	    {
	      
	      /* New break point found.  */
	      
	      *(sbreak++) = tprev;
	      kbreak++;
	    }
	  tprev = *st;
	  kmult = 1;
	}
    }
  
  /* Find break point in end of interval.  */
  
  if (et[in] != tprev && kmult >= ik-1)
    {
      
      /* Remember last internal break point.  */
      
      *(sbreak++) = tprev;
      kbreak++;
    }
  *(sbreak++) = et[in];
  kbreak++;
  
  /* Reduce break point array to correct size.  */
  
  if (kbreak < in+2)
    if ((*gbreak = increasearray(*gbreak,kbreak,double)) == SISL_NULL) goto err101;
  
  /* Break points found.  */
  
  *jnbreak = kbreak;
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: 
  *jstat = -101;
  s6err("s1235",*jstat,kpos);
  goto out;
  
 out: 
  return;
}
