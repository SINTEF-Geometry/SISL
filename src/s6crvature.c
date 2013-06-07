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
 * $Id: s6crvature.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CURVATURE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      s6curvature(double eder[],int idim,double ecurv[],int *jstat)
#else	 
void s6curvature(eder,idim,ecurv,jstat)
     int idim,*jstat;
     double eder[],ecurv[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given position, first and second derivative of a curve at
*              a point, compute curvature vector.
*
*
*
* INPUT      : eder    - Array containing position, 1. and 2. derivative of
*                        curve. Dimension is 3*idim.
*              idim    - Dimension of geometry space.
*                       
*
* OUTPUT     : ecurv   - Curvature vector.
*              jstat   - status messages  
*                                         = 1      : Tangent vector zero.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Express the curve using cord length parametrisation. Then the
*              curvature is equal to the 2. derivative of the curve.
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : s6length - Length of vector.   
*              s6scpr   - Scalar product between two vectors.  
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  int kstat = 0;    /* Status variable.  */
  int ki;           /* Counter.   */
  double tleng;     /* Length of first derivative.  */
  double tleng2;    /* Square of length.  */
  double tdot;      /* Scalar product between 1. and 2. derivative. */
  
  /* Compute length of 1. derivative. */

  tleng = s6length(eder+idim,idim,&kstat);
  
  if (kstat == 0)
    {
      /* The first derivative is zero. */

      for (ki=0; ki<idim; ki++) ecurv[ki] = (double)0.0;
      goto warn1;
    }
  else
    {
      tleng2 = tleng*tleng;
      tdot = s6scpr(eder+idim,eder+2*idim,idim);
      
      for (ki=0; ki<idim; ki++)
	ecurv[ki] = (eder[2*idim+ki] - eder[idim*ki]*tdot/tleng2)/tleng2;
    }
  
  *jstat = 0;
  goto out;
  
  warn1 :
    *jstat = 1;
  goto out;
  
  out :
    return;
}

  
