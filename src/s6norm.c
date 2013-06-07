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
 * $Id: s6norm.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */
#define S6NORM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6norm(double e1[],int idim,double e2[],int *jstat)
#else
double s6norm(e1,idim,e2,jstat)
     double e1[];
     int    idim;
     double e2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To normalize a vector to length 1
*
* INPUT      : e1      - The vector to be normalized
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : jstat   - Status message
*                         0 - The length of the vector is zero
*                         1 - The length of the vector is one
*              s6norm  - The actual length of the vector
*              e2      - The nomalized version of the vector
*
* METHOD     : The length of the input vector is calulated, and the
*              output is assigned the values of the input vector and 
*              divided by the length.
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*                                  
*********************************************************************
*/
{
  register int ki;      /* Running variable in loop */
  register double tsum=DZERO; /* Dummy variables in summing loop */
  
  /* If the dimension is 1 the length of the vector is the same as the
   *  absolute value of the number */
  
  if (idim == 1)
    tsum = fabs(e1[0]);
  else
    {
      for (ki=0;ki<idim;ki++)
	tsum += (e1[ki]*e1[ki]);
      
      tsum = sqrt(tsum);
    }
  
  if (DNEQUAL(tsum,DZERO))
    for (ki=0;ki<idim;ki++)
      e2[ki] = e1[ki]/tsum;
  else
    {
      for (ki=0;ki<idim;ki++)
        e2[ki] = DZERO;
      goto mes00;
    }
  goto mes01;

  /* Length of vector is zero    */

 mes00: 
  *jstat = 0;
       goto out;

  /* Length of vector different from zero   */

 mes01: 
  *jstat = 1;
  goto out;
out: return(tsum);
}
