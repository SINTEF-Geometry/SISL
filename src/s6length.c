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
 * $Id: s6length.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */
#define S6LENGTH

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6length(double e1[],int idim,int *jstat)
#else
double s6length(e1,idim,jstat)
     double e1[]; 
     int    idim;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the length of a vector
*
* INPUT      : e1      - The vector 
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : jstat   - Status message
*                         0 - The length of the vector is zero
*                         1 - The length of the vector is not zero
*              s6norm  - The actual length of the vector
*
* METHOD     : The length of the input vector is calulated. 
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-june 1988
*                                  
*********************************************************************
*/
{
  register int ki;            /* Running variable in loop */
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
    goto mes01;

  /* Length of vector is zero    */

  *jstat = 0;
  goto out;

  /* Length of vector different from zero   */

 mes01: *jstat = 1;
        goto out;

 out: return(tsum);
}
