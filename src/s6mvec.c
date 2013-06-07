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
 * $Id: s6mvec.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6MVEC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6mvec(double emat[],double evec1[],int inbvec,double evec2[])
#else
void s6mvec(emat,evec1,inbvec,evec2)
     double emat[];
     double evec1[];
     int    inbvec;
     double evec2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To post multiply a matrix by inbvec 3-D vectors.
*
* INPUT      : emat    - The matrix (16 elements, Homogenous coordinates)
*              evec1   - The input vectors (3-D coordinates)
*              inbvec  - The number of input vectors
*
* OUTPUT     : evec2   - The resulting vectors (3-D coordinates)
*
*-  
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-24
*                                  
*********************************************************************
*/
{
  double svec2[3];           /* Temporary vector                        */
  double *svec1;             /* Pointer to current vector               */
  double *svec3;             /* pointer to current result vector        */
  register double tdum;      /* Temporary storage of real               */
  register int ki,kj,kl,kp;  /* Control variables in loops              */
  int kstop;                 /* Stop condition for loop                 */
  
  /* Multiply matrix by all vectors */
  
  kstop = 3*inbvec;
  
  for (kl=0;kl<kstop;kl=kl+3)
    {
      /* Multiply rotational part of the matrix by the vector */
      
      svec1 = evec1 + kl;
      svec3 = evec2 + kl;
      
      for (ki=0;ki<3;ki++)
        {
	  kp = ki;
	  
	  tdum = DZERO;
	  for (kj=0;kj<3;kj++)
            {
	      tdum += emat[kp]*svec1[kj];
	      kp += 4;
            }
	  
	  /* Add translation part */
	  svec2[ki]  = tdum + emat[kp];
        }
      /*  Check if the bottom row is 0,0,0,1 */                     
      
      if (DNEQUAL(emat[3],DZERO) || DNEQUAL(emat[7],DZERO) ||
	  DNEQUAL(emat[11],DZERO) || DNEQUAL(emat[15],(double)1.0))
        {
	  /* Compute last element of vector */
	  
	  tdum = evec1[0]*emat[3] + evec1[1]*emat[7] + evec1[2]*emat[11];
	  if (DNEQUAL(tdum,DZERO))
            {
	      for (ki=0;ki<3;ki++)
		svec2[ki] /= tdum;
            }
	}
      svec3[0] = svec2[0];    
      svec3[1] = svec2[1];
      svec3[2] = svec2[2];
    }
  
  return;
}
