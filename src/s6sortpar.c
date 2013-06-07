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
 * $Id: s6sortpar.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6SORTPAR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      s6sortpar(double evec1[],double epar1[],int ipar,int idim,
		double evec2[],double epar2[],int *jstat)
#else
void s6sortpar(evec1,epar1,ipar,idim,evec2,epar2,jstat)
     double evec1[];
     double epar1[];
     int ipar;
     int idim;
     double evec2[];
     double epar2[];
     int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : Sort vector of parameter values and a corresponding vector
*              after increasing parameter value. 
* 
* 
* INPUT      : evec1    - Vector corresponding to vector of parameter 
*                         values. To each parameter value idim elements
*                         correspond. Dimension of array is ipar*idim.
*              epar1    - Parameter values to be sorted. Dimension is ipar.
*              ipar     - Number of parameter value.
*              idim     - Dimension of geometry space.
*
* 
* OUTPUT     : evec2    - Vector sorted according to corresponding parameter
*                         values. Dimension is ipar*idim. This array may
*                         be the same as evec1.
*              epar2    - Sorted parameter value. Dimension is ipar. This
*                         array may be the same as epar1.
*              jstat    - status messages 
*                          > 0 : warning 
*                          = 0 : ok 
*                          < 0 : error 
* 
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 08.90.
*
*********************************************************************
*/
{
  int ki,kj;            /* Counters.  */
  double tpar;          /* Intermediate storage of double.  */
  double *svec = SISL_NULL;  /* Intermediate storage of double vector.  */
  
  /* Allocate scratch for local array.  */

  if ((svec = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Move input arrays to output arrays.  */

  memcopy(evec2,evec1,idim*ipar,DOUBLE);
  memcopy(epar2,epar1,ipar,DOUBLE);
  
  /* Sort output arrays according to parameter value.  */

  for (ki=0; ki<ipar-1; ki++)
    for (kj=ki+1; kj<ipar; kj++)
      if (epar2[kj] < epar2[ki])
	{
	  /* Interchange parameter values.  */

	  tpar = epar2[kj];
	  epar2[kj] = epar2[ki];
	  epar2[ki] = tpar;

	  /* Interchange elements of corresponding vector.  */

	  memcopy(svec,evec2+kj*idim,idim,DOUBLE);
	  memcopy(evec2+kj*idim,evec2+ki*idim,idim,DOUBLE);
	  memcopy(evec2+ki*idim,svec,idim,DOUBLE);
	}
  
  /* Arrays sorted.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

 err101:
  *jstat = -101;
  goto out;
  
 out:
  
  /* Free space occupied by local array.  */

  if (svec != SISL_NULL) freearray(svec);
  
  return;
}
