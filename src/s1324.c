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
 * $Id: s1324.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1324

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1324(double ecentr[],double aradiu,double enorm[],int idim,
	   double carray[],int *jstat)
#else
void s1324(ecentr,aradiu,enorm,idim,carray,jstat)
     double ecentr[];
     double aradiu;
     double enorm[];
     int    idim;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make two matrix of dimension 4x4
*              describing a 3-D circle as two implicit functions.
*
*
* INPUT      : ecentr - Center of the circle
*              aradiu - Radius of the circle
*              enorm  - Normal vector of circle plane
*              idim   - The dimension of the space the cirle lies in
*
*
*
* OUTPUT     : carray - The description of the circle. Outside
*                       this function the space for this array must be
*                       allocated. The need is 32 double variables.
*                       First the matrix for the sphere is stored,
*                       then the matrix of the plane.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The circle is described as an intersection between a
*              cylinder and the plane. The matrix describing the
*              cylinder is put first in the output array, the matrix
*              describing the plane follows then.
*              
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 29-June-1988
*                                    
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kdimp2;         /* idim + 2                                      */
  int kstop;          /* Stop condition for for loop                   */
  int ki;             /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  int kstat;          /* Status variable                               */
  
  
  
  /* Test i legal input */
  if (idim != 3) goto err104;
  
  kdimp1 = idim + 1;
  kdimp2 = idim + 2;
  kstop  = 2*kdimp1*kdimp1;
  
  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = (double)0.0;
    }
  
  /* Make description of cylinder */
  
  s1322(ecentr,enorm,aradiu,idim,1,carray,&kstat);
  if (kstat<0) goto error;
  
  
  /* Make description of plane, element (1,4), (2,4) and (3,4) */
  
  carray[28] = enorm[0];
  carray[29] = enorm[1];
  carray[30] = enorm[2];
  
  /* Make element (4,4) */
  
  carray[31] = -s6scpr(enorm,ecentr,idim);
  
  *jstat = 0;
  goto out;
  
  /* Dimension not 3 */
 err104: *jstat = -104;
  s6err("s1324",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine */
 error: *jstat = kstat;
  goto out;
  
  
 out:
  return;
}
