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
 * $Id: s6chpar.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CHPAR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6chpar(double ecoef1[],int in1,int in2,int idim,double ecoef2[])
#else
void s6chpar(ecoef1,in1,in2,idim,ecoef2)
     double ecoef1[];
     int    in1;
     int    in2;
     int    idim;
     double ecoef2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Change parameter directions of vertices of surface.
*
*
*
* INPUT      : ecoef1 - Vertices of original surface.
*              in1    - Number of vertices in first parameter direction.
*              in2    - Number of vertices in second parameter direction.
*              idim   - Dimension of the space in which the surfac lies.
*
*
* OUTPUT     : ecoef2 - Vertices after the changing of parameter directions.
*
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  register int ki,kj,kk;  /* Counters.  */
  
  for (ki=0; ki<in1; ki++)
    for (kj=0; kj<in2; kj++)
      for (kk=0; kk<idim; kk++)
	ecoef2[(ki*in2+kj)*idim+kk] = ecoef1[(kj*in1+ki)*idim+kk];
}
