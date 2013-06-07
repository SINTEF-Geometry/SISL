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
 * $Id: sh6getnext.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETNEXT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt* 
      sh6getnext(SISLIntpt *pt,int index)
#else
SISLIntpt* sh6getnext(pt,index)
   SISLIntpt *pt;
   int       index;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt and an index, fetch the next point
*              given by index.
*              If error, return SISL_NULL.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              index    - Index of link at pt.
*
*
* OUTPUT     : 
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{

   SISLIntpt *nextpt = SISL_NULL;

   /* check if index is within range */

   if(pt != SISL_NULL &&
      index >= 0 &&
      index < pt->no_of_curves) nextpt = pt->pnext[index];

   goto out;

   
   out :
      return nextpt;
}
