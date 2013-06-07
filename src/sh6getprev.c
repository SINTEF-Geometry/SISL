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
 * $Id: sh6getprev.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETPREV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6getprev(SISLIntpt *pt1,SISLIntpt *pt2)
#else
int sh6getprev(pt1,pt2)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt pt1 and a pointer to another Intpt pt2,
*              fetch the index of the pt1 array corresponding
*              to pt2. If no such index exists return -1.
*
*
* INPUT      : pt1       - Pointer to the Intpt.
*              pt2     - Pointer to another Intpt.
*
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. May 91.
*
*********************************************************************
*/
{
   int       ncurv;   /* number of curves pt1 is connected to       */
   int       index;   /* index number for pnext array              */

   index = -1;

   if(pt1 == SISL_NULL || pt2 == SISL_NULL) goto out;

   ncurv = pt1->no_of_curves;  /* note ncurv can be zero */

   index=0;
   while(index < ncurv && pt1->pnext[index] != pt2) index++;
   if(index == ncurv) index = -1;  /* no index found */

   goto out;

   out :
      return index;
}
