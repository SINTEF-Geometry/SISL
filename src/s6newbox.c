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
 * $Id: s6newbox.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6NEWBOX

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s6newbox(SISLbox *pbox,int inum,int itype,
	      double aepsge,int *jstat)
#else
void s6newbox(pbox,inum,itype,aepsge,jstat)
     SISLbox *pbox;
     int    inum;
     int    itype;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Create a particular box exist within an existing
*              box instance.
*
*
*
* INPUT      : pbox   - Box to modify.
*              inum   - Number of elements in min- and max-arrays.
*              itype  - Kind of box to create.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*              aepsge - Geometry resolution.
*
* OUTPUT     :  jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
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
* WRITTEN BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/                                     
{
   int knum = (inum == 1) ? inum : 2*inum;  /* If the geometry space has
					       dimension larger than 1,
					       a double set of min- and
					       max-arrays is to be made. */

   if (itype < 0 || itype > 2) goto err126;
   
   /* Test no such box exist, create the necessary arrays.  */
   
   if (pbox->e2min[itype] == SISL_NULL)
   {
      if ((pbox->e2min[itype] = newarray(knum,DOUBLE)) == SISL_NULL) goto err101;
      if ((pbox->e2max[itype] = newarray(knum,DOUBLE)) == SISL_NULL) goto err101;
   }
  
   /* Set the tolerance. */
   
   if (itype != 0) pbox->etol[itype] = aepsge;
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in input.  Kind of box do not exist.  */
   
   err126 : *jstat = -126;
   goto out;
   
   out :
      return;
}
