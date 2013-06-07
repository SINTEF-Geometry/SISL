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
 * $Id: s6existbox.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6EXISTBOX

#include "sislP.h" 

#if defined(SISLNEEDPROTOTYPES)
int s6existbox(SISLbox *pbox,int itype,double aepsge)
#else
int s6existbox(pbox,itype,aepsge)
     SISLbox *pbox;
     int    itype;
     double aepsge;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if a particular box exist within an existing
*              box instance.
*
*
*
* INPUT      : pbox   - Box to test.
*              itype  - Kind of box to test existance of.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*              aepsge - Geometry resolution.
*
* OUTPUT     : s6existbox -  Status.
*                            -1 : Kind of box exist, but is expanded
*                                 with another tolerance.
*                             0 : Kind of box do not exist.
*                             1 : Requested box exist.
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
   if (pbox->e2min[itype] == SISL_NULL) return(0);  /* No box is made. */
   
   if (itype != 0 && DNEQUAL(pbox->etol[itype],aepsge))
      return(-1);  /* Box exist, but with another size of the expansion. */
   
   return(1);
}
   
