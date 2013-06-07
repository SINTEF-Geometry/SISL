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
 * $Id: sh6iscnect.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6ISCONNECT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6isconnect(SISLIntpt *pt0, SISLIntpt *pt1, SISLIntpt *pt2)
#else
int sh6isconnect(pt0, pt1, pt2)
    SISLIntpt *pt0;
    SISLIntpt *pt1;
    SISLIntpt *pt2;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Boolean. Test if there exist a connection between
*              the two intersection points pt1 and pt2.
*
*
* INPUT      : pt0       - Pointer to previous Intpt, possibly SISL_NULL.
*              pt1       - Pointer to first Intpt.
*              pt2       - Pointer to second Intpt.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, 93-02.
*
*********************************************************************
*/
{
   int kstat = 0;   /* Status on wether a connection is found.  */
   int ki;          /* Counter.                                 */
   SISLIntpt *qt;   /* Intersection point.                      */
   int been_here = -199;

   /* Test if the points are equal.  */

   if (pt1 == pt2) return 1;
   
   /* UJK, aug 93, oo-loop in sh6isconn. */
   if (pt1->marker == been_here) return 0;
   pt1->marker = been_here;

   /* Traverse the intersection points connected to pt1.  */
   
   for (ki=0; ki<pt1->no_of_curves; ki++)
   {
      qt = sh6getnext(pt1,ki);
      if (qt == pt0) continue;
      
      kstat = sh6isconnect(pt1,qt,pt2);
      
      if (kstat) return 1;
   }
   
   /* No connection is found.  */
   
   return 0;
}
