/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"



#define SH6CONDIR

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void sh6condir(SISLIntdat *rintdat, SISLIntpt *pt1, SISLIntpt *pt2, int *jstat)
#else
  void sh6condir(rintdat, pt1, pt2, jstat)
     SISLIntdat *rintdat;
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Connect the two points, make sure that they are main points
*              and set direction.
*
*
* INPUT      : rintdat  - Pool of intersection point
*              pt1      - Pointer to first Intpt.
*              pt2      - Pointer to second Intpt.
*              jstat    - Error flag.
*                        jstat < 0: Error
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 06-2023
*********************************************************************
*/
{
  int kstat = 0;
  *jstat = 0;
  sh6tomain (pt1, &kstat);
  sh6tomain (pt2, &kstat);
	      
  sh6idcon (&rintdat, &pt1, &pt2, &kstat);
  if (kstat < 0)
    *jstat = kstat;
  else
    sh6setdir (pt1, pt2, &kstat);
  if (kstat < 0)
    *jstat = kstat;
}
