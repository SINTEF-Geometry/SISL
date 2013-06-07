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
 * $Id: s6idnpt.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDNPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idnpt(SISLIntdat **pintdat,SISLIntpt **pintpt,int itest,int *jstat)
#else
void s6idnpt(pintdat,pintpt,itest,jstat)
     SISLIntdat **pintdat;
     SISLIntpt  **pintpt;
     int    itest;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To insert a new intersection point into pintdat.
*              If pintdat is SISL_NULL a new pintdat is also made.
*              If pintpt is close to an other intersection point
*              the object pintpt is pointing to is freed, and
*              pintpt is set to point to the already inserted point.
*
*
*
* INPUT      : pintpt   - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection data.
*              itest    - Indikate testing equalety.
*                               = 1      : Testing.
*                               = 0      : No testing.
*
*
* OUTPUT     : jstat  - status messages  
*                               = 2      : Already existing.
*                               = 1      : Already inserted.
*                               = 0      : Intersection point inserted.
*                               < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              newIntdat  - Create new intdat structure.
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  register int ki,kj;              /* Counters.    */
  
  /* We have to be sure that we have an intdat structure. */
  
  if ((*pintdat) == SISL_NULL)
    {
      if (((*pintdat) = newIntdat()) == SISL_NULL) goto err101;
    }
  
  
  /* Than we have to be sure that we do not have the intersection point
     before or an equal point. */
  
  for (ki=0; ki<(*pintdat)->ipoint; ki++)
    if ((*pintdat)->vpoint[ki] == (*pintpt))
      {
	*jstat = 1;
	goto out;
      }
    else if (itest && (*pintpt)->iinter != 2)
      {
	for (kj=0; kj<(*pintpt)->ipar; kj++)
	  if (DNEQUAL((*pintpt)->epar[kj],(*pintdat)->vpoint[ki]->epar[kj]))
	    break;
	
	if (kj == (*pintpt)->ipar)
	  {
	    freeIntpt(*pintpt);
	    (*pintpt) = (*pintdat)->vpoint[ki];
	    *jstat = 2;
	    goto out;
	  }
      }
  
  
  /* Than we have to be sure that the array vpoint is great enought. */
  
  if (ki == (*pintdat)->ipmax)
    {
      (*pintdat)->ipmax += 20;
      
      if (((*pintdat)->vpoint = increasearray((*pintdat)->vpoint,
					      (*pintdat)->ipmax,SISLIntpt *)) == SISL_NULL) 
	goto err101;
    }
  
  
  /* Now we can insert the new point. */
  
  (*pintdat)->vpoint[ki] = (*pintpt);
  (*pintdat)->ipoint++;
  *jstat = 0;
  goto out;
  

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6idnpt",*jstat,0);
        goto out;

 out: ;
}
