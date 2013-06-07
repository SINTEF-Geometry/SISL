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
 * $Id: s6idkpt.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDKPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idkpt(SISLIntdat **pintdat,SISLIntpt **pintpt,SISLIntpt **rtpt,SISLIntpt **rfpt,int *jstat)
#else
void s6idkpt(pintdat,pintpt,rtpt,rfpt,jstat)
     SISLIntdat **pintdat;
     SISLIntpt  **pintpt;
     SISLIntpt  **rtpt;
     SISLIntpt  **rfpt;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To remove an intersection point pintpt from pintdat.
*              Pintpt is than killed. If there is a point pointing on
*              pintdat ptpt is set to point to this point otherwise
*              ptpt is SISL_NULL. If pintdat is pointing to a point pfpt is
*              set to point to this point otherwise pfpt is SISL_NULL.
*              If pintdat is empty pintdat is killed and set to SISL_NULL.
*              pintpt is set to SISL_NULL.
*
*
*
* INPUT/OUTPUT:pintpt   - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT     : rtpt     - Pointer to a pointer to a point pointing to pintpt.
*              rfpt     - Pointer to a pointer to a point pintpt points to.
*              jstat    - status messages  
*                               = 1      : Pintpt is not in pintdat.
*                               = 0      : OK!
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
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  int ki;              /* Counters.    */
  int knum;
  
  (*rtpt) = (*rfpt) = SISL_NULL;
  *jstat = 0;
  
  /* We have to be sure that we have an intdat structure. */
  
  if ((*pintdat) == SISL_NULL)
    goto out;
  
  if ((*pintpt) == SISL_NULL)
    {
      *jstat = 1;
      goto out;
    }
  
  
  /* Than we have to be sure that we do not have the intersection point
     before or an equal point. */
  
  for (knum = -1,ki=0; ki<(*pintdat)->ipoint; ki++)
    {
      if ((*pintdat)->vpoint[ki] == (*pintpt))
	knum = ki;
      
      if ((*pintdat)->vpoint[ki] == (*pintpt)->pcurve)
	(*rfpt) = (*pintdat)->vpoint[ki];
      
      if ((*pintdat)->vpoint[ki]->pcurve == (*pintpt))
	(*rtpt) = (*pintdat)->vpoint[ki];
    }
  
  
  if (knum == -1)
    *jstat = 1;
  else
    {
      (*pintdat)->vpoint[knum] = (*pintdat)->vpoint[(*pintdat)->ipoint-1];
      ((*pintdat)->ipoint)--;
      (*pintdat)->vpoint[(*pintdat)->ipoint] = SISL_NULL;
      
      if ((*rtpt) != SISL_NULL)
	(*rtpt) ->pcurve = SISL_NULL;
      
      if ((*pintdat)->ipoint == 0)
	{
	  freeIntdat(*pintdat);
	  (*pintdat) = SISL_NULL;
	}
    }
  
  freeIntpt(*pintpt);
  (*pintpt) = SISL_NULL;
  
 out: ;
}
