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
 * $Id: s6idcpt.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDCPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idcpt(SISLIntdat *pintdat,SISLIntpt *pintpt,SISLIntpt **rintpt)
#else
void s6idcpt(pintdat,pintpt,rintpt)
     SISLIntdat *pintdat;
     SISLIntpt  *pintpt;
     SISLIntpt  **rintpt;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To find the point which is closest to pintpt
*              in the parametric space. If pintpt is the only
*              point in pintdat *rintpt is SISL_NULL.
*
*
*
* INPUT       :pintpt   - Pointer to an intersection point.
*              pintdat  - Pointer to intersection data.
*
*
* OUTPUT     : rintpt   - Pointer to a pointer to a point closest
*                         to pintpt.
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
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  if (pintdat == SISL_NULL)
    *rintpt = SISL_NULL;
  else
    {
      int ki,knr;                /* Counters.          */
      double tdist,td;           /* To store distanse. */
      
      if (pintpt == pintdat->vpoint[0])
        tdist = HUGE;
      else
        tdist = s6dist(pintdat->vpoint[0]->epar,pintpt->epar,pintpt->ipar);
      
      for (knr=0,ki=1; ki<pintdat->ipoint; ki++)
        {
	  if (pintpt == pintdat->vpoint[ki])
	    td = HUGE;
	  else
	    td = s6dist(pintdat->vpoint[ki]->epar,pintpt->epar,pintpt->ipar);
	  
	  if (td < tdist)
	    {
	      knr = ki;
	      tdist = td;
	    }
        }
      
      if (tdist == HUGE)
        *rintpt = SISL_NULL;
      else
        *rintpt = pintdat->vpoint[knr];
    }
}

