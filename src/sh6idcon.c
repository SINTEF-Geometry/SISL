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
 * $Id: sh6idcon.c,v 1.2 2001-03-19 15:59:08 afr Exp $
 *
 */


#define SH6IDCON

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined (SISLNEEDPROTOTYPES)
void
sh6idcon (SISLIntdat ** pintdat, SISLIntpt ** pintpt1, SISLIntpt ** pintpt2, int *jstat)
#else
void 
sh6idcon (pintdat, pintpt1, pintpt2, jstat)
     SISLIntdat **pintdat;
     SISLIntpt **pintpt1;
     SISLIntpt **pintpt2;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To connect two intersection points in pintdat into a list.
*              If pintdat is SISL_NULL a new pintdat is also made.
*              If  one of pintpt is close to an other intersection point
*              the object pintpt is pointing to is removed, and
*              pintpt is set to point to the already inserted point.
*              The direction of the connection is set to go 
*              from pintpt1 to pintpt2.
*
*
*
* INPUT      : pintpt1  - Pointer to a pointer to new intersection point.
*              pintpt2  - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection date.
*
*
* OUTPUT     : jstat  - status messages
*                               = 0      : Connection done.
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
*              sh6connect - Do the connection on low level.
*              sh6idnpt   - Insert a new intpt structure.
*
* WRITTEN BY : Ulf J. Krystad 06.91
* REVISED BY:  Micael Floater July 91. Set direction of connection.
* REVISED BY:   Ulf J. Krystad 07.91  Removed direction of connection.
*********************************************************************
*/
{
  int kstat;			/* Local status variable.                     */

  /* First we have to be sure that pintdat contain the two points. */

  sh6idnpt (pintdat, pintpt1, 1, &kstat);
  if (kstat < 0)
    goto error;

  sh6idnpt (pintdat, pintpt2, 1, &kstat);
  if (kstat < 0)
    goto error;

    /* Connect */
  sh6connect (*pintpt1, *pintpt2, &kstat);
  if (kstat < 0)
    goto error;

    /* Set direction of connection. */
    /*  sh6setdir(*pintpt1, *pintpt2, &kstat);
    if (kstat < 0)
    goto error; */


  *jstat = 0;
  goto out;

  /* Error from lower function */
error:
  *jstat = kstat;
  s6err ("sh6idcon", *jstat, 0);
  out:
     ;
}
