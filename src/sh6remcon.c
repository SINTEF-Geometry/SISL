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
 * $Id: sh6remcon.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6REMCON

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6remcon(SISLIntdat **pintdat,SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6remcon(pintdat,pt1,pt2,jstat)
   SISLIntdat **pintdat;
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove connection between pt1 and pt2.
*
* INPUT      : pintdat  - Pointer to pointer to the SISLIntdat data.
*              pt1      - First point.
*              pt2      - Second point.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat = -2  => error in space allocation
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
* CALLS      : s6err      - Gives error message.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
  int kstat;                 /* Local status variable.                     */
  
   *jstat = 0;
  

  sh6disconnect(pt1,pt2,&kstat);
  if (kstat < 0) goto error;

  /* Do we want to remove either point from pintdat? */

  
  goto out;  
  

  /* Error in lower level routine. */

  error : *jstat = kstat;
  s6err("sh6remcon",*jstat,0);
  goto out;                       
  
  
  out: ;
}
