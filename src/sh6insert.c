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
 * $Id: sh6insert.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6INSERT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6insert(SISLIntdat **pintdat,SISLIntpt *pt1,SISLIntpt *pt2,SISLIntpt **ptnew,int *jstat)
#else
void sh6insert(pintdat,pt1,pt2,ptnew,jstat)
   SISLIntdat **pintdat;
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   SISLIntpt **ptnew;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Place a new intersection point between two (connected)
*              points. 
*              NOTE that if there exists a point in pintdat which is
*              close to ptnew, then the input point is killed and
*              no insertion is done.
*              pt1 and pt2 must lie contiguously in a list.
*
* INPUT      : pintdat  - Pointer to pointer to the SISLIntdat data.
*              pt1      - First point.
*              pt2      - Second point.
*              ptnew    - Point to be placed between pt1 and pt2.
*              jstat    - Error flag.
*                         jstat =  1  => Warning, point existing in
*                                        pintdat, nothing done.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat = -2  => error in space allocation
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
* CALLS      : s6err      - Gives error message.
*              sh6idnpt    - Insert a new intpt structure.
*              copyIntpt  - Copy an intpt structure.
*              newIntdat  - Create new intdat structure.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* MODYFIED BY: UJK, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
  int kstat;                /* Local status variable.                     */
  
   *jstat = 0;
  
  /* First we have to be sure that pintdat contains ptnew. */
  
  sh6idnpt(pintdat,ptnew,1,&kstat);
  if (kstat < 0) goto error;
  if (kstat > 0) 
    { 
       /* Point already existing in data structure, point killed, 
	  no insertion */
       *jstat = 1;
       goto out;
    }
  
  /* UJK, aug. 92 insert always mainpts if one of the neighbour is a main */
  /* if (sh6ismain(pt1) && sh6ismain(pt2)) */
   if (sh6ismain(pt1) || sh6ismain(pt2))
     sh6tomain(*ptnew,&kstat);
  else
     sh6tohelp(*ptnew,&kstat);
  if (kstat < 0) goto error;

  /* Then insert the point. */
  sh6insertpt(pt1,pt2,*ptnew,&kstat);
  if (kstat < 0) goto error;
  
  
  goto out;
  

/* Error. pt1 and pt2 are not properly connected.  */


/* Error in sub function.  */

error:  *jstat = kstat;
        s6err("sh6insert",*jstat,0);
        goto out;

   out:
      return;
}
