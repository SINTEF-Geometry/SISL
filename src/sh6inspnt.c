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
 * $Id: sh6inspnt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6INSERTPT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6insertpt (SISLIntpt * pt1, SISLIntpt * pt2, SISLIntpt * ptnew, int *jstat)
#else
void 
sh6insertpt (pt1, pt2, ptnew, jstat)
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     SISLIntpt *ptnew;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Place a new intersection point between two (connected)
*              points.
*              pt1 and pt2 must be linked.
*
* INPUT      : pt1      - First point.
*              pt2      - Second point.
*              ptnew    - Point to be placed between pt1 and pt2.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*              UJK, Connect before disconnect, to
*              keep the index order
*********************************************************************
*/
{
  int kstat;			/* Local status variable.                  */
  int index1=0,index2=0;
  int crv_dir1=0,crv_dir2=0;
  
  *jstat = 0;

  sh6getlist (pt1, pt2, &index1, &index2, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */
  if (kstat == 1)
    goto err1;			/* pt1 and pt2 are not linked. */

  /* Save info in curve_dir */
  crv_dir1 = pt1->curve_dir[index1];
  crv_dir2 = pt2->curve_dir[index2];


  /* Check pt1,pt2,ptnew. */

  sh6connect (pt1, ptnew, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */

  /* Set values in curve_dir */
  sh6getlist (pt1, ptnew, &index1, &index2, &kstat);
  pt1->curve_dir[index1]   = crv_dir1;
  ptnew->curve_dir[index2] = crv_dir2;

  sh6connect (pt2, ptnew, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */

  /* Set values in curve_dir */
  sh6getlist (pt2, ptnew, &index1, &index2, &kstat);
  pt2->curve_dir[index1] = crv_dir2;
  ptnew->curve_dir[index2] = crv_dir1;


  sh6disconnect (pt1, pt2, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */
  if (kstat == 1)
    goto err1;			/* pt1 and pt2 are not linked. */


  goto out;


/* Error. pt1 and pt2 are not linked.  */

err1:*jstat = -1;
  s6err ("sh6insertpt", *jstat, 0);
  goto out;

/* Error in sub function.  */

error:*jstat = kstat;
  s6err ("sh6insertpt", *jstat, 0);
  goto out;

out:
  return;
}
