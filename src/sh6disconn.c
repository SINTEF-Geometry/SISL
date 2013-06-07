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
 * $Id: sh6disconn.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6DISCONNECT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6disconnect(SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6disconnect(pt1,pt2,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove the connection between pt1 and pt2 (which
*              must be unique) if any.
*
* INPUT      : pt1      - First point.
*              pt2      - Second point.
*              jstat    - Error flag.
*                         jstat =  0  => Successful.
*                         jstat =  1  => pt1 and pt2 are not connected
*                         jstat = -1  => Error in data structure.
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
  int kstat;                 /* Local status variable.            */
  int index1,index2;         /* Indices for pt1 and pt2.          */
  
   *jstat = 0;
  

   /* Check if pt1 and pt2 are connected. */

   sh6getlist(pt1,pt2,&index1,&index2,&kstat);
   if(kstat < 0) goto err1;
   if(kstat == 1)
   {
       *jstat = 1;
       goto out;
   }


   /* Disconnect. */

   pt1->no_of_curves--;
   pt1->pnext[index1] = pt1->pnext[pt1->no_of_curves];
   pt1->curve_dir[index1] = pt1->curve_dir[pt1->no_of_curves];

   pt2->no_of_curves--;
   pt2->pnext[index2] = pt2->pnext[pt2->no_of_curves];
   pt2->curve_dir[index2] = pt2->curve_dir[pt2->no_of_curves];

  
  goto out;  
  

  
  /* No connection exists. */

  err1 : *jstat = -1;
  s6err("sh6disconnect",*jstat,0);
  goto out;                       
  
  out: ;
}
