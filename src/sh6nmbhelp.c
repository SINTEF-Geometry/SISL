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

/*
 *
 * $Id: sh6nmbhelp.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6NMBHELP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6nmbhelp(SISLIntpt *pt,int *jstat)
#else
int sh6nmbhelp(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt, return the number of help points
*              it is linked to.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*********************************************************************
*/
{
   int num; /* Number of lists. */
   int ki; /* Loop variable.  */

   num=0;

   /* Count number of main lists pt lies in. */

   for(ki=0; ki<pt->no_of_curves; ki++)
   {
       if(pt->pnext[ki] == SISL_NULL) goto err1;
       if(sh6ishelp(pt->pnext[ki])) num++;
   }

   goto out;
   

err1:
   /* Error in data structure. */
   
   *jstat = -1;
   s6err("sh6nmbhelp",*jstat,0);
   goto out;
   
   
   out :
      return num;
}






