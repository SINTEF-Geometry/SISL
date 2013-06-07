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
 * $Id: sh6tohelp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TOHELP

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6tohelp(SISLIntpt *pt,int *jstat)
#else
void sh6tohelp(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if pt is a mai point. If it is, transform to
*              help point, if not give a message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => pt was a main point, now help.
*                         jstat =  1  => pt is not a main point.
*                         jstat = -1  => Error in pt.
*                         jstat = -2  => Illegal to convert status
*                         jstat = -3  => Error in data structure.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
   int kstat; /* Local status */
   int num; 

   *jstat=0;

   if(pt == SISL_NULL) goto err1;

   if(sh6ismain(pt))  /* If pt is a help point. */
   {
      /* ??????????? */
      /* if(pt->no_of_curves > 2) goto err2; */
      
      num=sh6nmbmain(pt,&kstat);
      /* Problem in sh6edgred when starting reduction */
      /* if(num > 1) goto err2; */

       pt->iinter = -pt->iinter;  /* Convert status to main point. */
   }
   else
   {
       *jstat=1;
   }

   goto out;
   

err1:
   /* Error in input. pt is null. */
   
   *jstat = -1;
   s6err("sh6tohelp",*jstat,0);
   goto out;
   
   /* Error, Illegal to change status. */
   
   /* err2:
    *jstat = -2;
   s6err("sh6tohelp",*jstat,0);
   goto out; */
   
   
   out :
      return;
}
