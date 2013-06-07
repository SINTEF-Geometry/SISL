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
 * $Id: sh6idsplit.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define S6IDSPLIT


#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    sh6idsplit (SISLIntdat ** pintdat, SISLIntpt * psource, int *jstat)
#else
void
   sh6idsplit (pintdat, psource, jstat)
     SISLIntdat **pintdat;
     SISLIntpt *psource;
     int *jstat;
#endif


/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To split an intersection point into several instanses so
*              that the each of the instances has exactly one (main) neighbour.
*              
*
*
* INPUT:       psource  - Pointer to an intersection point.
* 
* 
* INPUT/OUTPUT:pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT  :    jstat    - status messages
*                               = 0      : OK
*                               < 0      : error
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
*
* WRITTEN BY : Ulf J. Krystad, 04.92.
*
*********************************************************************
*/
{
  int ki;			/* Counters.                         */
  int no_main;			/* No of neighbours (main points)    */
  int test= FALSE;              /* No equality testing when inserted
				   in pintdat                        */
  int kstat = 0;                /* Local status.                     */
  SISLIntpt *pneighb = SISL_NULL;	/* Current neighbour                 */
  SISLIntpt *pshadow = SISL_NULL;	/* Current copy of source point      */
  /* ------------------------------------------------*/
  
  *jstat = 0;
  
  if (psource == SISL_NULL)
    {
       *jstat = 1;
       goto out;
    }
  
  /* Get number of neighbours */
  no_main = sh6nmbmain (psource, &kstat);
  if (kstat < 0)
    goto error;
  
  for (ki=psource->no_of_curves - 1; no_main > 1; ki--)
    {
       pneighb = sh6getnext(psource, ki);
       if (!pneighb) goto error;
       if (sh6ismain(pneighb))
	 {
	    pshadow = hp_copyIntpt(psource);
	    sh6idnpt(pintdat, &pshadow, test=FALSE, &kstat);
	    if (kstat < 0) goto error;
	    
	    sh6insertpt(psource, pneighb, pshadow, &kstat);
	    if (kstat < 0) goto error;
	    
	    sh6disconnect(psource, pshadow, &kstat);
	    if (kstat < 0) goto error;
	    no_main--;
	 }
    }
  goto out;
  
  
error:
  *jstat = kstat;
  goto out;

out:;
}
