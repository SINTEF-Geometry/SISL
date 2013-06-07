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
 * $Id: refine_all.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define REFINE_ALL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
refine_all (SISLIntdat ** pintdat,
	    SISLObject * po1,
	    SISLObject * po2,
	    double eimpli[],
	    int ideg,
	    double aepsge,
	    int *jstat)

#else
void
refine_all (pintdat,
	    po1,
	    po2,
	    eimpli,
	    ideg,
	    aepsge,
	    jstat)

     SISLIntdat **pintdat;
     SISLObject *po1;
     SISLObject *po2;
     double eimpli[];
     int ideg;
     double aepsge;
     int *jstat;

#endif
/*
*********************************************************************
*
* PURPOSE    : An empty function, (to ensure similarity with other
*              versions).
*
*
* INPUT      : pintdat     - Pointer to pointer to the SISLIntdat data.
*              po1         - Pointer surface object.
*              po2         - Pointer surface object.
*              eimpli      - Array containing descr. of implicit surf
*	       ideg        - Type of impl surf.
              ang_tol     - Angle control tolerance ie ??
*              aepsge      - Absolute tolerance
*
*
* OUTPUT     :  jstat  - status messages
*                       = ?      : ?
*                       = 0      : ok
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway, July-1990
*
*********************************************************************
*/
{
  *jstat = 0;
}
