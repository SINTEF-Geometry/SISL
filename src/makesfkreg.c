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
 * $Id: makesfkreg.c,v 1.6 1994-11-30 12:53:02 pfu Exp $
 *
 */


#define MAKE_SF_KREG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void make_sf_kreg (SISLSurf * ps, SISLSurf ** rsnew, int *jstat)
#else
void
   make_sf_kreg (ps, rsnew, jstat)
     SISLSurf *ps;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a surface to a k-regular basis.
*
*
*
* INPUT      : ps	- Surface to be made k-regular.
*
*
*
* OUTPUT     : rsnew	- The new surface on a k-regular basis.
*              jstat	- status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
*
* WRITTEN BY : Ulf J. Krystad, SI, 04.92.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-08. Added error propagation.
*
**********************************************************************/
{
  int kn1=ps->in1;	/* Number of vertices in 1. par. dir.  */
  int kn2=ps->in2;	/* Number of vertices in 2. par. dir.  */
  int kk1=ps->ik1;	/* Order in 1. par. dir.               */
  int kk2=ps->ik2;	/* Order in 2. par. dir.               */
  /* --------------------------------------------------------- */

  s1001 (ps, ps->et1[kk1-1], ps->et2[kk2-1],
		ps->et1[kn1], ps->et2[kn2], rsnew, jstat);
  if (*jstat < 0)  goto error;

  goto out;

  /* Error in lower level routine */
error:
  s6err ("make_sf_kreg", *jstat, 0);

out:;

}
