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
 * $Id: s1399.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1399

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1399(SISLCurve *pc,double astart,double astop)
#else
void s1399(pc,astart,astop)
     SISLCurve  *pc;
     double astart;
     double astop;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Change the knotvector to go from astart to astop.
*
* INPUT      : pc      - The curve.
*              astart  - Parametervalue at new startpoint.
*              astop   - Parametervalue at new endpoint.
*
*-
* CALLS      :
*
* WRITTEN BY : Morten Daehlen, SI, 88-09.
*
********************************************************************/
{
  int  kk= pc->ik;             /* Order of the input curve.             */
  int  kn= pc->in;             /* Number of vertices in the input curve.*/
  double *st=SISL_NULL;             /* Pointers used in loop.                */ 
  double a,b;
  int ii, kpos=0, kstat=0;
  if (!pc) goto out;
  
  if((st = newarray(kk+kn,DOUBLE)) == SISL_NULL) goto err101;
  
  a = pc -> et[kk-1];
  b = pc -> et[kn];
  
  for (ii=0;ii<kn+kk;ii++)
    st[ii] = astart+(((pc -> et[ii])-a)/(b-a))*(astop-astart);
  for (ii=0;ii<kn+kk;ii++)
    pc -> et[ii] = st[ii];
  goto out;
  
  /* Error in scratch allocation */
  
  err101: 
    kstat = -101;
    s6err("s1399",kstat,kpos);
    goto out;
      
  out:
    if (st != SISL_NULL) freearray(st);
    return;
}
