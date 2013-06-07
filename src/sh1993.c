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
 * $Id: sh1993.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */

#define SH1993

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1993(SISLCurve *c1,double aepsge,int *jstat)
#else
void sh1993(c1,aepsge,jstat)
     SISLCurve *c1;
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-curve intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : c1    - SISLCurve in the intersection problem.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                         = 1      : simpel case.
*                                         = 0      : not simpel case.
*                                         < 0      : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-06.
* Revised by : ALA and UJK 01.11.90, totale rewritten, same strategy
*              as in s1994 (surface case).
* REWISED BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
  register int ki,kj;

  int kk,kn;
  int kbez;
  double tmax;
  double tmin;
  double tdiff;
  double *scoef=SISL_NULL;
  /* ----------------------------------------------------------- */
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmax = - HUGE;
  tmin =   HUGE;
  
  /* Get curve attributes. */
  kk  = c1->ik;
  kn  = c1->in;
  kbez = (kk == kn);
  
  /* Run through vertices to find
     intervall of first derivative. */
  
  for (tdiff=DZERO, ki=1, scoef=c1->ecoef; ki<kn; ki+=kj, scoef+=kj)
  {
     for (kj=1; ki+kj<=kn; kj++)
     {
	if (tdiff*(*(scoef+kj) - *(scoef+kj-1)) < DZERO)
	   {
	      scoef += (kj-1);
	      ki += (kj-1);
	      kj = 1;
	   }
	   tdiff = *(scoef + kj) - *scoef;
	   if (fabs(tdiff) >= aepsge) break;
     }
     if (ki+kj > kn) break;
     
     tmin = min(tmin,tdiff);
     tmax = max(tmax,tdiff);
  }
  
  
  /* Simple case when no genuin zero's of first derivative. */
  if (kbez && (tmin*tmax >=DZERO)) 
    *jstat = 1;
  else if (tmin*tmax > DZERO) 
    *jstat = 1;
  else if (tmin == tmax)
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;

}



