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
 * $Id: s1993.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S1993

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1993(SISLCurve *c1,int *jstat)
#else
void s1993(c1,jstat)
     SISLCurve *c1;
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
*********************************************************************
*/
{
  register int ki;

  int kk,kn;
  int kbez;
  double tmax;
  double tmin;
  double tdiff;
  double *scoef=SISL_NULL;
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
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
  
  for (ki = 1,scoef= c1->ecoef;ki < kn; ki++,scoef++ )
      {
	tdiff = *(scoef + 1) - *scoef;
	tmin = min(tmin,tdiff);
	tmax = max(tmax,tdiff);
      }
  
  if (fabs(tmin) < noice) tmin = DZERO; 
  if (fabs(tmax) < noice) tmax = DZERO; 
  
  
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



