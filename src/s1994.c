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
 * $Id: s1994.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S1994

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1994(SISLSurf *s1,int *jstat)
#else
void s1994(s1,jstat)
     SISLSurf *s1;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-surface intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : s1    - Surface in the intersection problem.
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
* WRITTEN BY : TDO SI, 89-08.
*
*********************************************************************
*/
{
  register int ki,kj;
  int kk1, kk2, kn1, kn2;
  int kbez;
  
  double tmaxt, tmaxs;
  double tmint, tmins;
  double tdiff;
  double *scoef=SISL_NULL;
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmaxt = tmaxs = - HUGE;
  tmint = tmins =   HUGE;
  
  /* Get surface attributes. */
  kk1  = s1->ik1;
  kk2  = s1->ik2;
  kn1  = s1->in1;
  kn2  = s1->in2;
  kbez = (kk1 == kn1) && (kk2 == kn2); 
  
  
  /* If the surface is linear in some direction it is simpel case. */
  if ((kk1 == 2 && kn1 == 2) || (kk2 == 2 && kn2 == 2)) goto out;
  
  
  /* Run through vertices in first parameter direction to find
     intervall of first derivative. */
  
  for (kj = 0, scoef = s1->ecoef;kj < kn2; kj++,scoef++)
    for (ki = 1; ki < kn1; ki++,scoef++ )
      {
	tdiff = *(scoef + 1) - *scoef;
	tmint = min(tmint,tdiff);
	tmaxt = max(tmaxt,tdiff);
      }
  
  /* Run through vertices in second parameter direction to find
     intervall of first derivative. */
  for (ki = 0 ;ki < kn1; ki++)
    for (kj = 1 , scoef = s1->ecoef + ki; kj < kn2; kj++, scoef +=kn1 )
      {
	tdiff = *(scoef + kn1) - *scoef;
	tmins = min(tmins,tdiff);
	tmaxs = max(tmaxs,tdiff);
      }

  /* ALA and UJK 30.10.90, remove noice near by zero */
  
  if (fabs(tmint) < noice) tmint = DZERO; 
  if (fabs(tmaxt) < noice) tmaxt = DZERO; 
  if (fabs(tmins) < noice) tmins = DZERO; 
  if (fabs(tmaxs) < noice) tmaxs = DZERO; 


  
  
  /* The first derivatives decide directions of possible intersection curves. */
  if (kbez && (tmint*tmaxt >=DZERO || tmins*tmaxs >=DZERO)) 
    *jstat = 1;
  else if (tmint*tmaxt > DZERO || tmins*tmaxs > DZERO) 
    *jstat = 1;
  else if (tmint == tmaxt  || tmins == tmaxs) 
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;
  
  goto out;
 out: ;
}

