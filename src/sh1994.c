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
 * $Id: sh1994.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH1994

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1994(SISLSurf *s1,double aepsge,int *jstat)
#else
void sh1994(s1,aepsge,jstat)
     SISLSurf *s1;
     double aepsge;
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
* INPUT      : s1     - Surface in the intersection problem.
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
* WRITTEN BY : TDO SI, 89-08.
* REWISED BY : Vibeke Skytt, SI, 91-02.
*              UJK, SI,91-10 Bug in first direction when a row
*                            ends with 2 or more equal numbers.
*                            Also added a test to include
*                            tmax < tmin (both HUGE)
*********************************************************************
*/
{
  register int ki,kj,kh;
  int kk1, kk2, kn1, kn2;
  int kbez;
  
  double tmaxt, tmaxs;
  double tmint, tmins;
  double tdiff;
  double *scoef=SISL_NULL;
  
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
  
  /* UJK, 91-10 */
  /* for (kj=0, scoef=s1->ecoef; kj<kn2; kj++,scoef++) */
  for (kj=0, scoef=s1->ecoef; kj<kn2; kj++,scoef=s1->ecoef+kn1*kj)
     for (tdiff=DZERO, ki=1; ki<kn1; ki+=kh, scoef+=kh)
     {
	for (kh=1; ki+kh<=kn1; kh++)
	{
	   if (tdiff*(*(scoef+kh) - *(scoef+kh-1)) < DZERO)
	      {
		 scoef += (kh-1);
		 ki += (kh-1);
		 kh = 1;
	      }
	      tdiff = *(scoef + kh) - *scoef;
	      if (fabs(tdiff) >= aepsge) break;
	}
	if (ki+kh > kn1) break;
	
	tmint = min(tmint,tdiff);
	tmaxt = max(tmaxt,tdiff);
     }
  
  /* Run through vertices in second parameter direction to find
     intervall of first derivative. */
  
  for (ki=0; ki<kn1; ki++)
     for (tdiff=DZERO, kj=1, scoef=s1->ecoef+ki; kj<kn2; kj+=kh, scoef+=kh*kn1)
     {
	for (kh=1; kj+kh<=kn2; kh++)
	{
	   if (tdiff*(*(scoef+kh*kn1) - *(scoef+(kh-1)*kn1)) < DZERO)
	      {
		 scoef += (kh-1)*kn1;
		 kj += (kh-1);
		 kh = 1;
	      }
	      tdiff = *(scoef + kh*kn1) - *scoef;
	      if (fabs(tdiff) >= aepsge) break;
	}
	if (kj+kh > kn2) break;
	
	tmins = min(tmins,tdiff);
	tmaxs = max(tmaxs,tdiff);
     }

  /* UJK, 91-10, maybe parameters not set */
  if (tmint > tmaxt || tmins > tmaxs)
  {
     *jstat = 1;
     goto out;
  }
  
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

