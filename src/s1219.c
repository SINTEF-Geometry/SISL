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
 * $Id: s1219.c,v 1.2 2000-05-23 08:44:05 vsk Exp $
 *
 */


#define S1219

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1219(double *et,int ik,int in,int *ileft,double ax,int *jstat)
#else
void s1219(et,ik,in,ileft,ax,jstat)
     double *et;
     int    ik;
     int    in;
     int    *ileft;
     double ax;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To localize the point ax in the array et.
*              The output ileft should satisfy the relations
*                          
*                    et[ileft] <= ax < et[ileft+1].
* 
*              There are two exceptions to this. If (ax >= et[in])
*              then ileft should be in-1 (this corresponds to extending
*              the polynomial piece between et[in-1] and et[in] to the
*              right of the natural parameter interval.
*              Similarly, if (ax < et[ik-1]) then ileft should still be
*              ik-1.
*
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax     - The point at which the B-spline values and derivatives
*                       are to be computed.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                       where ax is located, check the relations above.
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The aim is to do as little work as possible in the cases
*              where ileft has the right or almost the right value.
*              First of all we make sure that ileft has a legal value
*              (a value in the range ik-1 to in-1). Then we check
*              if the current value is OK.
*              If it is not we check that ax is in the interior of et
*              or if the right value is obtained by either increasing
*              or decreasing ileft by 1. If the right value still has
*              not been found we do a binary search.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
*
*********************************************************************
*/                                     
{
  int kpos=0;         /* The position of the error.                      */
  int kleft;          /* Local version of ileft to avoid the pointer.    */
  
  /* Check the input. */
  
  if (ik < 1) goto err110;
  
  if (in < ik) goto err111;
  
  if (et[ik-1] == et[ik] || et[in-1] == et[in]) goto err112;
  
  /* Make sure that kleft is in the legal range. */
  
  kleft = min(max(ik-1,*ileft),in-1);
  
  /* Check if the current value of kleft is acceptable. */
  
  if (et[kleft] <= ax && ax < et[kleft+1]) ;
  
  /* Check if ax is outside (et[ik-1],et[in]). */
  
  else if (ax >= et[in-1])
    kleft = in - 1;
  else if (ax <= et[ik-1])
    kleft = ik - 1;
  
  /* Check if it is sufficient to increase or decrease kleft by one. */
  
  else if (et[kleft+1] <= ax && ax < et[kleft+2])
    kleft += 1;
  else if (kleft > 0 && et[kleft-1] <= ax && ax < et[kleft])
    kleft -= 1;
  
  /* Last resort - a binary search. */
  else
    {
      
      /* kmin and kmax gives the upper and lower limits on the possible values
	 of kleft.                                                       */
      
      int kmin,kmax;
      
      kmin = ik - 1; kmax = in - 1;
      kleft = (kmin+kmax)/2;
      
      while (ax < et[kleft] || et[kleft+1] <= ax)
	{
	  if (ax < et[kleft])
	    kmax = kleft;
	  else
	    kmin = kleft;
	  
	  kleft = (kmin+kmax)/2;
	}
    }
  
  *ileft = kleft;
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  /* Polynomial order less than 1. */
 err110: *jstat = -110;
  s6err("s1219",*jstat,kpos);
  goto out;
  
  /* Fewer B-splines than the order. */
 err111: *jstat = -111;
  s6err("s1219",*jstat,kpos);
  goto out;
  
  /* Error in knot vector.
     (The first or last interval of the knot vector is empty.) */
 err112: *jstat = -112;
  s6err("s1219",*jstat,kpos);
  goto out;
  
 out: return;
}
