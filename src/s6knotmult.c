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
 * $Id: s6knotmult.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6KNOTMULT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s6knotmult(double et[],int ik,int in,int *ileft,double ax,int *jstat)
#else
int s6knotmult(et,ik,in,ileft,ax,jstat)
     double et[];
     int    in;
     int    ik;
     int    *ileft;
     double ax;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To determine the multiplicity of a knot at a specified
*              parameter value.
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
* INPUT/OUTPUT :
*              ileft - Pointer to the interval in the knot vector
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
* CALLS      : s1219
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
*
*********************************************************************
*/                                     
{
  int kpos=0;         /* The position of the error.                      */
  int kstat;          /* Local status variable                           */
  int kmult=0;        /* Multiplicity of knot                            */
  int ki;             /* Loop variable                                   */
  
  /* Localize knot interval */
  
  s1219(et,ik,in,ileft,ax,&kstat);
  if (kstat<0) goto error;
  
  if (et[*ileft] == ax)
    {
      kmult = 1;
      ki    = *ileft-1;
      for (ki=(*ileft)-1; 0 <= ki; ki--)
        if (et[ki] == ax) kmult++;
    }
  if (et[in] == ax)
    {
      for (ki=in ; ki<in+ik;ki++)
        if (et[ki] == ax) kmult++;
    }
  
  *jstat = 0;
  goto out;

/* Error in lower level function */

error:  *jstat = kstat;
        s6err("s6knotmult",*jstat,kpos);
        goto out;

out:
return(kmult);
}
