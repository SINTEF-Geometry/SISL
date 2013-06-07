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
 * $Id: s1232.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1232

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1232(double et1[],int in,int ik,
	   double afak1,double afak2,double et2[],int *jstat)
#else
void s1232(et1,in,ik,afak1,afak2,et2,jstat)
     double et1[];
     int    in;
     int    ik;
     double afak1;
     double afak2;
     double et2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Produce a knot vector that is stretched in the start
*              and/or end compared to an input knot vector.
*
*
*
* INPUT      : et1    - Original knot vector.
*              in     - Number of vertices of the curve which the
*                       knot vector belongs to.
*              ik     - Order of the curve which the knot vector
*                       belongs to.
*              afak1  - How much the new knot vector is to be stretched
*                       at the start compared with the original one. 
*                       afak1 >= 0. The parameter interval is extended 
*                       afak1 times the length of the parameter interval
*                       in the start of the interval.
*              afak2  - How much the new knot vector is to be stretched
*                       at the end compared with the original one. 
*                       afak2 >= 0.
*
*
*
* OUTPUT     : et2    - The stretched knot vector.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Replace the ik first and last knots with knots at the
*              new endpoints of the knot vector. The interior knots are
*              left unchanged.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/                   
{
  int kpos = 0;   /* Position of error.                */
  int ki;         /* Counter.                          */
  double tleng;   /* Length of parameter interval.     */
  double tstart;  /* New start value of knot vector.   */
  double tend;    /* New end value of knot vector.     */
  
  /* Test input.  */
  
  if (ik < 1) goto err110;
  if (in < ik) goto err111;
  
  /* Test if knot vector degenerated.  */
  
  tleng = et1[in] - et1[ik-1];
  if (tleng <= (double)0.0) goto err112;  
  
  /* Copy input knot vector to output knot vector.   */
  
  memcopy(et2,et1,in+ik,double);
  
  if (afak1 > (double)0.0)
    {
      
      /* Extend basis at start of curve.  */
      
      tstart = et1[ik-1] - tleng*afak1;
      for (ki=0; ki<ik; ki++) et2[ki] = tstart;
    }
  
  if (afak2 > (double)0.0)
    {
      
      /* Extend basis at end of curve.  */
      
      tend = et1[in] + tleng*afak2;
      for (ki=in; ki<in+ik; ki++) et2[ki] = tend;
    }
  
  /* New knot vector produced.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Order less than 1.  */
  
 err110: *jstat = -110;
  s6err("s1232",*jstat,kpos);
  goto out;
  
  /* Error in input. Number of vertices less than order.  */
  
 err111: *jstat = -111;
  s6err("s1232",*jstat,kpos);
  goto out;
  
  /* Error in input. Knot vector degenerated.  */
  
 err112: *jstat = -112;
  s6err("s1232",*jstat,kpos);
  goto out;
  
 out: return;
}
