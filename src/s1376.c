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
 * $Id: s1376.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1376

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1376(double et[],int in,int ik,double **gt,int *jkn,int *jkk,int *jstat)           
#else
void s1376(et,in,ik,gt,jkn,jkk,jstat)
     double et[];
     int    in;
     int    ik;
     double **gt;
     int    *jkn;
     int    *jkk;
     int    *jstat;           
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make the knot vector for the representing a spline
*              basis of order 4*(ik-1)+1, with the same knot values as et.
*              This basis is used for representing a curve or surface
*              put into a conic equation.
*
* INPUT      : et     - Knots of input spline basis
*              in     - Number of vertices in input basis
*              ik     - Order of input basis
*
* OUTPUT     : gt     - Pointer to array of knots. The array is allocated
*                       inside this routine.
*              jkn    - Number of vertices
*              jkk    - Order of B-spline basis produced
*
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
* METHOD     : 
*
* REFERENCES :
*
*-                                   
* CALLS      : 
*
* WRITTEN BY : Tor Dokken, SI, 88-11.
*
*********************************************************************
*/                                                               
{                                                                     
  double tval;     /* Value of knot                                 */
  double *sdum;    /* Pointer to knot array                         */
  int ki,kl;       /* Variable in loop                              */
  int knumb;       /* Number of intervals                           */
  int kstop;       /* Loop stop variable                            */
  int kpos=0;      /* Position of error                             */
  
  /* Run through the knot vector to decide how many intervals exist */
  
  knumb = 0;       
  tval = et[ik-1];
  
  for (ki=ik ; ki<=in ; ki++)
    {
      if (tval < et[ki])
        {
	  /*      New knot value found */
	  knumb = knumb + 1;
	  tval = et[ki];
        }
    }
  
  *jkk = 4*(ik-1) + 1;
  *jkn = (*jkk-1)*(knumb-1) + *jkk;
  
  sdum = newarray(*jkn+*jkk,DOUBLE);
  if (sdum == SISL_NULL) goto err101;
  
  *gt  = sdum; 
  
  /* Make knot values */
  
  tval = et[ik-1];
  
  /* Make jkk first knot values */
  
  for (kl=0;kl<*jkk;kl++)
    {
      sdum[kl] = tval;
    }
  
  /* kl points to the array entry where the next knot value is to be stored
   */
  
  for (ki=ik ; ki<=in ; ki++)
    {
      if (tval < et[ki])
        {
	  /* New knot value, remember this and make knots */
	  tval = et[ki];
	  kstop = kl + *jkk-1;
	  for (;kl<kstop;kl++)
            sdum[kl] = tval;
        }   
    }
  
  /* Make last knot value */
  
  sdum[kl] = tval;
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1376",*jstat,kpos);
  goto out;
 out:
  
  return;
}                                               
