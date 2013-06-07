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
 * $Id: newknots.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define NEWKNOTS

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
  newknots(double et[],int in,int ik,double epar[],int inpar,double aeps,
	   double **ginsert,int *jinsert,int *jstat)
#else
void newknots(et,in,ik,epar,inpar,aeps,ginsert,jinsert,jstat)
   double et[];
   int in;
   int ik;
   double epar[];
   int inpar;
   double aeps;
   double **ginsert;
   int *jinsert;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To set up an array of parameter values to insert into
*              a B-spline curve when the curve is to have knot 
*              multiplisity equal to the order at specified parameter
*              values.
*
*
*
* INPUT      : et      - Knot vector of the B-spline curve.
*              in      - Number of vertices of the curve.
*              ik      - Order of the curve.
*              epar    - Parameter values at which the refined curve
*                        is going to have an order-multiple knot.
*              inpar   - Number of parameter values in epar.
*              aeps    - The smallest legal distance between new and
*                        current knots.
*
*
*
* OUTPUT     : ginsert - Array containing new knots to the curve.
*              jinsert - Number of knots in ginsert.
*              jstat   - status messages  
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
* CALLS      : s1219  - Find position of parameter value in the knot
*                       vector of a curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 11.90.
*
*********************************************************************
*/
{
   int kstat = 0;        /* Local status variable.               */
   int ki,kj;            /* Counters.                            */
   int knpar;            /* Number of new knots found currently. */
   int kleft = 0;        /* Position in knot vector.             */
   double tpar;          /* Value of current new knot.           */
   
   /* Allocate maximum needed scratch for output array.  */
   
   *jinsert = 0;
   if ((*ginsert = newarray(ik*inpar,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Traverse the array of parameter values and compute how
      may knots is to be inserted at each value.  */
   
   for (ki=0; ki<inpar; ki++)
   {
      tpar = epar[ki];
      
      /* Find position in knot vector.  */
      
      s1219(et,ik,in,&kleft,tpar,&kstat);
      if (kstat < 0) goto error;
      
      /* Check if the parameter value is close to a knot.  */
      
      if (tpar-et[kleft] < aeps)
      {
	 tpar = et[kleft];
	 for (kj=1, knpar=ik-1; 
	  kj<=kleft && DEQUAL(et[kleft],et[kleft-kj]);
	  kj++, knpar--);
      }
      else if (et[kleft+1]-tpar < aeps)
      {
	 tpar = et[kleft+1];
	 for (kj=2, knpar=ik-1;
	  kleft+kj<in+ik && DEQUAL(et[kleft+1],et[kleft+kj]);
	  kj++, knpar--);
      }
      else knpar = ik;
      
      /* Register the new knots.  */
      
      for (kj=0; kj<knpar; kj++)
	 (*ginsert)[*jinsert + kj] = tpar;
      (*jinsert) += knpar;
   }
   
   /* Set correct size of output array.  */
   
   if (*jinsert != ik*inpar)
   {
      if ((*ginsert = increasearray(*ginsert,MAX(1,*jinsert),DOUBLE)) == SISL_NULL)
	 goto err101;
   }
   
   /* New knots found.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in lower level routine.  */
    
   error : *jstat = kstat;
   goto out;
   
   out :
      return;
}
