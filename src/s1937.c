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
 * $Id: s1937.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1937

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1937 (double et[], int iordr, int ref, int left, double alfa[], double etref[])
#else
void
s1937 (et, iordr, ref, left, alfa, etref)
     double et[];
     int iordr;
     int ref;
     int left;
     double alfa[];
     double etref[];

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To calculate the discrete B-spline values number iref
*	   of the refinement etref of the array et.
*
*
* INPUT:   et	 - The original knot vector.
*	   iordr - The order of the discrete B-spline to be
*		   calculated.
*	   ref	 - The index of the discrete B-spline to be
*		   calculated.
*	   left	 - Arrayindex, satisfying:
*		   et[left-1] <= etref[ref-1] < et[left]
*	   etref - Refined knot vector.
*
*
* OUTPUT:  alfa	 - The values of the calculated discrete B-splines.
*		   alfa[0]    - Corresponds to number left-iordr+1.
*		   alfa[1]    - Corresponds to number left-iordr+2.
*		   alfa[left] - Corresponds to number left.
*
* METHOD: We use the Oslo-algorithm developped by Cohen, Lyche and
*         Riesenfeld.
*
* REFERENCES: Cohen, Lyche, Riesenfeld: Discrete B-splines and subdivision
*	      techniques in computer aided geometric design, computer
*	      graphics and image processing, vol 14, no.2 (1980)
*
* CALLS: No.
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07
*
*********************************************************************
*/
{
  int ki, kl, kr, low;		/* Loop control parameters. 	*/
  int stop, start;		/* and array indicators.	*/
  double tj, td1, td2;		/* Parameters used to improve.	*/
  double beta1, beta;		/* algorithm.			*/


  /* We have et[left-1] <= etref[ref-1] < et[left]
     So the discrete B-splines can be calculated. */

  low = left - iordr;
  start = left - 1;
  stop = iordr - 1;
  alfa[stop] = 1;

  for (kr = 0; kr < stop; kr++)
    {
      beta1 = (double) 0.0;
      tj = etref[ref + kr];
      if (start < 0)
	start = 0;

      for (ki = start; ki < left; ki++)
	{
	  kl = ki - low;
	  td1 = tj - et[ki];
	  td2 = et[ki + kr + 1] - tj;
	  beta = alfa[kl] / (td1 + td2);
	  alfa[kl - 1] = td2 * beta + beta1;
	  beta1 = td1 * beta;
	}
      alfa[iordr - 1] = beta1;
      start--;
    }
  return;
}
