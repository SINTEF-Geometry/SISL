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
 * $Id: s6schoen.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */



#define S6SCHOEN

#include "sislP.h"
#if defined(SISLNEEDPROTOTYPES)
double
s6schoen(double et[], int ik, int index)
#else
double s6schoen(et,ik,index)
     double et[];
     int    ik;
     int    index;
#endif

/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To determine the knot value of a specified vertice.
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              index  - The vertice index at where the knot values are to 
*                       be computed.
*
*                
*
* INPUT/OUTPUT :
*              s6schoen - The knot value at index.
*
*
* METHOD     : The aim is to calculate the knot value of a vertice, using the
*              Schoenberg spline expression:
*
*               *
*              t  = (t    + .............. +t     )/k-1.
*               i     i+1                    i+k-1
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Per Evensen,SI, August 1991.
*
*********************************************************************
*/                                     
{
  int i;             /* Loop variable                                   */
  double kval=DZERO; /* knot value variable                             */
  
  for (i=index+1;i<index+ik;i++) kval+=et[i];
  kval = kval/(ik-1);

return(kval);
}

