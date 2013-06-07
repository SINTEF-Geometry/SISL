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
 * $Id: s1934.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1934

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1934 (double *et, int in, int ik, double start, double end, int *jstat)
#else
void
s1934 (et, in, ik, start, end, jstat)
     double *et;
     int in;
     int ik;
     double start;
     double end;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To map the support of a knot vector into a specified
*		interval.
*
*
* INPUT      :	et	- The original knot vector
*		in	- The number of degrees of freedom in the
*			  B-basis given by the knot vector.
*		ik	- The order of the basis.
*		start	- Start of the interval into which the knot
*			  vector is to be mapped.
*		end	- End of the interval into which the knot
*			  vector is to be mapped.
*
*
* OUTPUT     :	et	- The changed knot vector.
*               jstat   - Output status:
*                          < 0: Error.
*                          = 0: Ok.
*                          > 0: Warning.
*
* METHOD     :
*
* REFERENCES :  Fortran version:
*               T.Dokken, SI, 1981-10
*
* CALLS      : s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland
* REWISED BY :  Vibeke Skytt, SI, 92-10. The output knot vector will
*                                        have k-tupple knots in the ends.
*
*********************************************************************
*/
{
  int ii, stop;			/* Loop control parameters 		*/
  int kpos = 0;			/* Error position indicator		*/
  double store1;		/* Parameters used to decrease execution
				 * time					*/
  double fac;			/* Factor used in computation of new
				 * knot vector				*/

  *jstat = 0;


  /* Test if legal input */

  if ((ik < 1) || (in <ik))
    goto err112;

  if (start == end)
    goto err124;


  /* Perform normalization */

  store1 = et[ik - 1];
  fac = (end - start) / (et[in] -store1);
  stop = in +ik;

  for (ii=0; ii<ik; ii++) et[ii] = start;
  
  for (ii = ik; ii < in; ii++)
    et[ii] = fac * (et[ii] - store1) + start;

  for (ii = in; ii < stop; ii++) et[ii] = end;


  /* Normalization performed */

  goto out;


  /* Error in description of B-spline */

err112:
  *jstat = -112;
  s6err ("s1934", *jstat, kpos);
  goto out;

  /* The parameter interval is of zero length */

err124:
  *jstat = -124;
  s6err ("s1934", *jstat, kpos);
  goto out;

out:
  return;
}
