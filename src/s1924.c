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
 * $Id: s1924.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1924

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1924 (int id1, int id2, int id3, int id4, int in1, int in2,
       double **ew, int *jstat)
#else
void
s1924 (id1, id2, id3, id4, in1, in2, ew, jstat)
     int id1;
     int id2;
     int id3;
     int id4;
     int in1;
     int in2;
     double **ew;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To compute the weight matrix to blend together two
*		surfaces.
*
*
* INPUT      :	id1	- Derivative on edge one up to (id1-1).
*		id2	- Derivative on edge two up to (id2-1).
*		id3	- Derivative on edge three up to (id3-1).
*		id4	- Derivative on edge four up to (id4-1).
*		in1	- Size of matrix in first direction
*		in2	- Size of matrix in second direction
*
*
* OUTPUT     :  ew	- Matrix containing the weights.
*            :  jstat   - output status:
*                         < 0: error.
*                         = 0: ok.
*                         > 0: warning.
*
*
* REFERENCES :  Fortran version:
*		Morten Daehlen, SI, 1987-01
*
* CALLS: s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-08
*
*********************************************************************
*/
{
  int ki, kj;			/* Loop control parameters 		*/
  int in1m1, in2m1;		/* Equals in1-1 and in2-1		*/
  int kpos = 0;

  *jstat = 0;


  /* Initialization */

  *ew = new0array (in1 * in2, double);
  if (*ew == SISL_NULL)
    goto err101;

  in1m1 = in1 - 1;
  in2m1 = in2 - 1;
  (*ew)[0] = (double) 0.5;
  (*ew)[in1m1] = (double) 0.5;
  (*ew)[in2m1 * in1] = (double) 0.5;
  (*ew)[in2m1 * in1 + in1m1] = (double) 0.5;

  for (ki = 1; ki < in1m1; ki++)
    {
      (*ew)[ki] = (double) 1.0;
      for (kj = 1; kj < id1; kj++)
	if (ki < id4 || ki >= in1 - id2)
	  (*ew)[kj * in1 + ki] = (double) 0.5;
	else
	  (*ew)[kj * in1 + ki] = (double) 1.0;

      (*ew)[in2m1 * in1 + ki] = (double) 1.0;
      for (kj = in2 - id3; kj < in2m1; kj++)
	if (ki < id4 || ki >= in1 - id2)
	  (*ew)[kj * in1 + ki] = (double) 0.5;
	else
	  (*ew)[kj * in1 + ki] = (double) 1.0;
    }

  for (kj = id1; kj < in2 - id3; kj++)
    for (ki = id4; ki < in1 - id2; ki++)
      {
	if ((double) (kj + 1) <= (double) (in2 + 1) / (double) 2.0)
	  {
	    if ((double) (ki + 1) <= (double) (in1 + 1) / (double) 2.0)
	      (*ew)[kj * in1 + ki] = (double) (ki + 1) / (ki + kj + 2);
	    else
	      (*ew)[kj * in1 + ki] = (double) (in1 - ki) / (in1 - ki + kj + 1);
	  }
	else
	  {
	    if ((double) (ki + 1) <= (double) (in1 + 1) / (double) 2.0)
	      (*ew)[kj * in1 + ki] = (double) (ki + 1) / (ki + 1 + in2 - kj);
	    else
	      (*ew)[kj * in1 + ki] = (double) (in1 - ki) / (in1 + in2 - ki - kj);
	  }
      }

  /* Ok. */

  goto out;


  /* Memory error. */

err101:
  *jstat = -101;
  s6err ("s1924", *jstat, kpos);
  goto out;

out:

  return;
}
