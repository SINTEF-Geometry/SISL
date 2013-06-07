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
 * $Id: s1935.c,v 1.4 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1935

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void s1935 (double et1[], int in1, double et2[], int in2,
	    double *knt[], int *in, int ik, int *jstat)
#else
void
s1935 (et1, in1, et2, in2, knt, in, ik, jstat)
     double *et1;
     int in1;
     double *et2;
     int in2;
     double *knt[];
     int *in;
     int ik;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To produce a knot vector that is the union of two knot
*	   vectors. The original knot vectors and the produced knot
*	   vector are of the same order.
*
*
* INPUT:   et1	- The first knot vector.
*	   in1	- The number of degrees of freedom in the
*	          B-basis given by the first knot vector.
*	   et2	- The second knot vector.
*	   in2	- The number of degrees of freedom in the
*		  B-basis given by the second knot vector.
*	   ik	- The order of the knot vector to be produced.
*
* OUTPUT:  jstat - Output status:
*                   < 0: Error.
*		    = 0: Ok.
*                   > 0: Warning.
*
* METHOD:
*
* REFERENCES:  Fortran version:
*              T.Dokken, SI, 1981-11
*
* CALLS: s6err.
*
* WRITTEN BY:  Christophe R. Birkeland, SI, 1991-07
* CORRECTED BY: Vibeke Skytt, SI, 92-10. Removed exact test on equality
*                                        of knots.
* CORRECTED BY: Christophe R. Birkeland, SI, 93-05.
*         Reinstalled exact test on equality of knots. Necessary
*         for correct working of routine s1931-s1937.
* CORRECTED BY: Paal Fugelli, SINTEF, 1994-07.
*         Changed the setting of 'curr', for the while loop, to avoid
*         over running array bounds (address error) and added DEQUAL in
*         tests for knot equality.
*
*********************************************************************
*/
{
  int kpos = 0;			/* Error position indicator.		*/
  int pek1, pek2;		/* Index for array et1 and et2. 	*/
  int stop1;			/* Length of et1, equals in1+ik.	*/
  int stop2;			/* Length of et2, equals in2+ik.	*/
  double curr;			/* Parameter used in calculation of
				   new knots.				*/

  *jstat = 0;

  /* Test if legal input */

  if (ik < 1) goto err110;
  if ((in1 < ik) || (in2 < ik)) goto err111;


  /* Allocate array for new knot vector */

  if((*knt = newarray (2 * ik + in1 + in2, DOUBLE))==SISL_NULL) goto err101;

  /* Test if input knot vectors degenerate */

  if (et1[ik - 1] >= et1[in1]) goto err112;

  if (et2[ik - 1] >= et2[in2]) goto err112;

  /* PRODUCTION OF KNOTS */

  *in = 0;
  /* PFU 05/07-94  curr = MIN (et1[0], et2[0]); */
  pek1 = 0;
  pek2 = 0;
  stop1 = in1 + ik;
  stop2 = in2 + ik;

  while ((pek1 < stop1) && (pek2 < stop2))
    {
      /* Test if error in knot vector */

      curr = MIN (et1[pek1], et2[pek2]);

      if ((et1[pek1] < curr) || (et2[pek2] < curr)) goto err112;

      if ( DEQUAL(et1[pek1],curr) ) pek1++;
      if ( DEQUAL(et2[pek2],curr) ) pek2++;
      (*knt)[*in] = curr;
      (*in) ++;
    }

  /* Some knots may remain in one of the arrays */

  if ((pek1 < stop1) || (pek2 < stop2))
    {
      if (pek1 >= stop1)
	{
	  for (; pek2 < stop2; pek2++, (*in) ++)
	    (*knt)[*in] = et2[pek2];
	}
      else
	for (; pek1 < stop1; pek1++, (*in) ++)
	  (*knt)[*in] = et1[pek1];
    }

  /* Knots produced */

  *in -=ik;
  *knt = increasearray (*knt, *in +ik, DOUBLE);
  if (*knt == SISL_NULL) goto err101;
  goto out;


  /* Memory error */

  err101:
    *jstat = -101;
    s6err ("s1935", *jstat, kpos);
    goto out;

  /* Error in description of B-spline curve: Order less than 1.*/

  err110:
    *jstat = -110;
    s6err ("s1935", *jstat, kpos);
    goto out;

  /* No. of vertices less than order. */

  err111:
    *jstat = -111;
    s6err ("s1935", *jstat, kpos);
    goto out;

  /* Error in knot vector */

  err112:
    *jstat = -112;
    s6err ("s1935", *jstat, kpos);
    goto out;

  out:
    return;
}
