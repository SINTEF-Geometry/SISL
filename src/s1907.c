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
 * $Id: s1907.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1907

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1907(double *epoint, int *ntype, double *epar, int iopen, int icnsta,
	   int icnend, int inbpnt, int idim, double *opoint[], int *otype[],
	   double *opar[], int *knbpnt, int *jstat)
#else
void s1907(epoint, ntype, epar, iopen, icnsta, icnend, inbpnt, idim,
	   opoint, otype, opar, knbpnt, jstat)
     double *epoint;
     int *ntype;
     double *epar;
     int iopen;
     int icnsta;
     int icnend;
     int inbpnt;
     int idim;
     double *opoint[];
     int *otype[];
     double *opar[];
     int *knbpnt;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To transform from one type of interpolation conditions to
*	   another. Parameter values are also transformed.
*
* INPUT: epoint - Array (length idim*inbpnt) containing the 'old'
*                 interpolation conditions.
*        ntype  - Array (length inbpnt) containing type indicator for
*                 points/derivatives/second-derivatives:
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Derivative to next point.
*                  4 - Derivative to prior point.
*                  5 - Second derivative to next point.
*                  6 - Second derivative to prior point.
*                 13 - Start-point of tangent to next point.
*                 14 - End-point of tangent to prior  point.
*	 epar   - Array containing wanted paramerization for the
*		  interpolation points.
*        iopen  - Open / closed flag.
*        icnsta - Additional condition at the start of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at start.
*        icnend - Additional condition at the end of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at end.
*        inbpnt - No. of points/derivatives in the epoint array.
*        idim   - The dimension of the space in which the points lie.
*
* OUTPUT:
*	 opoint - The new interpolation conditions. (length idim*inbpnt)
*        otype  - Array containing kind of condition. Dimension
*                 is inbpnt.
*                 =  0 : A point is given.
*                 =  d : The d'th derivatative condition to the
*                        previous point is given.
*                 = -d : The d'th derivatative condition to the
*                        next point is given.
*	 opar   - The new parametrization. The derivative conditions are
*		  given the same parameter value as the point they belong to.
*        knbpar - no. of points.
*	 stat   - Status variable.
*                                     > 0     : warning
*                                     = 0     : ok
*                                     < 0     : error
:
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s6dist,s6err.
*
* WRITTEN BY :  Trond Vidar Stensby, SI, 1991-07
*
*********************************************************************
*/
{
  int kpos = 0;
  int count1, count2;		/* Loop control variables. */
  int dummy;
  int start, start2;
  int cpv;			/* Current parameter value. */

  *jstat = 0;


  /* Allocate output arrays. */

  if (icnsta && icnend)
    *knbpnt = inbpnt + 2;
  else if (icnsta || icnend)
    *knbpnt = inbpnt + 1;
  else
    *knbpnt = inbpnt;

  *opoint = newarray ((*knbpnt) * idim, DOUBLE);
  if (*opoint == SISL_NULL)
    goto err101;
  *otype = newarray (*knbpnt, INT);
  if (*otype == SISL_NULL)
    goto err101;

  if (iopen == SISL_CRV_OPEN)
    *opar = newarray (*knbpnt, DOUBLE);
  else
    *opar = newarray (*knbpnt + 1, DOUBLE);
  if (*opar == SISL_NULL)
    goto err101;


  /* Insert additional interpolation conditions. */

  if (icnsta != 0)
    {
      for (count1 = 0; count1 < idim; count1++)
	(*opoint)[count1] = (double) 0.0;

      (*otype)[0] = -2;
      (*opar)[0] = epar[0];
    }

  if (icnend != 0)
    {
      dummy = (*knbpnt) * idim;
      for (count1 = ((*knbpnt) - 1) * idim; count1 < dummy; count1++)
	(*opoint)[count1] = (double) 0.0;

      (*otype)[(*knbpnt) - 1] = 2;
    }

  /* Copy the rest of the points. */

  if (icnsta != 0)
    start = 1;
  else
    start = 0;

  cpv = -1;

  for (count1 = 0; count1 < inbpnt; count1++)
    {
      /* Transfor interpolation conditions. */

      if (ntype[count1] == 13)
	{
	  start2 = (count1 + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[(count1 + 1) * idim + count2] -
		epoint[count1 * idim + count2];
	    }
	}
      else if (ntype[count1] == 14)
	{
	  start2 = (count1 + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[count1 * idim + count2] -
		epoint[(count1 - 1) * idim + count2];
	    }
	}
      else
	{
	  start2 = (count1 + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[count1 * idim + count2];
	    }
	}

      /* Transform derivative indicators. */

      if (ntype[count1] == 1 || ntype[count1] == 2)
	{
	  (*otype)[count1 + start] = 0;
	  cpv++;
	  (*opar)[count1 + start] = epar[cpv];
	}
      else if (ntype[count1] == 3)
	{
	  (*otype)[count1 + start] = -1;
	  (*opar)[count1 + start] = epar[cpv + 1];
	}
      else if (ntype[count1] == 4)
	{
	  (*otype)[count1 + start] = 1;
	  (*opar)[count1 + start] = epar[cpv];
	}
      else if (ntype[count1] == 5)
	{
	  (*otype)[count1 + start] = -2;
	  (*opar)[count1 + start] = epar[cpv + 1];
	}
      else if (ntype[count1] == 6)
	{
	  (*otype)[count1 + start] = 2;
	  (*opar)[count1 + start] = epar[cpv];
	}
      else if (ntype[count1] == 13)
	{
	  (*otype)[count1 + start] = -1;
	  (*opar)[count1 + start] = epar[cpv + 1];
	}
      else if (ntype[count1] == 14)
	{
	  (*otype)[count1 + start] = 1;
	  (*opar)[count1 + start] = epar[cpv];
	}
    }

  if (icnend != 0)
    (*opar)[(*knbpnt) - 1] = epar[cpv];

  if (!(iopen == SISL_CRV_OPEN))
    {
     /* UJK, Following calculations not necessary, 
      * Parameter value given as input.
        
        Calculate distance between first and last point.
        
        for (count1 = 0; count1 < (*knbpnt) && (*otype)[count1] != 0; 
	     count1++);
        for (count2 = (*knbpnt) - 1; count2 >= 0 && (*otype)[count2] != 0; 
	     count2--);
        if (count1 > (*knbpnt) || count2 < 0) goto err164;
        
        (*opar)[*knbpnt] = (*opar)[(*knbpnt) - 1] +
        s6dist (&(*opoint)[count1 * idim], &(*opoint)[count2 * idim], idim); 

      */

     cpv++;
     (*opar)[*knbpnt] = epar[cpv];     
   }

  /* OK */

  goto out;

/*  UJK + CBI, not used:  
    No point conditions specified.

    err164:
      *jstat = -164;
      s6err ("s1907", *jstat, kpos);
      goto out;
*/

  /* Allocation error. */

err101:
  *jstat = -101;
  s6err ("s1907", *jstat, kpos);
  goto out;

out:
  return;
}
