/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
