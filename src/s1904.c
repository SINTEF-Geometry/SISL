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
 * $Id: s1904.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1904

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1904 (double epar[], int in, int ik, int cuopen, double *eknots[], int *jstat)

#else
void
s1904 (epar, in, ik, cuopen, eknots, jstat)
     double epar[];
     int in;
     int ik;
     int cuopen;
     double *eknots[];
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To produce the knot vector of a B-spline basis satisfying
*              the interpolation requirements reflected in the epar
*              array. The knots is placed in such a way that the
*              interpolation conditions will lie close to the points
*              corresponding to the top point of the basis functions,
*              using the Laksaa method. This knot placement in particularily
*              suited when the coefficients are given, and a knot vector
*              is to be found.
*
* INPUT      : epar   - Array containing a parametrization of the
*                       interpolation conditions. Each interpolation
*                       conditions has got a distinct parameter value
*                       exept from the cases where several conditions
*                       are conflicing. In that case a multiple parameter
*                       value indicates the need of a multiple knot. The
*                       parameter values are sorted in increasing order.
*                       The dimension of the array is in.
*              in     - Number of interpolation conditions.
*              ik     - Order of B-spline basis.
*	       cuopen - Open/closed curve.
*
* OUTPUT     : eknots - The produced knot vector. The dimension of
*                       the array is in+ik.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY :  Vibeke Skytt, SI, 91-04.
* REVISED BY :	Trond Vidar Stensby, SI, 91-07
*
*********************************************************************
*/
{
  int ki;			/* Counter used to traverse the knot vector.  	*/
  int kpar1, kpar2;		/* Counters used to traverse the parametrization
				   array 					*/
  int k1, k2;			/* Numers with wich to multiply interval size.	*/
  int kstop;			/* Control variable of loop.                  	*/
  int kmult;			/* Multiplisity of knot.                      	*/
  double tdegree;		/* Degree of wanted basis.  			*/
  double tprev;			/* Value of previous knot.           	        */
  double curr;			/* Value of current knot.           	        */
  double tval1;			/* Start parameter value.          		*/
  double tval2;			/* End parameter value.            		*/
  double tparint;		/* The parameter interval.			*/
  int kpos = 0;

  *jstat = 0;


  if (cuopen)
    {
      /* O P E N   C U R V E */

      *eknots = newarray (in +ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      tdegree = (double) (ik - 1);
      tval1 = epar[0];
      tval2 = epar[in -1];

      /* Store a knot of multiplisity equal to the order in the start
	 of the curve. The value of the knot is equal to the value
	 of the start parameter value.  */

      for (ki = 0; ki < ik; ki++)
	(*eknots)[ki] = tval1;

      if (in -3 * ik + 4 < 0)
	{
	  /* All internal knots are affected by both endpoints.
	     Place internal knots.  */

	  for (; ki < in; ki++)
	    (*eknots)[ki] = (*eknots)[ki - 1]
	      + ((double) (2 * ik - ki - 1) * (epar[ki - ik + 1] - epar[ki - ik])
		 + (double) (ki - in +ik - 1) *(epar[ki - 1] - epar[ki - 2])
		 + epar[ki - 2] - epar[ki - ik + 1]) / tdegree;
	}
      else
	{

	  /* Place the knots close to the startpoint.  */

	  for (; ki < 2 * ik - 2; ki++)
	    (*eknots)[ki] = (*eknots)[ki - 1]
	      + ((double) (2 * ik - ki - 1) * (epar[ki - ik + 1] - epar[ki - ik])
		 + epar[ki - 1] - epar[ki - ik + 1]) / tdegree;

	  /* Place the inner knots.  */

	  for (; ki < in -ik + 2; ki++)
	    (*eknots)[ki] = (*eknots)[ki - 1]
	      + (epar[ki - 1] - epar[ki - ik]) / tdegree;

	  /* Place the knots close to the endpoint. */

	  for (; ki < in; ki++)
	    (*eknots)[ki] = (*eknots)[ki - 1]
	      + ((double) (ki - in +ik - 1) *(epar[ki - 1] - epar[ki - 2])
		 + epar[ki - 2] - epar[ki - ik]) / tdegree;
	}

      /* Store a knot of multiplisity equal to the order at the end of
	 the curve. The value of the knot is equal to the value of the
	 end parameter value.  */

      for (ki = 0; ki < ik; ki++)
	(*eknots)[in +ki] = tval2;
    }
  else
    {
      /* C L O S E D   C U R V E */

      *eknots = newarray (in +2 * ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      kstop = in +2 * ik - 1;
      tparint = epar[in] -epar[0];
      tdegree = (double) (ik - 1);

      /* Set the knot starting the interval of full basis.  */

      (*eknots)[ik - 1] = (ik % 2 == 0) ? epar[0] : (double) 0.5 *(epar[0] + epar[1]);

      /* Compute the following knots. */

      for (ki = ik, kpar1 = ki - 1, kpar2 = ki - ik, k1 = k2 = 0;
	   ki < kstop; ki++, kpar1++, kpar2++)
	{
	  /* Check that the indexes into the parameter array is legal. */

	  while (kpar1 < 0)
	    kpar1 += in, k1--;
	  while (kpar1 > in)
	    kpar1 -= in, k1++;
	  while (kpar2 < 0)
	    kpar2 += in, k2++;
	  while (kpar2 > in)
	    kpar2 -= in, k2--;

	  /* Compute knot.  */

	  (*eknots)[ki] = (*eknots)[ki - 1]
	    + (epar[kpar1] - epar[kpar2] + (double) (k1 + k2) * tparint) / tdegree;
	}

      /* Compute the first ik-1 knots. */

      for (ki = ik - 1, kpar1 = ki - 1, kpar2 = ki - ik, k1 = k2 = 0;
	   ki > 0; ki--, kpar1--, kpar2--)
	{
	  /* Check that the indexes into the parameter array is legal. */

	  while (kpar1 < 0)
	    kpar1 += in, k1--;
	  while (kpar1 > in)
	    kpar1 -= in, k1++;
	  while (kpar2 < 0)
	    kpar2 += in, k2++;
	  while (kpar2 > in)
	    kpar2 -= in, k2--;

	  /* Compute knot.  */

	  (*eknots)[ki - 1] = (*eknots)[ki]
	    - (epar[kpar1] - epar[kpar2] + (double) (k1 + k2) * tparint) / tdegree;
	}
    }


  /* Check that the produced knots are in increasing order and that
     the multiplicity is not greater than ik.                       */

  if (cuopen)
    kstop = in +ik;

  for (ki = 1, tprev = (*eknots)[0], kmult = 0; ki < kstop; ki++, tprev = curr)
    {
       curr = (*eknots)[ki];
      kmult++;
      if (tprev > curr)
	goto err112;		/* Decreasing parameter value. */
      if (tprev < curr)
	kmult = 1;
      if (kmult > ik)
	goto err112;		/* Knot multiplisity greater than order. */
    }

  /* The knot vector is produced.  */

  goto out;


  /* Error in scratch allocation. */

err101:
  *jstat = -101;
  s6err ("s1904", *jstat, kpos);
  goto out;

  /* Error in the knot vector.  */

err112:
  *jstat = -112;
  s6err ("s1904", *jstat, kpos);
  goto out;

out:
  return;
}
