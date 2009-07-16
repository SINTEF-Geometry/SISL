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
