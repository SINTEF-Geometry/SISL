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
