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
