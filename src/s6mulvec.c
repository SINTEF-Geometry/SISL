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
 * $Id: s6mulvec.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6MULVEC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6mulvec (double ematrix[], double evect[], double eright[])
#else
void
s6mulvec (ematrix, evect, eright)
     double ematrix[];
     double evect[];
     double eright[];
#endif
/*
*************************************************************************
*
* Purpose: To multiply a 4*4 matrix by a 3-D vector.
*
* Input:
*        Ematrix - The matrix (length 16).
*        Evect   - The vector (length 3).
*
* Output:
*        Eright   - The resulting 3-D vector.
*
* Written by: A.M. Ytrehus, SI Oslo Nov.91.
* After FORTRAN (P6mvec), written by: S. Meen  SI.
*****************************************************************
*/
{
  double svect[4];
  double sright[4];
  double tdum;
  int ki;
  int kj;

  /* Store the 3D evect-array in the 4D svect-array. */

  for (ki = 0; ki < 3; ki++)
    svect[ki] = evect[ki];

  svect[3] = (double) 1.0;


  /* Multiply the matrix by the vector svect.
     The result is stored in the 4-D local vector sright.*/

  for (ki = 0; ki < 4; ki++)
    {
      tdum = (double) 0.0;

      for (kj = 0; kj < 4; kj++)
	tdum += ematrix[4 * ki + kj] * svect[kj];

      sright[ki] = tdum;
    }


  /* 
   * Check if the bottom line of the matrix is 0,0,0,1.
   * In that case, just store the first three elements of svect in evect.
   * Else, divide all 3 first elements of svect by svect[3]. 
   * --------------------------------------------------------------------
   */

  for (ki = 0; ki < 3; ki++)
    eright[ki] = sright[ki];

  if (!(ematrix[12] == (double) 0.0) && (ematrix[13] == (double) 0.0) &&
      (ematrix[14] == (double) 0.0) && (ematrix[15] == (double) 1.0))
    {
      tdum = sright[3];

      for (ki = 0; ki < 3; ki++)
	eright[ki] = sright[ki] / tdum;
    }

  goto out;

out:
  return;
}
