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
 * $Id: s6rotmat.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define S6ROTMAT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6rotmat (double eorigo[], double exaxis[], double enorm[],
	  double ematrix[], int *jstat)
#else
void 
s6rotmat (eorigo, exaxis, enorm, ematrix, jstat)
     double eorigo[];
     double exaxis[];
     double enorm[];
     double ematrix[];
     int *jstat;
#endif
/*
*************************************************************************
*
* Purpose: To compute the transformation matrix (dim = 4X4), given
*	   points (of dim = 3), in the local coordinate system.
*
* Input:
*        Eorigo  - Origo of the local coordinate system.
*        Exaxis  - Point on (local) positive x-axis.
*        Enorm   - Normal to the xy-plane. (z-axis).
*
* Output:
*        Ematrix - Transformation metrix to the local coordinate system.
*        Jstat   - Output status.
*
* Calls: s6scpr, s6crss, s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Nov.91.
* After FORTRAN (P6INN).
*****************************************************************
*/
{
  double syaxis[3], sxaxis[3];
  double teps = (double) 0.000001;
  double tdum, tdum1, tdum2, tdum3;
  int kp;
  int idim = 3;

  int kpos = 0;

  *jstat = 0;


  /* First column will be the relative X-axis, second column
     the relative Y-axis, third column the relative Z-axis, and
     fourth column the local origo. The fourth line will be 0,0,0,1. */


  /* Find the local X-axis. */

  for (kp = 0; kp < idim; kp++)
    sxaxis[kp] = exaxis[kp] - eorigo[kp];


  /* Calculate the local Y-axis = (Z) X (X). */

  s6crss (enorm, sxaxis, syaxis);


  /* Normalise the three axis-vectors, and put into the matrix. */

  tdum = s6scpr (sxaxis, sxaxis, idim);
  tdum1 = sqrt (tdum);
  if (tdum1 < teps)
    goto err166;

  tdum = s6scpr (syaxis, syaxis, idim);
  tdum2 = sqrt (tdum);
  if (tdum2 < teps)
    goto err166;

  tdum = s6scpr (enorm, enorm, idim);
  tdum3 = sqrt (tdum);
  if (tdum3 < teps)
    goto err166;


  for (kp = 0; kp < idim; kp++)
    {
      ematrix[4 * kp] = sxaxis[kp] / tdum1;
      ematrix[4 * kp + 1] = syaxis[kp] / tdum2;
      ematrix[4 * kp + 2] = enorm[kp] / tdum1;
      ematrix[4 * kp + 3] = eorigo[kp];
      ematrix[12 + kp] = (double) 0.0;
    }

  ematrix[15] = (double) 1.0;

  goto out;

  /* Impossible to create matrix. */

err166:
  *jstat = -166;
  s6err ("s6rotmat", *jstat, kpos);
  goto out;

out:
  return;
}
