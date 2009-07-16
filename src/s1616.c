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
 * $Id: s1616.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1616

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1616(double epoint[], int inbpnt, int idim, int eptyp[],
	   double econic[], int *jstat)
#else
void s1616(epoint, inbpnt, idim, eptyp, econic, jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    eptyp[];
     double econic[];
     int    *jstat;
#endif
/*
*************************************************************************
*
* PURPOSE: To find if the conic equation of the points in the array
*          epoint. The first and last points must be of type 1 or 2. The
*          intermediate points can be points or tangents, but no tangents
*          can be doubly defined. The routine s1614 should be run before this
*          routine to ensure that the description of the points/tangents in the
*          eptyp array are ok.
*
* INPUT:
*        Epoint - The points/tangents describing the conic.
*        Inbpnt - No. of points/tangents in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        Eptyp  - Type indicator for the points/tangents :
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Tangent to next point.
*                  4 - Tangent to prior point.
*
* Output:
*        Econic - The conic coefficients of the points in Epoint.
*        Jstat  - status messages:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method:
*	 If Inbpnt is 3, a circle is produced.
* 	 If Inbpnt is 4, a conic with the coeficient of x*y equal to
* 	 zero is produced.
* 	 If Inbpnt is 5, a general conic is produced.
* 	 The coeficients econic have the following meaning:
*	 econic(0)*x*x + 2*econic(1)*x*y + econic(2)*y*y +
* 	 econic(3)*x + 2*econic(4)*y + econic(5) = 0.0.
*-
* Calls: s6lufacp, s6lusolp, s1618, s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Oct.91.
* After FORTRAN, (P10634), written by: T. Dokken  SI.
* Revised by: J. Kaasa, SI, Aug. 1992 (Transposed smatrix and fixed the
*                         handling of singular matrix equations.)
* Revised by: J. Kaasa, SI, Aug. 1992 (Changed status 9005 to 105).
*****************************************************************
*/
{
  int ki, kj, kk;
  int kdim;
  int kn;
  double tmax = (double) HUGE;
  int kperm = 0;
  int ktyp;
  int kdir;
  double *smatrix = SISL_NULL;
  double *save = SISL_NULL;
  double *sright = SISL_NULL;
  double tx, ty, tdx, tdy;
  int *npiv = SISL_NULL;
  double solu[6];
  double tdiff;
  double tdum;
  int kpos = 0;
  int kstat = 0;

  *jstat = 0;


  /* Allocate local arrays. */

  smatrix = newarray (inbpnt * inbpnt, DOUBLE);
  if (smatrix == SISL_NULL)
    goto err101;

  sright = newarray (inbpnt, DOUBLE);
  if (sright == SISL_NULL)
    goto err101;

  save = newarray (inbpnt * inbpnt, DOUBLE);
  if (save == SISL_NULL)
    goto err101;

  npiv = newarray (inbpnt, INT);
  if (npiv == SISL_NULL)
    goto error;

  kn = inbpnt + 1;
  kdim = inbpnt * inbpnt;


  /* Build interpolation matrix to find the coeficients of the conic. */

  if (inbpnt == 3)
    {

      /* A circle is to be produced. We put econic(0) = econic(2) = 1.0,
	 and econic(1) = 0.0. (This is not done now but later.) */

      for (ki = 0; ki < inbpnt; ki++)
	{
	  kk = idim * ki;
	  ktyp = eptyp[ki];

	  if (ktyp < 3)
	    {
	      tx = epoint[kk];
	      ty = epoint[kk + 1];

	      smatrix[inbpnt*ki] = ((double) 2.0) * tx;
	      smatrix[inbpnt*ki + 1] = ((double) 2.0) * ty;
	      smatrix[inbpnt*ki + 2] = (double) 1.0;
	      sright[ki] = -(tx * tx + ty * ty);
	    }
	  else if (ktyp > 2)
	    {

	      /* Derivative condition. */

	      kdir = 1;
	      if (ktyp == 4)
		kdir = -1;

	      tx = epoint[kk + idim * kdir];
	      ty = epoint[kk + idim * kdir + 1];
	      tdx = epoint[kk];
	      tdy = epoint[kk + 1];

	      smatrix[inbpnt*ki] = ((double) 2.0) * tdx;
	      smatrix[inbpnt*ki + 1] = ((double) 2.0) * tdy;
	      smatrix[inbpnt*ki + 2] = (double) 0.0;
	      sright[ki] = -tdx * tx * ((double) 2.0) - tdy * ty * ((double) 2.0);
	    }
	}
    }
  else if (inbpnt == 4)
    {
      /* A conic with xy term equal to zero is to be produced.
         We put econic(1) = 0.0. (This is not done now but later.) */

      for (ki = 0; ki < inbpnt; ki++)
	{
	  kk = idim * ki;
	  ktyp = eptyp[ki];
	  if (ktyp < 3)
	    {
	      tx = epoint[kk];
	      ty = epoint[kk + 1];

	      smatrix[inbpnt*ki] = tx * tx;
	      smatrix[inbpnt*ki + 1] = ty * ty;
	      smatrix[inbpnt*ki + 2] = ((double) 2.0) * tx;
	      smatrix[inbpnt*ki + 3] = ((double) 2.0) * ty;
	      sright[ki] = -(double) 1.0;
	    }
	  else if (ktyp > 2)
	    {

	      /* Derivative condition. */

	      kdir = 1;
	      if (ktyp == 4)
		kdir = -1;

	      tx = epoint[kk + idim * kdir];
	      ty = epoint[kk + idim * kdir + 1];
	      tdx = epoint[kk];
	      tdy = epoint[kk + 1];

	      smatrix[inbpnt*ki] = ((double) 2.0) * tdx * tx;
	      smatrix[inbpnt*ki + 1] = ((double) 2.0) * tdy * ty;
	      smatrix[inbpnt*ki + 2] = ((double) 2.0) * tdx;
	      smatrix[inbpnt*ki + 3] = ((double) 2.0) * tdy;
	      sright[ki] = (double) 0.0;
	    }
	}
    }
  else if (inbpnt == 5)
    {

      /* A general conic is to be produced. */

      for (ki = 0; ki < inbpnt; ki++)
	{
	  kk = idim * ki;
	  ktyp = eptyp[ki];
	  if (ktyp < 3)
	    {
	      tx = epoint[kk];
	      ty = epoint[kk + 1];

	      smatrix[inbpnt*ki] = tx * tx;
	      smatrix[inbpnt*ki + 1] = ((double) 2.0) * tx * ty;
	      smatrix[inbpnt*ki + 2] = ty * ty;
	      smatrix[inbpnt*ki + 3] = ((double) 2.0) * tx;
	      smatrix[inbpnt*ki + 4] = ((double) 2.0) * ty;
	      sright[ki] = -(double) 1.0;
	    }
	  else if (ktyp > 2)
	    {

	      /* Derivative condition. */

	      kdir = 1;
	      if (ktyp == 4)
		kdir = -1;

	      tx = epoint[kk + idim * kdir];
	      ty = epoint[kk + idim * kdir + 1];
	      tdx = epoint[kk];
	      tdy = epoint[kk + 1];

	      smatrix[inbpnt*ki] = ((double) 2.0) * tdx * tx;
	      smatrix[inbpnt*ki + 1] = ((double) 2.0) * tdy * tx + ((double) 2.0) * tdx * ty;
	      smatrix[inbpnt*ki + 2] = ((double) 2.0) * tdy * ty;
	      smatrix[inbpnt*ki + 3] = ((double) 2.0) * tdx;
	      smatrix[inbpnt*ki + 4] = ((double) 2.0) * tdy;
	      sright[ki] = (double) 0.0;
	    }
	}
    }

  /* Solve the equation system. If the system is not solvable, interchange
     the r.h.side with one of the colomns on the l.h. side. */

  for (ki = 0; ki < kn; ki++)
    {
      /* Remember interpolation matrix. */

      for (kj = 0; kj < kdim; kj++)
	save[kj] = smatrix[kj];


      /* Find solution.
	 s6lusolp put the soltion into the r.h. side array given as input,
	 thus store the right hand side in econic. */

      for (kj = 0; kj < inbpnt; kj++)
	econic[kj] = sright[kj];

      s6lufacp (smatrix, npiv, inbpnt, &kstat);

      if (kstat >= 0 && kstat != 1)
        s6lusolp (smatrix, econic, npiv, inbpnt, &kstat);
      kstat = 0;

      /* If we are here, we have been able to find a solution, econic.
	 Test this. First, restore smatrix and sright. */

      for (kj = 0; kj < kdim; kj++)
	smatrix[kj] = save[kj];

      s1618 (smatrix, sright, econic, inbpnt, &tdiff);

      if (tdiff < tmax)
	{
	  /* If we are here, we have found the best solution until now. */

	  tmax = tdiff;
	  kperm = ki;

	  for (kj = 0; kj < inbpnt; kj++)
	    solu[kj] = econic[kj];
	  
          if (inbpnt == 3) break;
	}
	
      if (ki < (kn - 1))
	{
	  /* Change the right hand side with one column of the matrix. */

	  for (kj = 0; kj < inbpnt; kj++)
	    {
	      kk = inbpnt * kj;
	      tdum = -sright[kj];
	      sright[kj] = -smatrix[kk + ki];
	      smatrix[kk + ki] = tdum;
	    }
	}
    }

  /* If no solution found, make straight line. */

  if (tmax > 0.0001)
    {
      *jstat = 105;
      econic[0] = 0.;
      econic[1] = 0.;
      econic[2] = 0.;
      econic[3] = (epoint[1] - epoint[2*inbpnt - 1])/2.;
      econic[4] = (epoint[2*inbpnt - 2] - epoint[0])/2.;
      econic[5] = epoint[0]*epoint[2*inbpnt - 1] - 
	 epoint[1]*epoint[2*inbpnt - 2];
      goto out;
    }

  /* Get the remembered best solution. */

  for (kj = 0; kj < inbpnt; kj++)
    econic[kj] = solu[kj];

  econic[inbpnt] = (double) 1.0;


  /* If the order of the colomns were changed, substitute back again. */

  if (kperm != 0)
    {
      for (ki = 1; ki <= kperm; ki++)
	{
	  kj = kperm - ki;
	  tdum = econic[kj];
	  econic[kj] = econic[inbpnt];
	  econic[inbpnt] = tdum;
	}
    }

  /* Expand to full conic description when Inbpnt is less than 5. */

  if (inbpnt == 3)
    {
      /* Circle. Only x, y and constant coeff. calculated. */

      for (ki = 0; ki < inbpnt; ki++)
	econic[ki + inbpnt] = econic[ki];

      econic[0] = (double) 1.0;
      econic[1] = (double) 0.0;
      econic[2] = (double) 1.0;
    }
  else if (inbpnt == 4)
    {
      /* Conic with xy constant zero. */

      for (ki = 5; ki > 1; ki--)
	econic[ki] = econic[ki - 1];

      econic[1] = (double) 0.0;
    }

  goto out;

  /* Error in allocation. */

err101:
  *jstat = -101;
  s6err ("s1616", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1616", *jstat, kpos);
  goto out;

out:

  /* Free arrays. */

  if (smatrix != SISL_NULL)
    freearray (smatrix);
  if (save != SISL_NULL)
    freearray (save);
  if (sright != SISL_NULL)
    freearray (sright);
  if (npiv != SISL_NULL)
    freearray (npiv);

  return;
}
