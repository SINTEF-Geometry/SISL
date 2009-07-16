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
 * $Id: s1753.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1753

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1753 (double et[], double ecf[], int in, int ik, int idim, double etr[],
       double ecfr[], int inr, double ecc[], double ecw[], int *jstat)

#else
void
s1753 (et, ecf, in, ik, idim, etr, ecfr, inr, ecc, ecw, jstat)
     double et[];
     double ecf[];
     int in;
     int ik;
     int idim;
     double etr[];
     double ecfr[];
     int inr;
     double ecc[];
     double ecw[];
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : 	To raise the description of a B-spline curve one order.
*
*
* INPUT      : 	et	- Description of knot vector of original description
*		ecf	- Coefficients of original description
*		in	- Number of vertices of original description
*		ik	- Order of original description
*		idim	- The dimension of the space in which the curve lies
*		etr	- Knot vector of the raised basis
*		inr	- Number of vertices in the raised curve
*		ecc	- Array for internal use only
*		ecw	-        ---- " ----
*
* OUTPUT     : 	ecfr	- Knots of the raised curve
*		jstat	- Status variable:
*				< 0  	: error
*				= 0	: OK.
*
* METHOD     : 	The order raising algorithm of Cohen, Lyche and Schumaker
*		is used.
*
*
* REFERENCES :	Fortran version:
*		T.Dokken, SI, 1984-06
*
*
* CALLS      :  s6err.
*
*
* WRITTEN BY : 	Christophe R. Birkeland, SI, 1991-07
* REWRITTEN BY :
* REVISED BY :
*
*********************************************************************
*/
{
  int ki, kj, kk, kl, kr, kstop;/* Loop control variables 		*/
  int kjmid, ikmid;		/* kjmid=(kj-1)*idim  ikmid=(ik-1)*idim */
  int kpos = 0;			/* Error position indicator		*/
  double ty1, ty2, tyi, tyik;	/* Parameters used in Main Loop		*/
  double dummy;
  double tden;

  *jstat = 0;


  /* Check input values. */

  if ((ik < 1) || (in <ik) ||(inr < (ik + 1)))
    goto err112;


  /* Initiate local variables. */

  kr = 1;
  for (kj = 1; kj <= inr; kj++)
    {

      /* Find kr, such that et[kr-1]<=etr[kj-1]<et[kr]	*/

      for (kr--; et[kr] <= etr[kj - 1]; kr++) ;


      /* Set ecc and ecw to zero. */

      for (ki = 0; ki < ik * idim; ki++)
	{
	  ecc[ki] = (double) 0.0;
	  ecw[ki] = (double) 0.0;
	}

      /* Initialize the remaining ecc and ecw entries. */

      kstop = MIN (ik, in +ik - kr);
      for (ki = MAX (0, ik - kr); ki < kstop; ki++)
	for (kl = 0; kl < idim; kl++)
	  {
	    dummy = ecf[(ki + kr - ik) * idim + kl];
	    ecc[ki * idim + kl] = dummy;
	    ecw[ki * idim + kl] = dummy;
	  }

      /* MAIN LOOP. */

      for (kk = ik - 1; kk > 0; kk--)
	{
	  ty1 = etr[kj + kk - 1];
	  ty2 = etr[kj + kk];
	  kstop = MAX (ik - kk, ik - kr);

	  for (ki = MIN (ik - 1, in +2 * ik - kk - kr - 1); ki >= kstop; ki--)
	    {
	      tyi = et[kr + ki - ik];
	      tyik = et[kr + ki + kk - ik];
	      tden = tyik - tyi;

	      for (kl = 0; kl < idim; kl++)
		{
		  ecc[ki * idim + kl] = ((ty2 - tyi) * ecc[ki * idim + kl] +
			   (tyik - ty2) * ecc[(ki - 1) * idim + kl]) / tden;
		  ecw[ki * idim + kl] = ((ty1 - tyi) * ecw[ki * idim + kl] +
			  (tyik - ty1) * ecw[(ki - 1) * idim + kl]) / tden +
		    ecc[ki * idim + kl];
		}
	    }
	}
      kjmid = (kj - 1) * idim;
      ikmid = (ik - 1) * idim;

      for (kl = 0; kl < idim; kl++)
	ecfr[kjmid + kl] = ecw[ikmid + kl] / ik;
    }

  goto out;


  /* Error in description of bases */

err112:
  *jstat = -112;
  s6err ("s1753", *jstat, kpos);
  goto out;

out:
  return;
}
