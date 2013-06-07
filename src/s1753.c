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
