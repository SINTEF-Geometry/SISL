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
 * $Id: s1932.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1932

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1932 (int inbcrv, SISLCurve ** crvarr, double start, double stop, double *et,
       int in, int iordr, double **iright, int *jstat)
#else
void
s1932 (inbcrv, crvarr, start, stop, et, in, iordr, iright, jstat)
     int inbcrv;
     SISLCurve **crvarr;
     double start;
     double stop;
     double *et;
     int in;
     int iordr;
     double **iright;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To express a set of B-spline curves described in different
*		B-spline bases using a basis reflecting the continuity in
*		all the B-spline bases used in the description of the
*		different curves.
*
*
* INPUT      :	inbcrv	- Number of curves in the curve-set.
*		crvarr	- Array (length inbcrv) of pointers to the
*			  curves in the curve-set.
*		start	- Start parameter-value of B-spline basis to be made.
*		stop	- Stop parameter-value of B-spline basis to be made.
*		et	- The knot vector of the basis in which the curves
*			  are to be expressed.
*		in	- The number of B-spline bases in which the curves
*			  are to be expressed.
*		iordr	- The order of the basis in which the curves are to
*			  be expressed.
*
* OUTPUT     :  iright	- Array containing the right hand side of the equation
*			  system. Sequence first curve, second curve etc...
*               jstat   - Output status:
*                          < 0: Error.
*                          = 0: Ok.
*                          > 0: Warning.
*
* METHOD     : 	The description of the curves are lifted to the order given
*		as input. The lifted knot vectors are then mapped into the
*		parameter range given by astart and astop. Then the curves
*		are expressed in the refined basis given et. This basis must
*		be calculated to be refinement of the lifted and mapped knot
*		vectors of the input curves.
*
* REFERENCES :  Fortran revised version:
*               T.Dokken, SI, 1989-02
*
* CALLS      :  s1750,s1934,s1936,s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07
*
*********************************************************************
*/
{
  int ki, kj, kp;		/* Loop control parameters 		*/
  int kpos = 0;			/* Error position indicator		*/
  int kordr;			/* Highest curve-order used in
				 * description				*/
  int idim;			/* Dimension of the space where the
				 * curves lie				*/
  int kstat;			/* Status variable from lower level
				 * routines				*/
  SISLCurve *crv = SISL_NULL;	/* SISLCurve, sent to s1936		*/
  double *kdumcf = SISL_NULL;	/* Contains curve-coefficients sent over
				 * from s1936				*/


  *jstat = 0;

  /* Initailzation of variables */

  idim = crvarr[0]->idim;
  kordr = 0;

  for (ki = 0; ki < inbcrv; ki++)
    if ((crvarr[ki]->in <crvarr[ki]->ik) ||(crvarr[ki]->ik < 1))
      goto err112;


  /* Find highest order used in description */

  for (ki = 0; ki < inbcrv; ki++)
    kordr = MAX (kordr, crvarr[ki]->ik);

  if (kordr > iordr)
    goto err151;


  /* Lift the order of the description of the curve, then refine the
     description of the curve, and copy the result into the right hand
     side of the equation ssystem. */

  /* Allocate array kdumcf and output array iright */

  kdumcf = newarray (idim * in, DOUBLE);
  if (kdumcf == SISL_NULL)
    goto err101;
  *iright = newarray (in *idim * inbcrv, DOUBLE);
  if (*iright == SISL_NULL)
    goto err101;

  kp = 0;

  for (ki = 0; ki < inbcrv; ki++)
    {
      /* Increase the order of the basis */

      s1750 (crvarr[ki], iordr, &crv, &kstat);
      if (kstat < 0)
	goto error;


      /* Normalize the lifted basis */

      s1934 (crv->et, crv->in, crv->ik, start, stop, &kstat);
      if (kstat < 0)
	goto error;


      /* Express the curve using the refined basis */

      s1936 (crv, et, in, kdumcf, &kstat);
      if (kstat < 0)
	goto error;

      if (crv != SISL_NULL)
	freeCurve (crv);
      crv = SISL_NULL;	


      /* Copy coefficients into right hand side of equation system */

      for (kj = 0; kj < (in *idim); kj++)
	(*iright)[kj + kp] = kdumcf[kj];
      kp += in *idim;
    }
  goto out;

  /* Memory error */

  err101:
    *jstat = -101;
    s6err ("s1932", *jstat, kpos);
    goto out;

  /* Error in curve description */ 
  
  err112:
     *jstat = -112;
     s6err("s1932",*jstat,kpos);
     goto out;
 
  /* Error in input */

  err151:  
    *jstat = -151;
    s6err ("s1932", *jstat, kpos);
    goto out;

  /* Error in lower level routines */

  error:
    *jstat = kstat;
    s6err ("s1932", *jstat, kpos);
    goto out;

  out:
   /* Free scratch occupied by local array.  */
   
   if (kdumcf != SISL_NULL) freearray(kdumcf);

   return;
}
