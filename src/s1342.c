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
 * $Id: s1342.c,v 1.2 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1342

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1342(double ep[],double ev[],int im,int idim,int ipar,double epar[],
	   double eeps[],int ilend,int irend,double aepsco,int itmax,
	   SISLCurve **rc,double emxerr[],int *jstat)
#else
void s1342(ep,ev,im,idim,ipar,epar,eeps,ilend,irend,aepsco,itmax,
           rc,emxerr,jstat)
     double ep[];
     double ev[];
     int    im;
     int    idim;
     int    ipar;
     double epar[];
     double eeps[];
     int    ilend;
     int    irend;
     double aepsco;
     int    itmax;
     SISLCurve  **rc;
     double emxerr[];
     int    *jstat;
#endif
/*
***********************************************************
*
* Purpose: To compute the approximation to the data given by the points
*          ep and the derivatives (tangents) ev, and represent it as a
*          B-spline curve with parametrization determined by the parameter
*          ipar. The approximation is determined by first forming the cubic
*          hermite interpolant to the data, and then performing knot
*          removal on this initial approximation.
*
* Input :
*         Ep     - Array (length idim*im) comtaining the points
*                  to be approximated.
*         Ev     - Array (length idim*im) containing the derivatives
*                  of the points to be approximated.
*         Im     - The no. of data points.
*         Idim   - The dimension of the euclidean space in which the
*                  curve lies.
*         Ipar   - Flag indicating the type of paramterization to be used.
*                   = 1 : Parametrize by accumulated cord length (arc
*                         length parameterization for the piecewise
*                         linear interpolant.)
*                   = 2 : Uniform parametrization.
*                   = 3 : parameterization given by epar.
*                  If ipar<1 or ipar>3, it will be set to 1.
*         Epar   - Array (length im) containing a parameterization
*                  of the given data.
*         Eeps   - Array (length idim) giving the desired accuracy of
*                  the spline-approximation in each component.
*         Ilend  - The no. of derivatives that are not allowed to change
*                  at the left end of the curve.
*                  The (0 - (ilend-1)) derivatives will be kept fixed.
*                  If ilend <0, this routine will set it to 0.
*                  If ilend<ik, this routine will set it to ik.
*         Irend  - The no. of derivatives that are not allowed to change
*                  at the right end of the curve.
*                  The (0 - (irend-1)) derivatives will be kept fixed.
*                  If irend <0, this routine will set it to 0.
*                  If irend<ik, this routine will set it to ik.
*         Aepsco - Two numbers differing by a relative amount < aepsco,
*                  will in some cases be considered equal. A suitable
*                  value is just above the unit roundoff in the machine.
*                  On a VAX 1.0E-6 is a reasonable choice.
*                  The computations are not guaranteed to have relative
*                  accuracy less than aepsco.
*         Itmax  - Max. no. of iteration.
* Output:
*         Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning:
*         Rc     - Pointer to curve.
*         Emxerr - Array (length idim) (allocated outside this routine.)
*                  containing an upper bound for the pointwise error
*                  in each of the components of the spline-approximation.
*-
* Method: First the cubic hermite interpolant is computed by a call
*         to s1379 or s1380, and then knots are removed from this
*         interpolant.
*
* Calls: s1379, s1380, s1340, s6err.
*
* The main routine, s1340, is written by Knut M|rken,  Si.
* Written by: C.R. Birkeland, Si, April 1993.
**********************************************************
*/
{
  int stat = 0;                /* Error control parameters         */
  int kpos = 0;
  SISLCurve *ocurve = SISL_NULL;    /* Local spline curve               */

  /* Check Input */

  if (im < 2 || idim < 1) goto err103;
  if (ipar < 1 || ipar > 3) ipar = 1;

  /* Compute hermite interpolant */

  if (ipar == 3)
    s1379(ep, ev, epar, im, idim, &ocurve, &stat);
  else
    s1380(ep, ev, im, idim, ipar, &ocurve, &stat);
  if (stat < 0) goto error;

  /* Perform data reduction on the cubic hermite interpolant */

  s1340(ocurve, eeps, ilend, irend, aepsco, itmax, rc,
	emxerr, &stat);
  if (stat < 0) goto error;

  /* Success */

  *jstat = 0;
  goto out;
  
  /* Error in input */

 err103: 
  *jstat = -103;
  s6err("s1342",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */

 error: 
  *jstat = stat;
  s6err("s1342",*jstat,kpos);
  goto out;

  /* Exit */

 out:
  if (ocurve != SISL_NULL) freeCurve(ocurve);
  return;
}
