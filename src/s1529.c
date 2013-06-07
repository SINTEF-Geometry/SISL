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

#define S1529

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1529(double ep[],double eder10[],double eder01[],double eder11[],
	   int im1,int im2,int idim,int ipar,
	   SISLSurf **rsurf,int *jstat)
#else
void s1529(ep,eder10,eder01,eder11,
	   im1,im2,idim,ipar,rsurf,jstat)
     double ep[];
     double eder10[];
     double eder01[];
     double eder11[];
     int    im1;
     int    im2;
     int    idim;
     int    ipar;
     SISLSurf  **rsurf;
     int    *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute the cubic Hermite interpolant to the data given.
*          More specifically, given positions, 10, 01, and 11
*          derivatives at points of a rectangular grid, the routine
*          computes a cubic tensor-product B-spline interpolant to
*          the given data with double knots at each data (the first
*          knot vector will have double knots at all interior points
*          in epar1, quadruple knots at the first and last points,
*          and similarly for the second knot vector).
*
* INPUT:
*          ep     - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          eder10 - Array of dimension idim*im1*im2 containing the
*                   first derivative in the first parameter direction.
*
*          eder01 - Array of dimension idim*im1*im2 containing the
*                   first derivative in the second parameter direction.
*
*          eder11 - Array of dimension idim*im1*im2 containing
*                   the cross derivative (twist vector).
*
*          ipar   - Flag showing the desired parametrization to be used:
*                   = 1: Mean accumulated cord-length parameterization.
*                   = 2: Uniform parametrization.
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     The interpolation is accomplished by using a one dimensional
*     routine for cubic Hermite spline interpolation. First, the data
*     is considered to be im2 positional and derivative vectors on
*     two curves in idim * im1 dimensional space sampled at the
*     points of epar2.
*     The first of these has position vectors given by ep and
*     derivative vectors given by eder01, the second position vectors
*     given by eder10 and derivative vectors given by eder11.
*     Running these curves through the one dimensional cubic Hermite
*     spline interpolation routine then produces two cubic splines rpos
*     and rder with coefficients of dimension (idim * im1) * (2 * im2)
*     on the knot vector et2 which is just the points of epar2 with
*     multiplicity 2 for the interior points and 4 for the
*     end points.
*     These coefficients are then considered to be im1 position vectors
*     and derivative vectors on a curve in 2*idim*im2 dimensional
*     space (after an appropriate tranposition) sampled at the
*     points of epar1. Running this data through the one dimensional
*     cubic Hermite spline routine results in a cubic spline
*     with coefficients of dimension (2 * idim * im2) * (2 * im1)
*     with knot vector et1 similar to et2.
*     A transposition of these coefficients yields the B-spline
*     coefficients of the bicubic Hermite tensor-product spline
*     interpolant.

* REFERENCES :
*
* CALLS      : s1528, s1530
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, May 1993.
*
*********************************************************************
*/                                                               
{
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *par1=SISL_NULL;    /* Transposed positions (in rpos)              */
  double *par2=SISL_NULL;    /* Transposed derivatives (in rder)            */
  
  
  /* Check input */        
  
  if (ipar < 1 || ipar > 3) goto err102;
  
  /* Generate parametrizations */
  
  s1528(idim, im1, im2, ep, ipar, SISL_CRV_OPEN, SISL_CRV_OPEN,
	&par1, &par2, &kstat);
  if(kstat < 0) goto error;

  /* Hermite interpolation */

  s1530(ep,eder10,eder01,eder11,par1,par2,
	im1,im2,idim,rsurf,&kstat);
  if(kstat < 0) goto error;

  /* Success */
  
  *jstat = 0;
  goto out;  
  
  /* Error in input data. */

  err102: *jstat = -102;
    s6err("s1530",*jstat,kpos);
    goto out;
    
  /* Error in lower level routine. */

  error:  *jstat =kstat;
    s6err("s1530",*jstat,kpos);
    goto out;
  
  out:
    if(par1 != SISL_NULL) freearray(par1);
    if(par2 != SISL_NULL) freearray(par2);
    return;
}
