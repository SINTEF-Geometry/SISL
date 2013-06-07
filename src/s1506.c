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
 * $Id:
 *
 */
#define S1506

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1506(SISLSurf *ps1,int ider,int m1,double *x,int m2,double *y,
      double eder[],double norm[],int *jstat)
#else
void s1506(ps1,ider,m1,x,m2,y,eder,norm,jstat)
     SISLSurf *ps1;
     int      ider;
     int      m1;
     double   *x;
     int      m2;
     double   *y;
     double   eder[];
     double   norm[];
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the surface pointed at by ps1 over an m1 * m2 grid
*              of points (x[i],y[j]). Compute ider derivatives and normals
*              if suitable.
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*              m1     - Number of grid points in first direction.
*               x     - Array of x values of the grid.
*              m2     - Number of grid points in first direction.
*               y     - Array of y values of the grid.
*
* OUTPUT     : eder   - Array where the derivatives of the surface
*                       are placed, dimension
*                         idim * ((ider+1)(ider+2) / 2) * m1 * m2.
*                       The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2)
*                       derivative, etc. at point (x[0],y[0]),
*                       followed by the same information at (x[1],y[0]),
*                       etc.
*              norm   - Normals of surface. Is calculated if ider >= 1.
*                       Dimension is idim*m1*m2.
*                       The normals are not normalized.
*              jstat  - status messages
*                          = 2      : Surface is degenerate
*                                     at some point, normal
*                                     has zero length.
*                          = 1      : Surface is close to
*                                     degenerate at some point
*                                     Angle between tangents,
*                                     less than angular tolerance.
*                          = 0      : Ok.
*                          < 0      : Error.
*
* METHOD     : We call s1504 to pre-evaluate the B-splines then call
*              s1505 to multiply them with the coefficients.
*
*-
* CALLS      : s1504, s1505.
*
* WRITTEN BY : Michael Floater, SINTEF, May 1998.
*********************************************************************
*/
{
  int kstat=0;         /* Local status variable.                          */
  int kpos=0;          /* The position of error.                          */
  int n1,n2;           /* The number of B-splines accociated with the knot
		 	 vectors st1 and st2.                            */
  int k1,k2;           /* The polynomial order of the surface in the two
		 	 directions.                                     */
  int kdim;            /* The space dimension of the surface. */
  double *ebder1=SISL_NULL; /* Triple array of dimension (ider+1)*k1*m1
                         containing dericatives of B-splines. */
  double *ebder2=SISL_NULL; /* Triple array of dimension (ider+1)*k2*m2
                         containing dericatives of B-splines. */
  int *ileft1=SISL_NULL;    /* Array of dimension m1 containing the left knots
                         of the B-splines in x. */
  int *ileft2=SISL_NULL;    /* Array of dimension m2 containing the left knots
                         of the B-splines in y. */
  double *et1=SISL_NULL;    /* x knot vector. */
  double *et2=SISL_NULL;    /* y knot vector. */


  n1 = ps1 -> in1;
  n2 = ps1 -> in2;
  k1 = ps1 -> ik1;
  k2 = ps1 -> ik2;
  et1 = ps1 -> et1;
  et2 = ps1 -> et2;
  kdim = ps1 -> idim;

  /* Check the input. */
  if (kdim < 1) goto err102;
  if (k1 < 1 || k2 < 1) goto err115;
  if (n1 < k1 || n2 < k2) goto err116;
  if (ider < 0) goto err178;

  /* Pre-evaluate B-splines in x. */
  ebder1 = newarray((ider+1)*k1*m1,DOUBLE);
  if(ebder1 == SISL_NULL) goto err101;

  ileft1 = newarray(m1,INT);
  if(ileft1 == SISL_NULL) goto err101;

  s1504(et1,k1,n1,x,m1,ider,ebder1,ileft1,&kstat);
  if(kstat < 0) goto error;

  /* Pre-evaluate B-splines in y. */
  ebder2 = newarray((ider+1)*k2*m2,DOUBLE);
  if(ebder2 == SISL_NULL) goto err101;
  ileft2 = newarray(m2,INT);
  if(ileft2 == SISL_NULL) goto err101;

  s1504(et2,k2,n2,y,m2,ider,ebder2,ileft2,&kstat);
  if(kstat < 0) goto error;

  /* Multiply out with the coefficients. */

  s1505(ps1,ider,m1,m2,ebder1,ebder2,
      ileft1,ileft2,eder,norm,&kstat);
  if(kstat < 0) goto error;

  *jstat = 0;
  goto out;

  /* Not enough memory. */
 err101: *jstat = -101;
  s6err("s1506",*jstat,kpos);
  goto out;

  /* kdim less than 1. */
 err102: *jstat = -102;
  s6err("s1506",*jstat,kpos);
  goto out;

  /* Polynomial order less than 1. */
 err115: *jstat = -115;
  s6err("s1506",*jstat,kpos);
  goto out;

  /* Fewer B-splines than the order. */
 err116: *jstat = -116;
  s6err("s1506",*jstat,kpos);
  goto out;

  /* Illegal derivative requested. */
 err178: *jstat = -178;
  s6err("s1221",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("s1506",*jstat,kpos);
  goto out;

 out:
  /* Free memory. */
  if(ebder1 != SISL_NULL) freearray(ebder1);
  if(ileft1 != SISL_NULL) freearray(ileft1);
  if(ebder2 != SISL_NULL) freearray(ebder2);
  if(ileft2 != SISL_NULL) freearray(ileft2);

    return;
}
