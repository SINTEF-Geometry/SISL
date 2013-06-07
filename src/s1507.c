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
 * $Id: s1507.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1507

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1507(SISLCurve **curves, int nc, int periodic,
      SISLCurve ***newcurves, int *jstat)
#else
void s1507(curves, nc, periodic, newcurves, jstat)
     SISLCurve **curves;
     int nc;
     int periodic;
     SISLCurve ***newcurves;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To reparametrize a sequence of B-spline curves so that
*              they are all parametrized over [0,1].
*              The newly parametrized curves are then "healed" (i.e. altered
*              near their end points) so that they connect with C1 continuity.
*
*              If periodic = 1, the healing will ensure that the
*              last and first curves also join with $C^1$ continuity.
*              Each B-spline curve must have at least four control points
*              and rational B-splines are not currently supported.
*
*              One application of this routine is in lofting a sequence
*              of piecewise G1 curves: a direct loft would only
*              guarantee a C0 lofted surface but by preprocessing (and
*              inevitably modifying) the curves to make them piecwise
*              C1 first will ensure that the lofted surface will be C1.
*
* INPUT      : curves   - An array of pointers to B-spline curves.
*              n        - Number of curves in the sequence.
*              periodic - Should the curves join periodically?
*
* OUTPUT     : newcurves - A pointer to a new array of curves.
*              jstat     - status messages
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : First the parameter intervals are scaled to [0,1] and
*              then the first two and last two control points in each
*              curve are altered to ensure both C1 continuity and that
*              the new curves are as close to the originals as possible
*              (in the sense of l2 in the control points).
*
* REFERENCES :
*
*
* CALLS      :
*
* WRITTEN BY : Michael Floater, SINTEF, September 1998.
*
**********************************************************************/
{
  int kpos = 0;			/* Position of error                        */

  int i,j;                      /* Loop variables.                          */

  int k,k1,k2;                  /* Order of B-spline curve                  */
  int n;                        /* Dimension of spline space                */
  int dim;                      /* Euclidean space dimension                */
  int inext;                    /* Usually i+1                              */
  int nheal;                    /* Number of joins between curves:
                                   = number of curves - 1 unless
                                   the periodic flag is 1.                  */

  SISLCurve **curves2;          /* Local array of SISL curves               */

  double startpar;              /* Start of parameter interval of a curve.  */
  double endpar;                /* Start of parameter interval of a curve.  */
  double len;                   /* Length of parameter interval of a curve. */
  double average;               /* average of coeff.s from adjacent curves. */
  double c1,c2,c3;              /* Three B-spline coefficients.             */
  double knotdiff1,knotdiff2;   /* numbers                                  */
  double lambda,denom;          /* more numbers                             */


  *jstat = 0;

  if( (curves2=newarray(nc,SISLCurve*)) == SISL_NULL) goto err101;
  memzero(curves2, nc, SISLCurve*);

  if (nc <= 0) goto err102;
  if (nc == 1) goto out; /* nothing to heal */

  dim = curves[0]->idim;

  /* Check input data. */
  for(i=0; i<nc; i++)
  {
    /* Check that the curves have the same dimension. */
    if(curves[i]->idim != dim) goto err102;

    /* Check that the curves are not rational. */
    if(curves[i]->ikind == 2 || curves[i]->ikind == 4) goto err102;

    /* Check that the curve has at least four control points
       (because we're going to tweak the first two to get C1 continuity
        with the previous curve and the last two to get C1 continuity
        with the next curve). */
    n = curves[i]->in;
    if(n < 4) goto err102;
  }

  /* Scale each parameter interval to [0,1]. */
  for(i=0; i<nc; i++)
  {
    /* Make a copy of the curve. */
    curves2[i] = copyCurve(curves[i]);

    /* Scale the knots to go from 0 to 1. */
    n = curves2[i]->in;
    k  = curves2[i]->ik;
    startpar = curves2[i]->et[k-1];
    endpar = curves2[i]->et[n];
    len = endpar - startpar;

    for(j=0; j<n+k; j++)
    {
      curves2[i]->et[j]  = (curves2[i]->et[j] - startpar) / len;
    }
  }

  /* How many healing operations are there? */

  nheal = (periodic == 0 ? nc - 1 : nc);

  /* Heal curves to get C^0 continuity. */

  for(i=0; i<nheal; i++)
  {
    inext = (i < nc-1 ? i+1 : 0);

    n = curves2[i]->in;
    for(k=0; k<dim; k++)
    {
      average = (curves2[i]->ecoef[dim*(n-1)+k] +
                 curves2[inext]->ecoef[k]) * 0.5;
      curves2[i]->ecoef[dim*(n-1)+k] = average;
      curves2[inext]->ecoef[k] = average;
    }
  }


  /* Heal curves to get C^1 continuity. */
  /* For each consecutive pair of curves, alter the three
     relevant control points c1, c2, c3
     (where c1 and c2 are the last two control points of the previous curve and
      c2 and c3 are the first two control points of the next curve:
      remember we have already fixed C0 continuity earlier in this routine)
     so that the three new points d1, d2, d3 are such that
     (c1-d1)^2 + (c2-d2)^2 + (c3-d3)^2 is minimized subject to
     the C1-continuity constraint which is of the form
       d2 = (1-lambda) d1 + lambda d3
     where lambda is a function of the two local knot intervals.
     The formula for d1,d2,d3 can be found by hand is used directly
     below. */

  for(i=0; i<nheal; i++)
  {
    inext = (i < nc-1 ? i+1 : 0);

    n = curves2[i]->in;
    k1 = curves2[i]->ik;
    k2 = curves2[inext]->ik;

    knotdiff1 = (curves2[i]->et[n+k1-2] - curves2[i]->et[n-1])
                    / (double)(k1-1);
    knotdiff2 = (curves2[inext]->et[k2] - curves2[inext]->et[1])
                    / (double)(k2-1);
    lambda = knotdiff2 / (knotdiff1 + knotdiff2);
    denom = 2.0 * (1.0 - lambda + lambda*lambda);

    for(k=0; k<dim; k++)
    {
      c1 = curves2[i]->ecoef[dim*(n-2)+k];
      c2 = curves2[inext]->ecoef[k];
            /* = curves2[i]->ecoef[dim*(n-1)+k] by C0 continuity. */
      c3 = curves2[inext]->ecoef[dim+k];

      curves2[i]->ecoef[dim*(n-2)+k]
        = ( (1.0 + lambda*lambda) * c1 + (1.0 - lambda) * c2
           - lambda * (1.0 - lambda) * c3 ) / denom;

      curves2[inext]->ecoef[dim+k]
        = ( - lambda * (1.0 - lambda) * c1 + lambda * c2
           + (1.0 + (1.0 - lambda)*(1.0 - lambda)) * c3 ) / denom;

      curves2[i]->ecoef[dim*(n-1)+k]
        = (1.0 - lambda) * curves2[i]->ecoef[dim*(n-2)+k]
          + lambda * curves2[inext]->ecoef[dim+k];

      curves2[inext]->ecoef[k]
        = (1.0 - lambda) * curves2[i]->ecoef[dim*(n-2)+k]
          + lambda * curves2[inext]->ecoef[dim+k];
    }
  }


  *newcurves = curves2;

  goto out;

  /* Error in space allocation */

err101:
  *jstat = -101;
  s6err ("s1507", *jstat, kpos);
  goto out;

  /* Error in input. */

err102:
  *jstat = -102;
  goto out;

out:

  return;
}
