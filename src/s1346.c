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
 * $Id: s1346.c,v 1.4 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1346

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1346(double ep[],int im1,int im2,int idim,int ipar,double epar1[],
	   double epar2[],double eeps[],int nend[],double edgeps[],
	   double afctol,double aepsco,int iopt,int itmax,int ik1,int ik2,
	   SISLSurf **rs,double emxerr[],int *jstat)
#else
void s1346(ep,im1,im2,idim,ipar,epar1,epar2,eeps,nend,
           edgeps,afctol,aepsco,iopt,itmax,ik1,ik2,
           rs,emxerr,jstat)
     double ep[];
     int    im1;
     int    im2;
     int    idim;
     int    ipar;
     double epar1[];
     double epar2[];
     double eeps[];
     int    nend[];
     double edgeps[];
     double afctol;
     double aepsco;
     int    iopt;
     int    itmax;
     int    ik1;
     int    ik2;
     SISLSurf   **rs;
     double emxerr[];
     int    *jstat;
#endif
/*
********************************************************************
*
* Purpose: To compute a tensor-product spline-approximation of order
*          (ik1,ik2) to the rectangular array of idim-dimensional
*          points given by ep.
*
* Input : Ep     - Array (length idim*im1*im2) containing the points
*                  to be approximated.
*         Im1    - The no. of points in first parameter-direction.
*         Im2    - The no. of points in second parameter-direction.
*         Idim   - The no. of components of each input-point.
*                  The approximation will be a parametric surface
*                  situated in the idim-dimensional euclidean space
*                  (usually 3).
*         Ipar   - Flag determining the parametrization of the data points:
*                   = 1: Mean accumulated cord-length parameterization.
*                   = 2: Uniform parametrization.
*                   = 3: Parametrization given by epar1 and epar2.
*         Epar1  - Array (length im1) containing a parametrization
*                  in the first parameter-direction. (Will only
*                  be used if ipar=3).
*         Epar2  - Array (length im2) containing a parametrization
*                  in the second parameter-direction. (Will only
*                  be used if ipar=3).
*         Eeps   - Array (length idim) containing the max. permissible
*                  deviation of the approximation from the given data
*                  points, in each of the components. More specifically,
*                  the approximation will not deviate more than eeps(kdim)
*                  in component no. kdim, from the bilinear approximation
*                  to the data.
*         Nend   - Array (length 4) giving the no. of derivatives to be
*                  kept fixed along each edge of the bilinear interpolant.
*                  The numbering of the edges is the same as for edgeps below.
*                  All the derivatives of order < (nend(i)-1) will be kept
*                  fixed along the edge i. Hence nend(i)=0 indicates that
*                  nothing is to be kpet fixed along edge i. (Used by the
*                  data reduction routine.)
*                  To be kept fixed here means to have error less than edgeps.
*                  In general, it is impossible to remove any knots and keep
*                  an edge completely fixed.
*         Edgeps - Array (length idim*4) containing the max. deviation from
*                  the bilinear interpolant which is acceptable along the
*                  edges of the surface.
*                  Edgeps(1,i):edgeps(idim,i) gives the tolerance along
*                  the edge corresponding to the i-th parameter having
*                  one of it`s extremal-values.
*                   i=1: min value of first parameter.
*                   i=2: max value of first parameter.
*                   i=3: min value of second parameter.
*                   i=4: max value of second parameter.
*                  (Used by the data-reduction routine.)
*                  Edgeps(kp,i) will only have significance if nend(i)>0.
*         Afctol - 0.0 >= afctol <= 1.0.
*                  Afctol indicates how the tolerance is to be shared
*                  between the two data-reduction stages. For the linear
*                  reduction, a tolerance of afctol*eeps will be used,
*                  while a tolerance of (1.0-afctol)*eeps will be used
*                  during the final data reduction (similarly for edgeps.)
*                  Default is 0.
*         Aepsco - Two numbers differing by a relative amount less than
*                  aepsco, will in some cases be considered equal.
*                  A suitable value is just above the unit roundoff
*                  of the computer.
*                  The computations are not guaranteed to have relative
*                  accuracy less than aepsco.
*          Iopt  - Flag indicating the order in which the data-reduction
*                  is to be performed:
*                   = 1: Remove knots inparameter-direction 1 only.
*                   = 2: Remove knots inparameter-direction 2 only.
*                   = 3: Remove knots first in parameter-direction 1 and
*                        then in parameter-direction 2.
*                   = 4: Remove knots first in parameter-direction 2 and
*                        then in parameter-direction 1.
*         Itmax  - Max. no. of iterations in the data-reduction.
*         Ik1    - The order of the approximation in first
*                  parameter-directon.
*         Ik2    - The order of the approximation in second
*                  parameter-directon.
*
* Output:
*         Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning:
*         Rs     - Pointer to surface.
*         Emxerr - Array (length idim) (allocated outside this routine.)
*                  containing the error in the approximation to the data.
*                  This is guaranteed upper bound on the max. deviation
*                  in each component, between the final approximation
*                  and the bilinear spline-approximation to the original data.
*
* Method:
*        First the bilinear interpolant to the data is computed, using the
*        parameterization given by ipar, and knots are removed from this
*        initial approximation by a call to the data-reduction routine for
*        surfaces. Then the order is raised to (ik1,ik2) and the final data
*        reduction is performed.
*-
* Calls: s1345, s1350, s6chpar, s6err.
*
* Written by: C.R.Birkeland, Si, April 1993.
* The main routine, s1345, is written by: Knut M|rken,  SI.
* Changed by: Per OEyvind, SINTEF, 1994-11.
*             Removed following memory leaks:
*              1) Improper use of copy flag to newSurf()
*              2) Forgetting to free temp array after using icopy == 1
**********************************************************************
*/
{
  int in1,in2;                /* Number of vertices                   */
  int newin1, newin2;
  int fouridim=4*idim;
  int i;                      /* Loop control parameters              */
  int stat=0, kpos=0;         /* Error message parameters             */
  double *par1 = SISL_NULL;
  double *par2 = SISL_NULL;
  double *knot1 = SISL_NULL;       /* Knot vectors in 1 and 2. par.dir.    */
  double *knot2 = SISL_NULL;
  double *error1 = SISL_NULL;      /* Arrays for error storage             */
  double *error2 = SISL_NULL;
  double *maxerr = SISL_NULL;
  double *newcoeff = SISL_NULL;    /* Coefficients array                   */
  SISLCurve *ocurve1 = SISL_NULL;  /* Used to store local curves           */
  SISLCurve *ocurve2 = SISL_NULL;
  SISLSurf *osurf1 = SISL_NULL;    /* Used to store local surfaces         */
  SISLSurf *osurf2 = SISL_NULL;

  /* Check Input */

  if (im1 < 2 || im2 < 2 || ik1 < 1 || ik2 < 1 || idim < 1)
    goto err103;
  if (ipar < 1 || ipar > 3) ipar = 1;

  if (ipar != 3)
    {
      /* Generate parametrization */

      s1528(idim, im1, im2, ep, ipar, SISL_CRV_OPEN, SISL_CRV_OPEN,
	    &par1, &par2, &stat);
      if (stat<0) goto error;
    }
  else
    {
      /* Parametrization is passed as parameter */

      par1 = epar1;
      par2 = epar2;
    }

  /* Represent input (points) as a surface of
   * order 2 (linear) in both directions.
   * First, generate knot vectors */

  knot1 = newarray(im1+2, DOUBLE);
  knot2 = newarray(im2+2, DOUBLE);
  if(knot1 == SISL_NULL || knot2 == SISL_NULL) goto err101;
  memcopy(&knot1[1],par1,im1,DOUBLE);
  memcopy(&knot2[1],par2,im2,DOUBLE);
  knot1[0] = knot1[1];
  knot2[0] = knot2[1];
  knot1[im1+1] = knot1[im1];
  knot2[im2+1] = knot2[im2];
  osurf1 = newSurf(im1, im2, 2, 2, knot1, knot2, ep,
		   1,idim, 1);
  if (osurf1 == SISL_NULL) goto err101;
  if (knot1) freearray(knot1);
  if (knot2) freearray(knot2);

  /* Compute tolerance vectors for linear reduction
   * Both max deviation of surface and max dev. of edges */

  maxerr = newarray(idim, DOUBLE);
  error1 = newarray(idim, DOUBLE);
  error2 = newarray(fouridim, DOUBLE);
  if (error1 == SISL_NULL || error2 == SISL_NULL || maxerr == SISL_NULL)
    goto err101;
  for (i=0; i<fouridim; i++)
    {
      edgeps[i] = MIN(edgeps[i], eeps[(i+idim)%idim]);
      error2[i] = afctol * edgeps[i];
    }
  for (i=0; i<idim; i++)
    error1[i] = afctol*eeps[i];

  /* Perform datareduction on the bilinear interpolant */

  s1345(osurf1, error1, nend, error2, aepsco, iopt, itmax,
	&osurf2, maxerr, &stat);
  if (stat<0) goto error;

  in1 = osurf2->in1;
  in2 = osurf2->in2;

  /* Free surface osurf1 */

  if(osurf1 != SISL_NULL)
    {
      freeSurf(osurf1);
      osurf1 = SISL_NULL;
    }

  /* Piecewise linear interpolant to the reduced
   * bilinear interpolant expressed as a surface
   * of orders ik1 and ik2 */

  /* Second parameter direction */

  s1350(osurf2->ecoef,&(osurf2->et2)[1], in2,
	in1 * idim, ik2, &ocurve1, &stat);
  if (stat<0) goto error;

  newin2 = ocurve1->in;

  /* Transpose result, store new coefficients in
   * array newcoeff */

  if( (newcoeff = newarray(idim * in1 * newin2, DOUBLE)) == SISL_NULL )
    goto err101;
  s6chpar(ocurve1->ecoef, in1, newin2, idim, newcoeff);

  /* First parameter direction */

  s1350(newcoeff, &(osurf2->et1)[1], in1,
	idim*newin2, ik1, &ocurve2, &stat);
  if (stat<0) goto error;
  newin1 = ocurve2->in;

  /* Free surface osurf2 */

  if(osurf2 != SISL_NULL)
    {
      freeSurf(osurf2);
      osurf2 = SISL_NULL;
    }

  /* Transpose back and get coefficients of bilinear
   * approximatoin surface of orders ik1 and ik2     */

  newcoeff = increasearray(newcoeff,
			   idim * newin1 * newin2, DOUBLE);
  if (newcoeff == SISL_NULL) goto err101;
  s6chpar(ocurve2->ecoef, newin2, newin1, idim, newcoeff);

  /* Store results as a surface */

  osurf1 = newSurf(newin1, newin2, ik1, ik2, ocurve2->et,
		   ocurve1->et, newcoeff, 1, idim, 1);

  free(newcoeff); newcoeff = SISL_NULL;

  if (osurf1 == SISL_NULL) goto err101;

  /* Set periodicity flag. */

  osurf1->cuopen_1 = ocurve2->cuopen;
  osurf1->cuopen_2 = ocurve1->cuopen;

  /* Compute tolerance for final datareduction */

  for (i=0; i<fouridim; i++)
    error2[i] = edgeps[i]-error2[i];
  for (i=0; i<idim; i++)
    error1[i] = eeps[i]-maxerr[i];

  /* Perform final datareduction step */

  s1345(osurf1, error1, nend, error2, aepsco, iopt, itmax,
	rs, emxerr, &stat);
  if (stat<0) goto error;

  /* Compute total (and final) error */

  for(i=0; i<idim; i++)
    emxerr[i] += maxerr[i];

  /* Success */

  *jstat = 0;
  goto out;

  /* Empty array. */

 err101:
  *jstat = -101;
  s6err("s1346",*jstat,kpos);
  goto out;

  /* Error in input */

 err103:
  *jstat = -103;
  s6err("s1346",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */

 error:
  *jstat = stat;
  s6err("s1346",*jstat,kpos);
  goto out;

  /* Exit */

 out:
  /* Free SISL-curves allocated in this routine */

  if(ocurve1 != SISL_NULL) freeCurve(ocurve1);
  if(ocurve2 != SISL_NULL) freeCurve(ocurve2);

  /* Free SISL-surfaces allocated in this routine */

  if(osurf1 != SISL_NULL) freeSurf(osurf1);

  /* Free arrays */

  if(error1 != SISL_NULL) freearray(error1);
  if(error2 != SISL_NULL) freearray(error2);
  if(maxerr != SISL_NULL) freearray(maxerr);

  if (ipar != 3)
    {
      freearray(par1);
      freearray(par2);
    }

  return;
}
