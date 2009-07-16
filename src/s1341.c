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
 * $Id: s1341.c,v 1.2 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1341

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1341(double ep[],int im,int idim,int ipar,double epar[],
	   double eeps[],int ilend,int irend,double afctol,double aepsco,
	   int itmax,int ik,SISLCurve **rc,double emxerr[],int *jstat)
#else
void s1341(ep,im,idim,ipar,epar,eeps,ilend,irend,afctol,aepsco,
           itmax,ik,rc,emxerr,jstat)
     double ep[];
     int    im;
     int    idim;
     int    ipar;
     double epar[];
     double eeps[];
     int    ilend;
     int    irend;
     double afctol;
     double aepsco;
     int    itmax;
     int    ik;
     SISLCurve  **rc;
     double emxerr[];
     int    *jstat;
#endif
/*
************************************************************
*
* Purpose: To compute a spline-approximation to the data given by the
*          points ep, and represent it as a B-spline curve with
*          parameterization determined by the parameter ipar.
*          The approximation is determined by first forming the piecewise
*          linear interpolant to the data, and then performing knot
*          removal on this initial approximation.
*
* Input :
*        Ep     - Array (length idim*im) containing the points to
*                 be approximated.
*        Im     - The no. of data points.
*        Idim   - The dimension of the euclidean space in which the data
*                 points lie, i.e. the number of components of each data point.
*        Ipar   - Flag indicating the type of parameterization to be used:
*                  = 1 : Paramterize by accumulated cord length.
*                        (Arc length parametrization for the piecewise
*                        linear interpolant.)
*                  = 2 : Uniform parameterization.
*                  = 3 : Parametrization given by epar.
*                 If ipar<1 or ipar>3, it will be set to 1.
*        Epar   - Array (length im) containing a parametrization
*                 of the given data.
*        Eeps   - Array (length idim) containing the tolerance to be
*                 used during the data reduction stage. The final
*                 approximation to the data will deviate less than eeps
*                 from the piecewise linear interpolant in each of the
*                 idim components.
*        Ilend  - The no. of derivatives that are not allowed to change
*                 at the left end of the curve.
*                 The (0 - (ilend-1)) derivatives will be kept fixed.
*                 If ilend <0, this routine will set it to 0.
*                 If ilend<ik, this routine will set it to ik.
*        Irend  - The no. of derivatives that are not allowed to change
*                 at the right end of the curve.
*                 The (0 - (irend-1)) derivatives will be kept fixed.
*                 If irend <0, this routine will set it to 0.
*                 If irend<ik, this routine will set it to ik.
*        Afctol = Number indicating how the tolerance is to be shared
*                 between the two data reduction stages. For the linear
*                 reduction, a tolerance of afctol*eeps will be used,
*                 while a tolerance of (1-afctol)*eeps will be used
*                 during the final data reduction. (Similarly for edgeps.)
*        Aepsco - Two numbers differing by a relative amount less than
*                 aepsco will in some cases be considered equal.
*        Itmax  - Max. no. of iterations in the data-reduction routine.
*        Ik     - The polynomial order of the approximation.
*
* Output:
*        Jstat  - Output staus:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > o : Warning.
*        Rc     - Pointer to curve.
*        Emxerr - Array (length idim) (allocated outside this routine.)
*                 containing for each component an upper bound on the
*                 max. deviation of the final approximation from the
*                 initial piecewise linear interpolant.
*
* Method: First the piecewise linear interpoolant is computed by a call
*         to s1350 or s1351, and then knots are removed from this
*         interpolant. The degree of the resulting linear spline is then
*         raised to ik, and a second data-reduction is performed.
*
* The main routine,s1340, is written by Knut M|rken,  Si.
* Written by: C.R.Birkeland  Si  Oslo,Norway April 1993.
********************************************************************
*/
{
  double *maxerr = SISL_NULL;    /* Arrays used to store error estimates */
  double *error1 = SISL_NULL;
  int i;
  int stat = 0;             /* Loop control variables               */
  int kpos = 0;
  SISLCurve *ocurve = SISL_NULL; /* Local spline curve                   */
  
  /* Check Input */
  
  if (im < 2 || idim < 1) goto err103;
  if (ipar < 1 || ipar > 3) ipar = 1;

  /* Set default value for afctol */

  if(afctol < 0 || afctol > 1) afctol =0;

  /* Piecewise linear interpolant to the data */

  if (ipar == 3)
    s1350(ep, epar, im, idim, 2, &ocurve, &stat);
  else
    s1351(ep, ipar, im, idim, 2, &ocurve, &stat);
  if (stat<0) goto error;

  /* Compute tolerance for linear reduction */

  error1 = newarray(idim, DOUBLE);
  maxerr = newarray(idim, DOUBLE);
  if (error1 == SISL_NULL || maxerr == SISL_NULL) goto err101;
  for (i=0; i<idim; i++)
    maxerr[i] = eeps[i]*afctol;

  /* Data reduction on piecewise linear interpolant */

  s1340(ocurve, maxerr, ilend, irend, aepsco, itmax, rc,
	error1, &stat);
  if (stat<0) goto error;
  freeCurve(ocurve);

  /* Piecewise linear interpolant to the reduced data 
     (coefficients of linear spline curve rc)
     expressed as B-spline of order ik */

  s1350((*rc)->ecoef, &((*rc)->et)[1], (*rc)->in, 
	idim, ik, &ocurve, &stat);
  if (stat<0) goto error;
  freeCurve(*rc);
  
  /* Compute tolerance for last reduction */

  for (i=0; i<idim; i++)
    maxerr[i] = eeps[i] - error1[i];

  /* Perform data reduction on the hermite interpolant */

  s1340(ocurve, maxerr, ilend, irend, aepsco, itmax, rc,
	emxerr, &stat);
  if (stat<0) goto error;

  /* Set periodicity flag.  */
  
  (*rc)->cuopen = ocurve->cuopen;
  
  /* Compute total error after reduction */

  for (i=0; i<idim; i++)
    emxerr[i] += error1[i];
  
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;

  /* Error in input */

 err103: 
  *jstat = -103;
  s6err("s1341",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */

 error:
  *jstat = stat;
  s6err("s1341",*jstat,kpos);
  goto out;

  /* Exit */

 out:
  if (maxerr != SISL_NULL) freearray(maxerr);
  if (error1 != SISL_NULL) freearray(error1);
  if (ocurve != SISL_NULL) freeCurve(ocurve);
  return;
}
