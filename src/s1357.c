/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"


#define S1357

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s1357(double epoint[],int inbpnt,int idim,int ntype[],double epar[],
	   int icnsta,int icnend,int iopen,int ik,double astpar,
	   double *cendpar,SISLCurve **rc,double **gpar,int *jnbpar,int *jstat)
#else
void s1357(epoint,inbpnt,idim,ntype,epar,icnsta,icnend,iopen,ik,astpar,
           cendpar,rc,gpar,jnbpar,jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    ntype[];
     double epar[];
     int    icnsta;
     int    icnend;
     int    iopen;
     int    ik;
     double astpar;
     double *cendpar;
     SISLCurve  **rc;
     double **gpar;
     int    *jnbpar;
     int    *jstat;
#endif
/*
*************************************************************************
*
* Purpose: To calculate a B-spline curve interpolating a set of points.
*          The points can be assigned a tangent (derivative).
*          The curve can be closed or open. If end-conditions are
*          conflicting, the condition closed curve rules out other
*          end conditions.
*          The parametrization is given by the array Epar.
*
* Input:
*        Epoint - Array (length idim*inbpnt) containing the points/
*                 derivatives to be interpolated.
*        Inbpnt - No. of points/derivatives in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        ntype  - Array (length inbpnt) containing type indicator for
*                 points/derivatives/second-derivatives:
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Derivative to next point.
*                  4 - Derivative to prior point.
*                ( 5 - Second derivative to next point. )
*                ( 6 - Second derivative to prior point. )
*                 13 - Start-point of tangent to next point.
*                 14 - End-point of tangent to prior  point.
*        Epar   - Array containing the wanted parametrization. Only parameter
*                 values corresponding to position points are given.
*                 For closed curves, one additional parameter value
*                 must be spesified. The last entry contains
*                 the parametrization of the repeted start point.
*                 (if the endpoint is equal to the startpoint of
*                 the interpolation the lenght of the array should
*                 be equal to inpt1 also in the closed case).
*        Icnsta - Additional condition at the start of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at start.
*        Icnend - Additional condition at the end of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at end.
*        Iopen  - Flag telling if the curve should be open or closed:
*                  1 : The curve should be open.
*                  0 : The curve should be closed.
*                 -1 : The curve should be closed and periodic.
*        Ik     - The order of the B-spline curve to be produced.
*        Astpar - Parameter-value to be used at the start of the curve.
*
* Output:
*        Jstat  - status variable:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*        cendpar - Parameter-value used at the end of the curve.
*        Rc     - Pointer to output-curve.
*        Gpar   - Pointer to the parameter-values of the points in
*                 the curve. Represented only once, although derivatives
*                 and second-derivatives will have the same parameter-
*                 value as the points.
*        Jnbpar - No. of different parameter-values.
*
* Method: First the parametrization of the curve and the type
*         specification of points/derivatives are checked and/or
*         corrected. Then the knots are calculated, and the
*         interpolation is performed.
*
* Calls: s1907,s1908,s1902,s1891,s1713,s1750,s6err.
*
* Written by: Trond Vidar Stensby, SI, 1991-07
* The Fortran routine, p19539, is written by: T. Dokken  SI.
* Bug fix:    Michael Floater, 82.9.93. The construction of gpar and *jnbpar
*                  at the end of the routine was incorrect.
* REWISED BY: Vibeke Skytt, 03.94. This routine corresponds to s1358,
*                                  but differ in the use of the parameter
*                                  iopen and the input array ntype is of
*                                  type int.
*****************************************************************
*/
{
  int kpos = 0;
  SISLCurve *qc = NULL;		/* Temporary SISLCurves.                    */
  SISLCurve *qc2 = NULL;
  int *ltype = NULL;		/* Type of interpolation condition 
				   (temporary)                              */
  int *ltype2 = NULL;		/* Type of interpolation condition (finial) */
  int *sder = NULL;		/* Vector of derivative indicators.         */
  int knpt;			/* Number of interpolation conditions.      */
  int kordr;			/* Local order.                             */
  int kstat;			/* Status variable.                         */
  int kn;			/* Number of coefficients of B-spline curve.*/
  int ki;			/* Loop control variable.                   */
  int kright = 1;		/* One equation system to solve in 
				   interpolation.                           */
  int knlr = 0;			/* Indicates shape of interpolation matrix. */
  int knrc = 0;			/* Indicates shape of interpolation matrix. */
  int kopen;                    /* Local open/closed parameter. Closed,
				   non-periodic is treated as an open curve.*/
  double *lpar = NULL;		/* Parameter values. (temporary)            */
  double *lcond = NULL;		/* Interpolation conditions. (temporary)    */
  double *sknot = NULL;		/* Knot vector.                             */
  double *spar = NULL;		/* Parameter valued. (finial)               */
  double *scond = NULL;		/* Interpolation conditions. (finial)       */
  double *scoef = NULL;		/* Coefficients of curve.                   */

  *jstat = 0;

  /* Set local open/closed parameter. */
  
  kopen = (iopen == SISL_CRV_PERIODIC) ? 0 : 1;
  
  /* Transform interpolation conditions. */

  s1907 (epoint, ntype, epar, iopen, icnsta, icnend, inbpnt,
	 idim, &lcond, &ltype, &lpar,&knpt, &kstat);
  if (kstat < 0) goto error;

  /* Test interpolation conditions, and adjust the input conditions
     if necessary.  */

  s1908 (lcond, ltype, lpar, knpt, ik, idim, iopen, &scond, &ltype2,
	 &spar, &knpt, &kstat);
  if (kstat < 0) goto error;

  /* Allocate scratch for derivative indicator. */

  sder = newarray (knpt, INT);
  if (sder == NULL) goto err101;

  for (ki = 0; ki < knpt; ki++)
    sder[ki] = abs (ltype2[ki]);

  kordr = MIN (ik, knpt);

  if (iopen == SISL_CRV_PERIODIC)
    {
      knlr = kordr / 2;
      knrc = kordr - knlr - 1;
      knpt--;
    }

  /* Produce knot vector. */

  s1902 (lpar, knpt, kordr, kopen, &sknot, &kstat);
  if (kstat < 0) goto error;

  /* Perform interpolation.  */

  s1891 (spar, scond, idim, knpt, kright, sder, kopen, sknot,
	 &scoef, &kn, kordr, knlr, knrc, &kstat);
  if (kstat < 0) goto error;

  /* Express the curve as a curve object.  */

  qc = newCurve (kn, kordr, sknot, scoef, 1, idim, 1);
  if (qc == NULL) goto err101;
  qc->cuopen = iopen;
  
  if (kordr < ik)
    {
      /* The order of the curve is less than expected. Increase the order. */

      qc2 = NULL;
      s1750 (qc, ik, &qc2, &kstat);
      if (kstat < 0) goto error;
	
      if (qc != NULL) freeCurve (qc);
      qc = qc2;
    }

  /* Set open/closed parameter of curve. */

  qc->cuopen = iopen;

  /* Set end of parameter interval.  */
  
  *cendpar = *(qc->et + qc->in);

  /* Interpolation performed. */

  /* Find distinct parameter values. */

  *gpar = lpar;

  /* Bug fix, MSF, 28.9.93. Change
      *jnbpar = 1;
     to:  */

  *jnbpar = 0;

  for (ki = 1; ki<knpt; ki++)
    {
      if (spar[ki - 1] < spar[ki])
	 
      /* Bug fix. MSF, 28.9.93 Change (*gpar)[(*jnbpar)++] = spar[ki]; to: */

	(*gpar)[(*jnbpar)++] = spar[ki-1];
    }

  /* Bug fix. MSF, 28.9.93 Change (*gpar)[(*jnbpar)++] = spar[ki]; to: */

  (*gpar)[(*jnbpar)++] = spar[ki-1];
  
  *rc = qc;
  goto out;


  /* Error in scratch allocation.  */

  err101:
    *jstat = -101;
    s6err ("s1357", *jstat, kpos);
    goto out;

  /* Error in lower level routine. */

  error:
    *jstat = kstat;
    s6err ("s1357", *jstat, kpos);
    goto out;

  out:
  /* Free scratch occupied by local arrays. */

    if (scond != NULL)    freearray (scond);
    if (scoef != NULL)    freearray (scoef);
    if (sder != NULL)     freearray (sder);
    if (ltype != NULL)    freearray (ltype);
    if (ltype2 != NULL)   freearray (ltype2);
    if (lcond != NULL)    freearray (lcond);
    if (sknot != NULL)    freearray (sknot);
    if (spar != NULL)     freearray (spar);

    return;
}

