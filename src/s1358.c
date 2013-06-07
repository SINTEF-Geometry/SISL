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
 * $Id: s1358.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1358

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1358(double epoint[],int inbpnt,int idim,double ntype[],double epar[],
	   int icnsta,int icnend,int iopen,int ik,double astpar,
	   double *cendpar,SISLCurve **rc,double **gpar,int *jnbpar,int *jstat)
#else
void s1358(epoint,inbpnt,idim,ntype,epar,icnsta,icnend,iopen,ik,astpar,
           cendpar,rc,gpar,jnbpar,jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     double    ntype[];
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
*****************************************************************
*/
{
  int kpos = 0;
  SISLCurve *qc = SISL_NULL;		/* Temporary SISLCurves.                    */
  SISLCurve *qc2 = SISL_NULL;
  int *ltype = SISL_NULL;		/* Type of interpolation condition 
				   (temporary)                              */
  int *ltype2 = SISL_NULL;		/* Type of interpolation condition (finial) */
  int *sder = SISL_NULL;		/* Vector of derivative indicators.         */
  int knpt;			/* Number of interpolation conditions.      */
  int kordr;			/* Local order.                             */
  int kstat;			/* Status variable.                         */
  int kn;			/* Number of coefficients of B-spline curve.*/
  int ki;			/* Loop control variable.                   */
  int kright = 1;		/* One equation system to solve in 
				   interpolation.                           */
  int knlr = 0;			/* Indicates shape of interpolation matrix. */
  int knrc = 0;			/* Indicates shape of interpolation matrix. */
  double *lpar = SISL_NULL;		/* Parameter values. (temporary)            */
  double *lcond = SISL_NULL;		/* Interpolation conditions. (temporary)    */
  double *sknot = SISL_NULL;		/* Knot vector.                             */
  double *spar = SISL_NULL;		/* Parameter valued. (finial)               */
  double *scond = SISL_NULL;		/* Interpolation conditions. (finial)       */
  double *scoef = SISL_NULL;		/* Coefficients of curve.                   */
  int *ityp = SISL_NULL;

/* -> new statement guen & ujk Thu Jul  2 14:59:05 MESZ 1992 */

/* make compatible to old use of s1604 */

  if (iopen==SISL_CRV_CLOSED) iopen = SISL_CRV_PERIODIC;

/* a new version with input-iopen == rc->cuopen
   should be made, with a name different from any
   other.
  NOTE: There is an error in this function when iopen = 0
        qc as input to s1713 (and s1750) then has wrong flag !!*/

/* <- new statement guen & ujk Thu Jul  2 14:59:05 MESZ 1992 */

  ityp = newarray (inbpnt, INT);
  if (ityp == SISL_NULL) goto err101;

  for (ki=0; ki < inbpnt; ki++) ityp[ki] = (int)ntype[ki];
  *jstat = 0;

  /* Transform interpolation conditions. */

  s1907 (epoint, ityp, epar, iopen, icnsta, icnend, inbpnt,
	 idim, &lcond, &ltype, &lpar,&knpt, &kstat);
  if (kstat < 0) goto error;

  /* Test interpolation conditions, and adjust the input conditions
     if necessary.  */

  s1908 (lcond, ltype, lpar, knpt, ik, idim, iopen, &scond, &ltype2,
	 &spar, &knpt, &kstat);
  if (kstat < 0) goto error;

  /* Allocate scratch for derivative indicator. */

  sder = newarray (knpt, INT);
  if (sder == SISL_NULL) goto err101;

  for (ki = 0; ki < knpt; ki++)
    sder[ki] = abs (ltype2[ki]);

  kordr = MIN (ik, knpt);

  if (!(iopen == SISL_CRV_OPEN))
    {
      knlr = kordr / 2;
      knrc = kordr - knlr - 1;
      knpt--;
    }

  /* Produce knot vector. */

  s1902 (lpar, knpt, kordr, iopen, &sknot, &kstat);
  if (kstat < 0) goto error;

  /* Perform interpolation.  */

  s1891 (spar, scond, idim, knpt, kright, sder, iopen, sknot,
	 &scoef, &kn, kordr, knlr, knrc, &kstat);
  if (kstat < 0) goto error;

  /* Express the curve as a curve object.  */

  qc = newCurve (kn, kordr, sknot, scoef, 1, idim, 1);
  if (qc == SISL_NULL) goto err101;
  qc->cuopen = (iopen == SISL_CRV_OPEN) ? iopen : SISL_CRV_PERIODIC;
  
  if (iopen == SISL_CRV_CLOSED)
    {
      /* A closed, non-periodic curve is expected. Pick the part of the
         interpolation curve that has got a full basis.  */

      s1713 (qc, sknot[kordr - 1], sknot[kn], &qc2, &kstat);
      if (kstat < 0) goto error;

      if (qc != SISL_NULL) freeCurve (qc);
      qc = qc2;
    }

  if (kordr < ik)
    {
      /* The order of the curve is less than expected. Increase the order. */

      qc2 = SISL_NULL;
      s1750 (qc, ik, &qc2, &kstat);
      if (kstat < 0) goto error;
	
      if (qc != SISL_NULL) freeCurve (qc);
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

   /*called from s1333: s1358(epoint,21,21,only 1,epar,0,0,1,3,,
           cendpar,rc,gpar,jnbpar,jstat) */
  /* Error: iopen = 1, order=3, hp has a case where the size of gpar is 21 and
     when leaving this loop *jnbpar (and knpt and ki) is also 21 !!!!!!!!!!!!? */
  /* The above error is now fixed, MSF.
     The following chunk of code merely removes double knots
     from the knot vector spar (length knpt) producing
     both the number *jnbpar of different knots and
     an array gpar with the different knots of length *jnbpar. */
 
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
    s6err ("s1358", *jstat, kpos);
    goto out;

  /* Error in lower level routine. */

  error:
    *jstat = kstat;
    s6err ("s1358", *jstat, kpos);
    goto out;

  out:
  /* Free scratch occupied by local arrays. */

    if (scond != SISL_NULL)    freearray (scond);
    if (scoef != SISL_NULL)    freearray (scoef);
    if (sder != SISL_NULL)     freearray (sder);
    if (ltype != SISL_NULL)    freearray (ltype);

/* -> added, guen & ujk Wed Jul  1 18:48:27 MESZ 1992 */
    if (ityp != SISL_NULL)     freearray (ityp);
/* <- added, guen & ujk Wed Jul  1 18:48:27 MESZ 1992 */

    if (ltype2 != SISL_NULL)   freearray (ltype2);
    if (lcond != SISL_NULL)    freearray (lcond);
    if (sknot != SISL_NULL)    freearray (sknot);
    if (spar != SISL_NULL)     freearray (spar);

    return;
}

