/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1356

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s1356(double epoint[],int inbpnt,int idim,int nptyp[],
      int icnsta,int icnend,int iopen,int ik,double astpar,
      double *cendpar,SISLCurve **rc,double **gpar,int *jnbpar,int *jstat)
#else
void s1356(epoint,inbpnt,idim,nptyp,icnsta,icnend,iopen,ik,astpar,
           cendpar,rc,gpar,jnbpar,jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     int    nptyp[];
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
*
* Input:
*        Epoint - Array (length idim*inbpnt) containing the points/
*                 derivatives to be interpolated.
*        Inbpnt - No. of points/derivatives in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        nptyp  - Array (length inbpnt) containing type indicator for
*                 points/derivatives/second-derivatives:
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Derivative to next point.
*                  4 - Derivative to prior point.
*                ( 5 - Second derivative to next point. )
*                ( 6 - Second derivative to prior point. )
*                 13 - Start-point of tangent to next point.
*                 14 - End-point of tangent to prior  point.
*        Icnsta - Additional condition at the start of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at start.
*        Icnend - Additional condition at the end of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at end.
*        Iopen  - Flag telling if the curve should be open or closed:
*                  1 : Open curve.
*                  0 : Closed, non-periodic curve.
*                 -1 : Periodic curve.
*        Ik     - The order of the B-spline curve to be produced.
*        Astpar - Parameter-value to be used at the start of the curve.
*
* Output:
*        cendpar - Parameter-value used at the end of the curve.
*        Rc     - Pointer to output-curve.
*        Gpar   - Pointer to the parameter-values of the points in
*                 the curve. Represented only once, although derivatives
*                 and second-derivatives will have the same parameter-
*                 value as the points.
*        Jnbpar - No. of different parameter-values.
*        Jstat  - status variable:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method: First the parametrization of the curve and the type
*         specification of points/derivatives are checked and/or
*         corrected. Then the knots are calculated, and the
*         interpolation is performed.
*
* Calls: s1906,s1901,s6err.
*
* Written by: A.M. Ytrehus  Si  Oslo, Norway  Sep. 1988.
* The main routine, p19501, is written by: T. Dokken  SI.
* Rewised by: Trond Vidar Stensby, SI, 1991-07
* REWISED BY: Vibeke Skytt, 03.94. This routine corresponds to s1334,
*                                  but differ in the use of the parameter
*                                  iopen and the input array ntype is of
*                                  type int.
*****************************************************************
*/
{
  int kpos = 0;
  int kstat = 0;
  int *ltype = SISL_NULL;		/* The kind of interpolation conditions. */
  int knpt;			/* Number of acepted interpolation conditions. */
  double *lcond = SISL_NULL;		/* The number of acepted interpolation conditions. */

  *jstat = 0;

  /* Transform interpolation conditions. */

  s1906 (epoint, nptyp, icnsta, icnend, inbpnt, idim, &lcond,
	 &ltype, &knpt, &kstat);
  if (kstat < 0)
    goto error;

  /* Interpolate. */

  s1912(s1909, s1902, lcond, ltype, knpt, astpar, ik, idim, iopen,
	cendpar, rc, gpar, jnbpar, &kstat);  
  if (kstat < 0) goto error;

  *jstat = 0;
  goto out;
  
  /* Error in lower level routine. */
 error: *jstat = kstat;
  s6err("s1356",*jstat,kpos);
  goto out;

 out:
  if (ltype != SISL_NULL)
    freearray (ltype);
  if (lcond != SISL_NULL)
    freearray (lcond);
    
  return;
}
