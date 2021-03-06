/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: s1334.c,v 1.2 2001-03-19 15:58:46 afr Exp $
 *
 */


#define S1334

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1334(double epoint[],int inbpnt,int idim,double nptyp[],
      int icnsta,int icnend,int iopen,int ik,double astpar,
      double *cendpar,SISLCurve **rc,double **gpar,int *jnbpar,int *jstat)
#else
void s1334(epoint,inbpnt,idim,nptyp,icnsta,icnend,iopen,ik,astpar,
           cendpar,rc,gpar,jnbpar,jstat)
     double epoint[];
     int    inbpnt;
     int    idim;
     double    nptyp[];
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
*                  0 : Closed curve.
*                  1 : Open curve.
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
*****************************************************************
*/
{
  int kpos = 0;
  int kstat = 0;
  int *ltype = SISL_NULL;		/* The kind of interpolation conditions. */
  int knpt;			/* Number of acepted interpolation conditions. */
  double *lcond = SISL_NULL;		/* The number of acepted interpolation conditions. */

  int *ityp = SISL_NULL;
  int ki;

  /* make compatible to old use of s1604 */

  if (iopen==SISL_CRV_CLOSED) iopen = SISL_CRV_PERIODIC;

  /* a new version with input-iopen == rc->cuopen
     should be made, with a name different from any
     other.
     NOTE: There is an error in this function when iopen = 0
           qc as input to s1713 (and s1750) then has wrong flag !!*/

  ityp = newarray (inbpnt, INT);
  if (ityp == SISL_NULL)
    goto err101;

  for (ki=0; ki < inbpnt; ki++) ityp[ki] = (int)nptyp[ki];

  /* <- guen & ujk Wed Jul  1 17:06:46 MESZ 1992 */

  *jstat = 0;

  /* Transform interpolation conditions. */

/* -> guen & ujk Wed Jul  1 17:06:46 MESZ 1992 */
/*  s1906 (epoint, nptyp, icnsta, icnend, inbpnt, idim, &lcond, */
/*	 &ltype, &knpt, &kstat);				*/
  s1906 (epoint, ityp, icnsta, icnend, inbpnt, idim, &lcond,
	 &ltype, &knpt, &kstat);
/* <- guen & ujk Wed Jul  1 17:06:46 MESZ 1992 */

  if (kstat < 0)
    goto error;

  /* Interpolate. */

  s1901(s1909, s1902, lcond, ltype, knpt, astpar, ik, idim, iopen,
	cendpar, rc, gpar, jnbpar, &kstat);  
  if (kstat < 0) goto error;

  *jstat = 0;
  goto out;
  
  /* Error in lower level routine. */
 error: *jstat = kstat;
  s6err("s1334",*jstat,kpos);
  goto out;

  /* allocation error */
 err101: *jstat = -101;
  s6err("s1334",*jstat,kpos);
  goto out;
  
 out:
  if (ltype != SISL_NULL)
    freearray (ltype);

/* -> added, guen & ujk Wed Jul  1 18:48:27 MESZ 1992 */
  if (ityp != SISL_NULL)
    freearray (ityp);
/* <- added, guen & ujk Wed Jul  1 18:48:27 MESZ 1992 */

  if (lcond != SISL_NULL)
    freearray (lcond);
    
  return;
}
