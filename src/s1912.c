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

#define S1912

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
typedef void (*fparamProc)(double[], int[], int, int, int, double, double *, 
			   double *[], double *[], int *);
typedef void (*fknotsProc)(double[], int, int, int, double *[], int *);
#else
typedef void (*fparamProc)();
typedef void (*fknotsProc)();
#endif


#if defined(SISLNEEDPROTOTYPES)
void
   s1912 (fparamProc fparam, fknotsProc fknots,
       double econd[], int ntype[], int inpt, double astpar, int ik,
       int idim, int iopen, double *cendpar, SISLCurve ** rcurve,
       double **gpar, int *jnbpar, int *jstat)
#else
void
   s1912 (fparam, fknots, econd, ntype, inpt, astpar, ik, idim, iopen,
       cendpar, rcurve, gpar, jnbpar, jstat)
     fparamProc fparam;
     fknotsProc fknots;
     double econd[];
     int ntype[];
     int inpt;
     double astpar;
     int ik;
     int idim;
     int iopen;
     double *cendpar;
     SISLCurve **rcurve;
     double **gpar;
     int *jnbpar;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Compute a B-spline curve interpolating a set of points.
*              The points can be assigned derivative conditions. The
*              curve can be open, closed, or closed and periodic.
*
* INPUT      : fparam - Routine computing parametrization of point set.
*	       fknots - Routine coputing knot vector of point set.
*              econd  - Array of interpolation conditions. Dimension
*                       is inpt*idim.
*              ntype  - Array containing kind of condition. Dimension
*                       is inpt.
*                       =  0 : A point is given.
*                       =  d : The d'th derivatative condition to the
*                              previous point is given.
*                       = -d : The d'th derivatative condition to the
*                              next point is given.
*              inpt   - Number of interpolation conditions.
*              astpar - Start parameter of parametrization.
*              ik     - Order of interpolating curve.
*              idim   - Dimension of geometry space.
*              iopen - Indicates if the curve is to be open, closed or
*                       periodic.
*
* OUTPUT     : cendpar - End parameter of parametrization.
*              rcurve  - Interpolating curve.
*	       gpar    - The distinct parameter values.
*	       jnbpar  - Number of distinct parameter values.
*              jstat   - status messages
*                        = 1      : Specified parametrization method
*                                   replaced by cord length parametrization.
*                                         = 0      : ok
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      :	s1905, s1891, s1713, s1750
*
* WRITTEN BY : Vibeke Skytt, SI, 91-04.
* REVISED BY : Trond Vidar Stensby, SI, 91-07
* REWISED BY : Vibeke Skytt, 94-03. Changed the concept of closed, 
*                                   non-periodic.
* REWISED BY : Johannes Kaasa, 95-11. Fixed error in output of the parameter
*              values (included the last parameter value).
*
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                             */
  int kpos = 0;
  int ki;			/* Counter.                                     */
  int knpt;			/* Number of accepted interpolation conditions. */
  int kn;			/* Number of coefficients of B-spline curve.    */
  int kordr;			/* Local order of curve.                        */
  int kright = 1;		/* One equation system to solve in interpolation. */
  int knlr = 0;			/* Indicates shape of interpolation matrix.     */
  int knrc = 0;			/* Indicates shape of interpolation matrix.     */
  int kopen;                    /* Local open/closed parameter. Closed,
				   non-periodic is treated as an open curve.*/
  int *ltype = SISL_NULL;		/* Type of accepted interpolation conditions.   */
  double *scond = SISL_NULL;		/* Array containing interpolation conditions.   */
  double *spar1 = SISL_NULL;		/* Parametrization array of interpolation conditions. */
  double *spar2 = SISL_NULL;		/* Parametrization array used to make knot vector. */
  double *sknot = SISL_NULL;		/* Knot vector of curve.                           */
  double *scoef = SISL_NULL;		/* Coefficients of curve.                          */
  int *sder = SISL_NULL;		/* Vector of derivative indicators.                */
  SISLCurve *qc = SISL_NULL;		/* Interpolation curve.                            */
  SISLCurve *qc2 = SISL_NULL;	/* Interpolation curve.                            */

  *jstat = 0;

  /* Set local open/closed parameter. */
  
  kopen = (iopen == SISL_CRV_PERIODIC) ? 0 : 1;
  
  /* Test interpolation conditions, and adjust the input conditions
     if necessary.  */

  s1905 (econd, ntype, inpt, ik, idim, iopen, &scond, &ltype, &knpt, &kstat);
  if (kstat < 0)
    goto error;

  /* Set local order.  */

  kordr = MIN (ik, knpt);

  /* Allocate scratch for derivative indicator. */

  if ((sder = newarray (knpt, INT)) == SISL_NULL)
    goto err101;

  for (ki = 0; ki < knpt; ki++)
    sder[ki] = (int) fabs ((double) ltype[ki]);

  /* Compute parametrization of point set. */

  (* fparam) (scond, ltype, knpt, idim, kopen, astpar, cendpar, 
	      &spar1, &spar2, &kstat);  
  if (kstat < 0) goto error;

  /* Make knot vector of curve.  */

  if (iopen == SISL_CRV_PERIODIC)
    {
      /* Closed, periodic curve. */

      knlr = kordr / 2;
      knrc = kordr - knlr - 1;
      knpt--;
    }

  /* Produce knot vector. */

  (* fknots) (spar1, knpt, kordr, kopen, &sknot, &kstat);
  if (kstat < 0) goto error;

  /* Perform interpolation.  */

  s1891 (spar1, scond, idim, knpt, kright, sder, kopen, sknot,
	 &scoef, &kn, kordr, knlr, knrc, &kstat);
  if (kstat < 0) goto error;

  /* Express the curve as a curve object.  */

  qc = newCurve (kn, kordr, sknot, scoef, 1, idim, 1);
  if (qc == SISL_NULL) goto err101;

  qc->cuopen = iopen;

  if (kordr < ik)
    {
      /* The order of the curve is less than expected. Increase the order. */

      qc2 = SISL_NULL;
      s1750 (qc, ik, &qc2, &kstat);
      if (kstat < 0) goto error;

      if (qc != SISL_NULL) freeCurve (qc);
      qc = qc2;
    }

  /* Interpolation performed. */

  /* Find distinct parameter values. */

  *gpar = spar1;
  *jnbpar = 0;

  for (ki = 1; spar1[ki] < *cendpar; ki++)
    {
      if (spar1[ki - 1] < spar1[ki])
	(*gpar)[(*jnbpar)++] = spar1[ki-1];
    }
  (*gpar)[(*jnbpar)++] = spar1[ki-1];
  (*gpar)[(*jnbpar)++] = spar1[ki];

  *gpar = increasearray (*gpar, *jnbpar, DOUBLE);

  *rcurve = qc;
  goto out;


  /* Error in scratch allocation.  */

err101:
  *jstat = -101;
  s6err ("s1912", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1912", *jstat, kpos);
  goto out;

out:
  /* Free scratch occupied by local arrays. */

  if (spar2 != SISL_NULL)
    freearray (spar2);
  if (scond != SISL_NULL)
    freearray (scond);
  if (scoef != SISL_NULL)
    freearray (scoef);
  if (sknot != SISL_NULL)
    freearray (sknot);
  if (sder != SISL_NULL)
    freearray (sder);
  if (ltype != SISL_NULL)
    freearray (ltype);

  return;
}
