/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1901.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1901

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
s1901 (fparamProc fparam, fknotsProc fknots,
       double econd[], int ntype[], int inpt, double astpar, int ik,
       int idim, int iopen, double *cendpar, SISLCurve ** rcurve,
       double **gpar, int *jnbpar, int *jstat)
#else
void
s1901 (fparam, fknots, econd, ntype, inpt, astpar, ik, idim, iopen,
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
  int *ltype = NULL;		/* Type of accepted interpolation conditions.   */
  double *scond = NULL;		/* Array containing interpolation conditions.   */
  double *spar1 = NULL;		/* Parametrization array of interpolation conditions. */
  double *spar2 = NULL;		/* Parametrization array used to make knot vector. */
  double *sknot = NULL;		/* Knot vector of curve.                           */
  double *scoef = NULL;		/* Coefficients of curve.                          */
  int *sder = NULL;		/* Vector of derivative indicators.                */
  SISLCurve *qc = NULL;		/* Interpolation curve.                            */
  SISLCurve *qc2 = NULL;	/* Interpolation curve.                            */

  *jstat = 0;

  /* Test interpolation conditions, and adjust the input conditions
     if necessary.  */

  s1905 (econd, ntype, inpt, ik, idim, iopen, &scond, &ltype, &knpt, &kstat);
  if (kstat < 0)
    goto error;

  /* Set local order.  */

  kordr = MIN (ik, knpt);

  /* Allocate scratch for derivative indicator. */

  if ((sder = newarray (knpt, INT)) == NULL)
    goto err101;

  for (ki = 0; ki < knpt; ki++)
    sder[ki] = (int) fabs ((double) ltype[ki]);

  /* Compute parametrization of point set. */

  (* fparam) (scond, ltype, knpt, idim, iopen, astpar, cendpar, 
	      &spar1, &spar2, &kstat);  
  if (kstat < 0) goto error;

  /* Make knot vector of curve.  */

  if (!(iopen == SISL_CRV_OPEN))
    {
      /* Closed curve. */

      knlr = kordr / 2;
      knrc = kordr - knlr - 1;
      knpt--;
    }

  /* Produce knot vector. */

  (* fknots) (spar1, knpt, kordr, iopen, &sknot, &kstat);
  if (kstat < 0) goto error;

  /* Perform interpolation.  */

  s1891 (spar1, scond, idim, knpt, kright, sder, iopen, sknot,
	 &scoef, &kn, kordr, knlr, knrc, &kstat);
  if (kstat < 0) goto error;

  /* Express the curve as a curve object.  */

  qc = newCurve (kn, kordr, sknot, scoef, 1, idim, 1);
  if (qc == NULL) goto err101;

  qc->cuopen = (iopen == SISL_CRV_OPEN) ? iopen : SISL_CRV_PERIODIC;

  if (iopen == SISL_CRV_CLOSED)
    {
      /* A closed, non-periodic curve is expected. Pick the part of the
	 interpolation curve that has got a full basis.  */

      s1713 (qc, sknot[kordr - 1], sknot[kn], &qc2, &kstat);
      if (kstat < 0) goto error;

      if (qc != NULL) freeCurve (qc);
      qc = qc2;
    }

  if (kordr < ik)
    {
      /* The order of the curve is less than expected. Increase the order. */

      qc2 = NULL;
      s1750 (qc, ik, &qc2, &kstat);
      if (kstat < 0) goto error;

      if (qc != NULL) freeCurve (qc);
      qc = qc2;
    }

  /* Interpolation performed. */

  /* Find distinct parameter values. */

  *gpar = spar1;

  *jnbpar = 1;
  for (ki = 1; spar1[ki] < *cendpar; ki++)
    {
      if (spar1[ki - 1] < spar1[ki])
	(*gpar)[(*jnbpar)++] = spar1[ki];
    }
  (*gpar)[(*jnbpar)++] = spar1[ki];

  *gpar = increasearray (*gpar, *jnbpar, DOUBLE);

  *rcurve = qc;
  goto out;


  /* Error in scratch allocation.  */

err101:
  *jstat = -101;
  s6err ("s1901", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1901", *jstat, kpos);
  goto out;

out:
  /* Free scratch occupied by local arrays. */

  if (spar2 != NULL)
    freearray (spar2);
  if (scond != NULL)
    freearray (scond);
  if (scoef != NULL)
    freearray (scoef);
  if (sknot != NULL)
    freearray (sknot);
  if (sder != NULL)
    freearray (sder);
  if (ltype != NULL)
    freearray (ltype);

  return;
}
