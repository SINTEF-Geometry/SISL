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
 * $Id:
 *
 */
#define S1542

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1542(SISLCurve *pc1,int m,double x[],double eder[],int *jstat)
#else
void s1542(pc1,m,x,eder,jstat)
     SISLCurve *pc1;
     int      m;
     double   x[];
     double   eder[];
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the curve pointed at by pc1 over a m grid
*              of points (x[i]). Only positions are evaluated.
*              This does not work for in the rational case.
*
* INPUT      : pc1    - Pointer to the curve to evaluate.
*              m      - Number of grid points.
*              x      - Array of parameter values of the grid.
*
* OUTPUT     : eder   - Array where the derivatives of the curve
*                       are placed, dimension
*                         idim * (ider+1) * m.
*                       The sequence is position at point x[0],
*                       followed by the same information at x[1],
*                       etc.
*              jstat  - status messages
*                          = 0      : Ok.
*                          < 0      : Error.
*
* METHOD     : We call s1540 to pre-evaluate the B-splines then call
*              s1541 to multiply them with the coefficients.
*
*-
* CALLS      : s1540, s1541.
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Nov. 1999.
*********************************************************************
*/
{
  int ider = 0;        /* Only evaluate position.                         */
  int kstat=0;         /* Local status variable.                          */
  int kpos=0;          /* The position of error.                          */
  int in;              /* The number of B-splines accociated with the knot
		 	 vector.                                          */
  int ik;              /* The polynomial order of the curve.              */
  int kdim;            /* The space dimension of the surface. */
  double *ebder=SISL_NULL;  /* Triple array of dimension (ider+1)*ik*m
                         containing dericatives of B-splines. */
  int *ileft=SISL_NULL;     /* Array of dimension m containing the left knots
                         of the B-splines. */
  double *et = SISL_NULL;   /* The knot vector. */

  in = pc1 -> in;
  ik = pc1 -> ik;
  et = pc1 -> et;
  kdim = pc1 -> idim;

  /* Check the input. */
  if (kdim < 1) goto err102;
  if (ik < 1) goto err115;
  if (in < ik) goto err116;
  if (ider < 0) goto err178;

  /* Pre-evaluate B-splines. */
  ebder = newarray((ider+1)*ik*m, DOUBLE);
  if(ebder == SISL_NULL) goto err101;

  ileft = newarray(m,INT);
  if(ileft == SISL_NULL) goto err101;

  s1540(et,ik,in,x,m,ider,ebder,ileft,&kstat);
  if(kstat < 0) goto error;

  /* Multiply out with the coefficients. */

  s1541(pc1,m,ebder,ileft,eder,&kstat);
  if(kstat < 0) goto error;

  /* Free memory. */
  if(ebder != SISL_NULL) freearray(ebder);
  if(ileft != SISL_NULL) freearray(ileft);

  *jstat = 0;
  goto out;

  /* Not enough memory. */
 err101: *jstat = -101;
  s6err("s1542",*jstat,kpos);
  goto out;

  /* kdim less than 1. */
 err102: *jstat = -102;
  s6err("s1542",*jstat,kpos);
  goto out;

  /* Polynomial order less than 1. */
 err115: *jstat = -115;
  s6err("s1542",*jstat,kpos);
  goto out;

  /* Fewer B-splines than the order. */
 err116: *jstat = -116;
  s6err("s1542",*jstat,kpos);
  goto out;

  /* Illegal derivative requested. */
 err178: *jstat = -178;
  s6err("s1221",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("s1542",*jstat,kpos);
  goto out;

 out:
    return;
}
