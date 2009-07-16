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
 * $Id: s1504.c,v 1.3 2006-11-17 20:55:25 vsk Exp $
 *
 */


#define S1504
#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
s1504(double *et,int ik,int in,
	   double *ax,int im,int ider,double ebder[],int ileft[],int *jstat)
#else
void s1504(et,ik,in,ax,im,ider,ebder,ileft,jstat)
     double *et;
     int    ik;
     int    in;
     double *ax;
     int    im;
     int    ider;
     double ebder[];
     int ileft[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the value and ider first derivatives of the
*              ik (possibly) nonzero B-splines associated with the knot
*              vector et at the array of points ax[0],...,ax[im-1].
*
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax     - The array of values at which the non-zero B-splines
*                       and derivatives are to be computed.
*              ider   - The number of derivatives to be computed.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*
* OUTPUT     : ebder  - Triple array of dimension [(ider+1)*ik*im] containing
*                       values of the ik nonzero B-splines and their
*                       derivatives at the points ax[0],...,ax[im-1].
*                       These numbers are stored in the following order:
*                       First the (ider+1) derivatives of the first nonzero
*                       B-spline at ax[0]. Then the (ider+1) derivatives of
*                       the second nonzero B-spline at ax[0], etc.
*                       Later we repeat for ax[1],... etc.
* OUTPUT     : ileft  - Array of dimension im containing the index of
*                       the first non-zero B-spline for each value in ax.
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : Call 1220 im times.
*
* CALLS      : s1220.
*
* WRITTEN BY : Michael Floater, SINTEF, May, 1998.
*
*********************************************************************
*/
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of error.                          */
  int kleft=0;        /* Local version of ileft.                         */
  int j,k,kk;         /* Control variables in for loops and for stepping
			 through arrays.                                 */
  double *eder = SISL_NULL;  /* B-spline evaluationas at a single value. */
  int size;           /* (ider+1) * ik.                                  */

  if (ider < 0) goto err178;

  size = (ider + 1) * ik;
  eder = newarray(size,double);
  if (eder == SISL_NULL) goto err101;

  kk = 0;
  for(k=0; k<im; k++)
  {
    s1220(et,ik,in,&kleft,ax[k],ider,eder,&kstat);
    if (kstat < 0) goto error;

    ileft[k] = kleft;
    for(j=0; j<size; j++)
    {
      ebder[kk+j] = eder[j];
    }
    kk += size;
  }

  *jstat = 0;
  goto out;

  /* Not enough memory. */
 err101: *jstat = -101;
  s6err("s1504",*jstat,kpos);
  goto out;

  /* Illegal derivative requested. */
 err178: *jstat = -178;
  s6err("s1504",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("s1504",*jstat,kpos);
  goto out;

 out: 
  if (eder != SISL_NULL)
      freearray(eder);
  return;
}
