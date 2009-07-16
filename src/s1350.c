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
 * $Id: s1350.c,v 1.2 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1350

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1350(double ep[],double epar[],
	   int im,int idim,int ik,
	   SISLCurve **rc, int *jstat)
#else
void s1350(ep,epar,im,idim,ik,rc,jstat)
     double ep[];
     double epar[];
     int im;
     int idim;
     int ik;	   
     SISLCurve **rc;
     int    *jstat;
#endif
/*
************************************************************
*
* Purpose: To compute the piecewise linear interpolant to a set 
*          of datapoints and express it as a linear combination 
*          of B-splines of order ik using the parametrization
*          given by epar.
*
* Input :
*        Ep     - Array [idim,im] containing the points to
*                 be approximated.
*        Epar   - Array (length im) containing a parametrization
*                 of the given data.
*        Im     - The no. of data points.
*        Idim   - The dimension of the euclidean space in which the data
*                 points lie, i.e. the number of components of each data point.
*        Ik     - The polynomial order of the approximation.
*
* Output:
*        Jstat  - Output staus:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > o : Warning.
*        Rc     - Pointer to curve.
*
* Method: The routine uses the parametrization given by the array
*         epar. The knotvector will have ik-multiple knots at both 
*         ends and (ik-1)-tuple knots at each interior points on the
*         knot vector. This makes it easy to determine the B-spline
*         coefficients of the piecewise linear interpolant.
*
*
* The fortran version was written by Knut M|rken,  Si.
* Written by: C.R.Birkeland  Si  Oslo,Norway April 1993.
********************************************************************
*/
{
  int i, j, k;                   /* Loop index                     */
  int kic, kit, kw1, kw2;        /* Used in calculations           */
  double ts, tw1, tw2;           /* Used in calculations           */
  int kpos = 0;                  /* Indicator of position of error */
  int in;
  int jidim;                     /*  j*idim                        */
  int jidimp1;                   /*  (j+1)*idim                    */
  double *et = SISL_NULL;             /* Array for knotvector           */
  double *ec = SISL_NULL;             /* Array for coefficients         */
  double ikinv;                  /*   1. / ik                      */
  int kclosed;                   /* Used to test if the curve is closed. */

  /* Check Input */
  
  if (im < 2 || idim < 1 || ik < 2) goto err103;

  /* Allocate matrices */

  in = (ik-1)*im + 2 - ik;
  et = newarray(in+ik, DOUBLE);
  ec = newarray(in*idim, DOUBLE);
  if (et==SISL_NULL || ec == SISL_NULL) goto err101;

  /* Perform the one and only division required 
     in this routine */

  ikinv = 1./(ik-1);

  /* Generate first knots and first coefficient */

  for(i=0; i<ik; i++)
      et[i] = epar[0];
  for(i=0; i<idim; i++)
      ec[i] = ep[i];
  
  /* Compute remaining knots and coefficients */

  kic = idim;
  kit = ik;
  for(j=0, jidim=0, jidimp1=idim; j<im-1; j++, jidim+=idim, 
       jidimp1+=idim)
    {
      ts = epar[j+1];
      
      /* Compute coefficients of the B-splines starting
	 at point j and set knots for next point */

      kw1 = ik-1;
      kw2 = 0;
      for (i=1; i<ik; i++)
	{
	  et[kit] = ts;
	  kit++;
	  kw1--;
	  kw2++;
	  tw1 = kw1*ikinv;
	  tw2 = kw2*ikinv;
	  for (k=0; k<idim; k++)	 
	    ec[kic + k] = tw1*ep[jidim + k] + 
	      tw2*ep[jidimp1 + k];
	  kic += idim;
	}
    }

  /* Set last knot */

  et[kit] = ts;
  if ((*rc = newCurve(in,ik,et,ec,1,idim,2)) == SISL_NULL)
        goto err101;

  /* Test if the input data is closed.  */
  
  for (kclosed=1, i=0; i<idim; i++)
     if (DNEQUAL(ep[i], ep[(im-1)*idim+i])) kclosed = 0;
  if (kclosed) (*rc)->cuopen = SISL_CRV_CLOSED;
     
  /* Success */
  
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
    if (et != SISL_NULL) freearray(et);  
    if (ec != SISL_NULL) freearray(ec);
    goto out;

  /* Error in input */

 err103: 
  *jstat = -103;
  s6err("s1350",*jstat,kpos);
  goto out;
  
  /* Exit */

 out:
  return;
}
