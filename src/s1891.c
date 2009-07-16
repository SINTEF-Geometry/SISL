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
 * $Id: s1891.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1891

#include "sislP.h"
#define MAX_SIZE  50


#if defined(SISLNEEDPROTOTYPES)
void
s1891 (double etau[], double epoint[], int idim, int inbpnt, int iright,
       int eder[], int iopen, double et[], double *ebcoef[], int *in,
       int ik, int inlr, int inrc, int *jstat)

#else
void
s1891 (etau, epoint, idim, inbpnt, iright, eder, iopen, et, ebcoef,
       in, ik, inlr, inrc, jstat)
     double etau[];
     double epoint[];
     int idim;
     int inbpnt;
     int iright;
     int eder[];
     int iopen;
     double et[];
     double *ebcoef[];
     int *in;
     int ik;
     int inlr;
     int inrc;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : 	To compute the B-spline coefficients required for
*		interpolation of epoint with a B-spline.
*
* INPUT      : 	etau	- Parameter values (i.e. parametrization of
*			  points in epoint) calculated in s1890.
*		epoint	- Point array. Contains points/derivatives.
*		idim	- Dimension of the space in which the points
*			  lie.
*		inbpnt	- Number of points/derivatives
*		iright	- Number of right hand sides of dimension
*			  (idim,inbpnt)
*		eder	- Parametrization of derivatives in epoint
*			  (calculated in s1890)
*		iopen	- Open or closed curve.
*		et	- Knot vector of length (in+ik)
*			  Open curve:   in = inbpnt
*			  Closed curve: in = inbpnt+ik-1
*		ik	- Order of B-spline basis to be used
*		inlr	- Indicates shape of matrix required for B-spline
*			  interpolation (see subroutines s1925 & s1926).
*		inrc	- Indicates shape of matrix required for B-spline
*			  interpolation (see subroutines s1925 & s1926).
*
* OUTPUT     : 	ebcoef	- Guiding points for B-spline curve.
*			  Dimension (1:idim*in*iright)
*		in	- Number of vertices in the B-spline curve
*			  produced.  Open curve:   in=inbpnt
*				     Closed curve: in=inbpnt+ik-1.
*               jstat	- Status variable:
*			  	< 0  	: error
*				> 0	: warning
*				= 0	: OK.
*
* METHOD     : 	See: 	Carl De Boor: "A practical guide to splines"
*
* REFERENCES :	Fortran version:
*		E. Aarn[s, CP, 1979-01
*
* CALLS      :  s1925,s6err.
*
* WRITTEN BY : 	Christophe R. Birkeland 1991-06
*
*********************************************************************
*/
{
  int kstat = 0;
  int kpos = 0;			/* Position of error			*/
  int ii;			/* Loop control parameter		*/
  int limit1, limit2;		/* Loop parameters			*/
  int kj, kl;
  int kdum, stop;
  int nur;			/* Number of upper rows in W		*/
  int inlx;			/* Equal to inlr if inlr>0, else=1 	*/
  int inrx;			/* Equal to inrc if inrc>0, else=1	*/
  int edarray[MAX_SIZE];        /* Array for ed below                   */
  int alloc_needed=FALSE;
  int *ed = SISL_NULL;		/* Arrays defining elements of W	*/
  double *ewarray=SISL_NULL;         /* Array for ew1, ew2 and ew3           */
  double *ew1 = SISL_NULL;		/* See subroutine s1926			*/
  double *ew2 = SISL_NULL;
  double *ew3 = SISL_NULL;

  *jstat = 0;


  /* Test if legal input. */

  if (ik < 1 || idim < 1) goto err112;

  /* Indicate dimension of B-spline. */

  *in = inbpnt;
  if (iopen != SISL_CRV_OPEN)    *in +=ik - 1;

  *ebcoef = new0array (*in *idim * iright, DOUBLE);
  if (*ebcoef == SISL_NULL) goto err101;

  if ((nur = inbpnt - inlr) > MAX_SIZE)
    alloc_needed = TRUE;

  /* Allocate arrays ew1, ew2, ew3, ed. */

  inlx = MAX (1, inlr);
  inrx = MAX (1, inrc);
  limit1 = (ik * nur) + (inrx * nur) + (inlx * inbpnt);
  
  if ((ewarray = new0array(limit1 + 1,DOUBLE)) == SISL_NULL) goto err101;
  
  ew1 = ewarray;
  ew2 = ew1 + (ik * nur);
  ew3 = ew2 + (inrx * nur);

  if (alloc_needed)
    {
       if ((ed = new0array(nur,INT)) == SISL_NULL)
	 goto err101;
    }
  else
    ed = edarray;
  
  s1925 (etau, epoint, inbpnt, eder, et, *ebcoef,*in, ik, iright, 
	 idim, ew1, nur, ed, ew2, inrc, ew3, inlr, &kstat);
  if (kstat < 0) goto error;

  /* For closed B-spline curves we have:
   * ebcoef(i) = ebcoef(i+inbpnt) ; i=1,...,ik-1. */

  if (iopen != SISL_CRV_OPEN)
    {
      stop = ik - 1;
      for (kl = 0; kl < iright; kl++)
	{
	  kdum = *in *kl;
	  for (kj = 0; kj < stop; kj++)
	    {
	      limit2 = (kj + kdum) * idim;
	      limit1 = inbpnt * idim + limit2;
	      for (ii = 0; ii < idim; ii++)
		(*ebcoef)[limit1 + ii] = (*ebcoef)[limit2 + ii];
	    }
	}
    }

  goto out;

  /* Error in lower level routine */

  error:
    *jstat = kstat;
    s6err ("s1891", *jstat, kpos);
    goto out;

  /* Error in array allocations */

  err101:
    *jstat = -101;
    s6err ("s1891", *jstat, kpos);
    goto out;

  /* Error in description of B-spline */

  err112:
    *jstat = -112;
    s6err ("s1891", *jstat, kpos);
    goto out;

  out:
    if (alloc_needed)    freearray (ed);
    if (ewarray)         freearray (ewarray);
    return;
}
#undef MAX_SIZE
