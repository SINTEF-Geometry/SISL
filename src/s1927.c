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
 * $Id: s1927.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1927

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1927 (double *w1, int nur, int ik, int *ed, double *w2, int nrc,
       double *w3, int nlr, double *ex[], double *ey, int *jstat)
#else
void
s1927 (w1, nur, ik, ed, w2, nrc, w3, nlr, ex, ey, jstat)
     double *w1;
     int nur;
     int ik;
     int *ed;
     double *w2;
     int nrc;
     double *w3;
     int nlr;
     double *ex[];
     double *ey;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : 	To solve the equation  W*ex = ey  where W is an almost
*		banded matrix whose LU-factorization is contained in
*		w1, w2, w3.
*
*
* INPUT      : 	w1		- Upper part of LU-factorization of W.
*				  Dimension (1:nur*ik-1).
*		nur		- Number of upper rows of W.
*		ik		- Number of non-zero entries in each upper
*				  row of W.
*		ed		- Pointers to the diagonal elements in w1.
*		w2		- Right part of LU-factorization of W.
*				  Dimension (1:nur*nrc-1).
*		nrc		- Number of right columns of W.
*		w3		- Lower part of LU-factorization of W.
*				  Dimension (1:nlr*(nur+nlr)-1).
*		nlr		- Number of lower rows.
*		ey		- Right side of equation.
*				  Dimension (1:nur+nlr-1).
*		dim		- Number of rows of the arrays ey and ex.
*
* OUTPUT     :	ex		- Solution to W*ex=ey. Dimension (1:nur+nlr-1).
*               jstat           - Output status:
*                                 < 0: Error.
*                                 = 0: Ok.
*                                 > 0: Warning.
*
* METHOD     :	A description of the form of the matrix W may be found in the
*		subroutine s1898.
*		Since the LU-factorization of W is contained in w1 and w3, the
*		solution of W*ex=ey proceeds as follows:
*		First solve L*z=ey, then solve U*ex=z.
*
* REFERENCES :	Fortran version:
*		E.Aarn[s, CP, 79-01
*
* CALLS      :
*
* WRITTEN BY : 	Christophe R. Birkeland, SI, 1991-06
* REWRITTEN BY :
* REVISED BY :
*
*********************************************************************
*/
{
  int kpos = 0;
  int ii, jj;			/* Loop control parameters 		*/
  int di;			/* Pointer to diagonal element of W 	*/
  int midi;			/* Parameter always equal: ii-di	*/
  int dim;			/* di minus 2:  di-2			*/
  int mur;			/* Used in calculation of index for w3  */
  int nn;			/* Number of rows/columns in w3		*/
  int nlc;			/* Number of left columns in W		*/
  double wii;			/* Used to store values from matrix W	*/
  double sum;			/* Stores values for calculation of ex 	*/

  *jstat = 0;


  /* Test if legal dimension of interpolatoin problem */

  if (nur < 1 || ik < 1 || nrc < 0 || nlr < 0)
    goto err160;
  nn = nur + nlr;
  nlc = nn - nrc;
  if (ik > nlc)
    goto err160;


  /* Allocate output array ex */

  *ex = new0array (nn, DOUBLE);
  if (*ex == SISL_NULL)
    goto err101;


  /* Solve L*z = ey */

  for (ii = 0; ii < nur; ii++)
    {
      di = ed[ii];
      wii = w1[(di - 1) * nur + ii];


      /* Test for errors */

      if (ii >= nlc)
	goto err163;
      if (di < 1 || ik < di || wii == (double) 0.0)
	goto err162;
      sum = ey[ii];
      if (di > 1)
	{
	  dim = di - 1;
	  midi = ii - di + 1;
	  for (jj = 0; jj < dim; jj++)
	    sum -= w1[jj * nur + ii] * ((*ex)[jj + midi]);
	}
      (*ex)[ii] = sum / wii;
    }

  /* Solve filled part of L*z = ey */

  for (; ii < nn; ii++)
    {
      mur = ii - nur;
      wii = w3[ii * nlr + mur];
      if (wii == (double) 0.0)
	goto err162;
      sum = ey[ii];
      if (ii >= 1)
	{
	  for (jj = 0; jj < ii; jj++)
	    sum -= w3[jj * nlr + mur] * ((*ex)[jj]);
	}
      (*ex)[ii] = sum / wii;
    }

  /* Solve U*ex = z   ; Jump if filled part of U is exhausted */

  for (ii = nn - 2; ii >= nur; ii--)
    {
      sum = (*ex)[ii];
      mur = ii - nur;
      for (jj = ii + 1; jj < nn; jj++)
	sum -= w3[jj * nlr + mur] * ((*ex)[jj]);
      (*ex)[ii] = sum;
    }

  /* Test if w2 contains diagonal elements */

  if (ii >= nlc)
    goto err163;
  if (nlc < nn)
    {
      for (; ii >= 0; ii--)
	{
	  sum = (*ex)[ii];
	  for (jj = nlc; jj < nn; jj++)
	    sum -= w2[(jj - nlc) * nur + ii] * ((*ex)[jj]);
	  (*ex)[ii] = sum;
	}
    }
  for (ii = nur - 1; ii >= 0; ii--)
    {
      di = ed[ii];
      if (di < ik)
	{
	  sum = (*ex)[ii];
	  midi = ii - di + 1;
	  for (jj = di; jj < ik; jj++)
	    sum -= w1[jj * nur + ii] * ((*ex)[jj + midi]);
	  (*ex)[ii] = sum;
	}
    }

  goto out;


  /* Memory error, array ex not allocated */

err101:
  *jstat = -101;
  s6err ("s1927", *jstat, kpos);
  goto out;

  /* error in dimension of interpolation problem */

err160:
  *jstat = -160;
  s6err ("s1927", *jstat, kpos);
  goto out;

  /* W is non-invertible */

err162:
  *jstat = -162;
  s6err ("s1927", *jstat, kpos);
  goto out;

  /* w2 contains diagonal element */

err163:
  *jstat = -163;
  s6err ("s1927", *jstat, kpos);
  goto out;

out:
  return;
}
