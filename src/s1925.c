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
 * $Id: s1925.c,v 1.3 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1925

#include "sislP.h"
#define s1925_MAX_ARRAY_SIZE 50

#if defined(SISLNEEDPROTOTYPES)
void
s1925 (double etau[], double epoint[], int inbpnt, int eder[],
       double et[], double ebcoef[], int in, int ik, int iright, int dim,
       double ew1[], int nur, int ed[], double ew2[], int inrc, double ew3[],
       int inlr, int *jstat)
#else
void
s1925 (etau, epoint, inbpnt, eder, et, ebcoef, in, ik, iright, dim, ew1, nur,
       ed, ew2, inrc, ew3, inlr, jstat)
     double etau[];
     double epoint[];
     int inbpnt;
     int eder[];
     double et[];
     double ebcoef[];
     int in;
     int ik;
     int iright;
     int dim;
     double ew1[];
     int nur;
     int ed[];
     double ew2[];
     int inrc;
     double ew3[];
     int inlr;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : 	To calculate B-spline coefficients ebcoef required for
*		interpolation of epoint(i) at etau(i); i=0,1,...,inbpnt-1
*
* INPUT      :	etau	- Vector containing the parametrization of the data
*			  points/derivatives/second-derivatives.
*		epoint	- Points/derivatives/second-derivatives to be inter-
*			  polated. Only entries from epoint are used.
*		inbpnt	- Number of points/derivatives/second-derivatives
*			  to be interpolated.
*		eder	- Indicator for the order of derivative of the data
*			  to be interpolated.
*		et	- Knot vector of B-spline. Length: in+ik
*		in	- Number of vertices in B-spline to be produced.
*		ik	- Order of B-spline to be produced
*		iright	- Number of right hand sides of dimension (dim*in)
*		dim	- The dimension of space in which the points lie
*		ew1,    - Arrays defining elements of W.
*		nur	- Number of upper rows of W.
*		ed      -
*		ew2     - Detailed description: see subroutine s1898.
*		inrc	- Number of right columns of W.
*		ew3     -
*		inlr	- Number of lower rows of W
*
* OUTPUT     :	ebcoef	- Coefficients of interpolating B-spline.
*			  Stored in array of dimension (1:dim*inbpnt).
*               ew1     - Arrays defining elements of LU-factorization of W.
*		ed      -
*               ew2     -
*		ew3 	- Detailed description: see subroutine s1898
*               jstat   - Output status:
*                         < 0: Error.
*                         = 0: Ok.
*                         > 0: Warning.
*
* METHOD     : 	The i-th equation of the set  W*ebcoef = B  enforces
*		interpolation of B(i) at etau(i). Thus B(i)=epoint(i).
*		W is generated row by row and stored in ew1,ew2 and ew3
*		as described in subroutine s1898. The system is solved by
*		one call to s1898 and then by 'dim' calls to s1899.
*
* REFERENCES : 	Fortran version:
*		E.Aarn[s, CP, 1979-01
*
* CALLS      : 	s1897,s1926,s1927,s6err.
*
* WRITTEN BY : 	Christophe R. Birkeland
* REVISED BY :  Johannes Kaasa, SI, June 1992 (Introduced iadd_save)
*
*********************************************************************
*/
{
  int kstat = 0;
  int kpos = 0;			/* Position of error			*/
  int open;			/* Used as a boolean parameter to
				 * indicate open or closed curve:
				 * open=TRUE   ; Open curve
			  	 * open=FALSE  ; Closed curve  		*/
  int left;			/* An integer chosen (usually) so that
				 * et(left-1)<=point<et(left)		*/
  int leftmax;
  int leftmin;
  int leftdel;
  int kmod;			/* Used in calculation of ew3 index  	*/
  int isum;
  int ii, jj, kl, stop;		/* Loop control parameters		*/
  int dim1;			/* Loop control parameter,
				 * values: 0..dim-1 			*/
  int imnur;			/* Equals ii-nur			*/
  int nn;			/* Equals nur+inlr
				 * i.e. Number of rows/columns in W	*/
  int nlc;			/* Number of lower columns
				 * Equals inbpnt-inrc			*/
  int kk;			/* Minimum: ik or nlc			*/
  double tk;			/* ik-th element of knot vector 	*/
  double taudel;		/* to left and taui			*/
  int lfmkm;			/* Equals left-ik			*/
  int iadd;
  int iadd_save;
  int kmiadd;
  int isub;
  int kmisub;			/* Equals ik-isub			*/
  int ish;
  int ideri;			/* Derivative order indicator 		*/
  int store;

  double taui;			/* Parametrization value = etau[ii] 	*/
  double *mcoef = SISL_NULL;		/* Arrays for internal use in 		*/
  double *ebder = SISL_NULL;		/* this subroutine			*/
  double sarray[s1925_MAX_ARRAY_SIZE];
  int alloc_needed=FALSE;
  
  *jstat = 0;

  nn = nur + inlr;

  /* Test if legal input */

  if (ik < 1)
    goto err109;

  if ((nur < 0) || ((nur + inlr) != inbpnt) || (inrc < 0) || (inlr < 0))
    goto err160;

  nlc = inbpnt - inrc;
  kk = MIN (ik, nlc);
  tk = et[ik - 1];


  /* Test if open or closed curve */

  if (inbpnt == in)
    open = TRUE;
  else
    open = FALSE;

  if (open == TRUE)
    {
      /* Open curve */

      taudel = (double) 0.0;
      leftdel = 0;
      leftmin = ik - 1;
      leftmax = in -1;
      if ((inrc != 0) || (inlr != 0))
	goto err160;
    }
  else
    {
      /* Closed curve, note: we assume that et(in) < etau(inbpnt) */

      if ((inbpnt + ik - 1) != in)
	goto err160;
      leftdel = in +1 - ik;
      iadd = (ik - 1) / 2;
      leftmax = leftdel + iadd - 1;
      leftmin = ik - iadd - 1;
      taudel = et[in] -tk;

      if ((iadd != inrc) || ((inrc + inlr) != (ik - 1)))
	goto err160;

      if (inrc > 0)
	{
	  stop = nur * inrc;
	  for (ii = 0; ii < stop; ii++)
	    ew2[ii] = (double) 0.0;
	}
    }

  /* Band part of W */

  /* Allocate array ebder */

  left = leftmin;
  if (ik > s1925_MAX_ARRAY_SIZE)
    {
       if ((ebder = newarray (ik, DOUBLE)) == SISL_NULL)
	 goto err101;
	alloc_needed = TRUE;
    }
  else
    ebder = sarray;
  
  for (ii = 0; ii < nur; ii++)
    {
      taui = etau[ii];
      ideri = eder[ii];


      /* Locate left so that  et[left] <= taui < et[left+1] */

      while (left < leftmax && et[left + 1] <= taui)
	left++;


      /* et(left-1) <= taui < et(left)  */

      ed[ii] = ii - (left - ik);


      /* Test if error in interpolation problem */

      if ((ed[ii] < 1) || (ed[ii] > ik))
	goto err165;


      iadd = MAX (0, ik - left - 1);
      iadd_save = iadd;
      kmiadd = ik - iadd;


      /* Compute the value and ideri first derivative of the
      ik (possibly) nonzero B-spline associated with the knot
      vector et at a point 				    */

      if (iadd > 0)
	{
	  s1897 (et, ik, taui + taudel, left + leftdel, ideri, ebder, &kstat);
	  if (kstat < 0)
	    goto error;

	  ed[ii] -= iadd;
	  ish = inrc - iadd;
	  if ((ish < 0) || (kmiadd < 0))
	    goto err160;

	  for (jj = 0; jj < iadd; jj++)
	    {
	      ew1[(jj + kmiadd) * nur + ii] = (double) 0.0;
	      ew2[(jj + ish) * nur + ii] = ebder[jj];
	    }
	}
      else
	{
	  s1897 (et, ik, taui, left, ideri, ebder, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      isub = MAX (0, ii - ed[ii] + kmiadd - nlc + 1);
      kmisub = ik - isub;
      if (isub > 0)
	{
	  ish = isub - (kmiadd - kk);
	  ed[ii] += ish;
	  iadd -= ish;
	  if (kmisub < 0)
	    goto err160;
	  stop = MIN (ik, kmisub + inrc);

	  for (jj = kmisub; jj < stop; jj++)
	    {
	      ew1[(jj - kmisub) * nur + ii] = (double) 0.0;
	      ew2[(jj - kmisub) * nur + ii] = ebder[jj];
	    }
	}
      for (jj = iadd_save; jj < kmisub; jj++)
	ew1[(jj - iadd) * nur + ii] = ebder[jj];


      /* Test if error in dimension of interpolation problem */

      if ((ideri < 0) || (ik <= ideri))
	goto err160;
    }

  /* Band part of W is now completed */

  if (ii < inbpnt)
    {
      /* Will compute lower, filled part of W for closed
	 curve interpolation */

      if ((inbpnt + ik - 1) != in)
	goto err160;
      store = inlr * inbpnt;
      for (jj = 0; jj < store; jj++)
	ew3[jj] = (double) 0.0;


      /* Repeat until filled part of W is completed */

      for (; ii < inbpnt; ii++)
	{
	  taui = etau[ii];


	  /* Locate left so that  et[left] <= taui < et[left+1] */

	  while (left < in -1 && et[left + 1] <= taui)
	    left++;


	  /* et(left-1) <= taui < et(left)  */

	  ideri = eder[ii];


	  /* Compute the value and the ideri first derivatives of the
	     ik (possibly) nonzero B-spline associated with the knot
	     vector et at the point (taui) */

	  s1897 (et, ik, taui, left, ideri, ebder, &kstat);
	  if (kstat < 0)
	    goto error;

	  imnur = ii - nur;
	  lfmkm = left - ik;
	  for (jj = 0; jj < ik; jj++)
	    {
	      isum = jj + lfmkm + 1;
	      if (isum >= 0)
		kmod = isum % inbpnt;
	      if (isum < 0)
		kmod = (isum + 1) % inbpnt + in -1;
	      ew3[kmod * inlr + imnur] = ebder[jj];
	    }
	  if ((ideri < 0) || (ik <= ideri))
	    goto err160;
	}
    }
  if (inlr != (inbpnt - nur))
    goto err160;


  /* W is now contained in ew1, ew2 and ew3 as required
     by the subroutine s1898  */

  s1926 (ew1, nur, kk, ed, ew2, inrc, ew3, inlr, &kstat);
  if (kstat < 0)
    goto error;

  store = iright * dim * inbpnt;
  for (jj = 0; jj < store; jj++)
    ebcoef[jj] = epoint[jj];


  /* epoint is now properly contained in ebcoef.
   * Solve interpolation equations 		 */

  if (nn > s1925_MAX_ARRAY_SIZE)
    {
       if (alloc_needed)
	 {
	    if ((ebder = increasearray(ebder,nn,DOUBLE)) == SISL_NULL)
	      goto err101;
	 }
       else
	 {
	    if ((ebder = newarray(nn,DOUBLE)) == SISL_NULL)
	      goto err101;
	    alloc_needed = TRUE;
	 }
    }

  for (kl = 0; kl < iright; kl++)
    for (dim1 = 0; dim1 < dim; dim1++)
      {
	store = inbpnt * dim * kl + dim1;
	for (jj = 0; jj < nn; jj++, store += dim)
	  ebder[jj] = ebcoef[store];

	s1927 (ew1, nur, kk, ed, ew2, inrc, ew3, inlr, &mcoef, ebder, &kstat);
	if (kstat < 0)
	  goto error;

	store = inbpnt * dim * kl + dim1;
	for (jj = 0; jj < nn; jj++, store += dim)
	  ebcoef[store] = mcoef[jj];

        if(mcoef != SISL_NULL)       /* KYS 200594: healed memory leak */
        {
          freearray(mcoef);
          mcoef = SISL_NULL;
        }
      }

  goto out;


  /* Error in array allocations */

err101:
  *jstat = -101;
  s6err ("s1925", *jstat, kpos);
  goto out;

  /* Order of B-spline zero or negative */

err109:
  *jstat = -109;
  s6err ("s1925", *jstat, kpos);
  goto out;

  /* Error in dimension of interpolation problem */

err160:
  *jstat = -160;
  s6err ("s1925", *jstat, kpos);
  goto out;

  /* Error in lower level routine */

error:
  *jstat = kstat;
  s6err ("s1925", *jstat, kpos);
  goto out;

  /* Error in interpolation problem */

err165:
  *jstat = -165;
  s6err ("s1925", *jstat, kpos);
  goto out;

out:
  if (alloc_needed)
    freearray (ebder);
  if (mcoef != SISL_NULL)
    freearray (mcoef);
  return;
}
#undef s1925_MAX_ARRAY_SIZE
