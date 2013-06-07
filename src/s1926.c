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
 * $Id: s1926.c,v 1.2 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1926

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1926 (double *w1, int nur, int ik, int *ed, double *w2, int nrc,
       double *w3, int nlr, int *jstat)
#else
void
s1926 (w1, nur, ik, ed, w2, nrc, w3, nlr, jstat)
     double *w1;
     int nur;
     int ik;
     int *ed;
     double *w2;
     int nrc;
     double *w3;
     int nlr;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :  To calculate the LU-factorization of an almost banded matrix W.
*
*
* INPUT      :  w1  	- Upper part of W. Dimension (1:nur*ik).
*		nur	- Number of rows of upper part.
*		ik	- Number of non-zero entries of each upper row.
*		ed	- Pointers do diagonal elements of W.
*		w2	- Right part of W. Dimension (1:nur*nrc).
*		nrc	- Number of right columns of W.
*		w3	- Lower part of W. Dimension (1:nlr*(nur+nlr)).
*		nlr	- Number of lower rows.
*
* OUTPUT     :  w1      -
*               w2      -
*               w3      - Contain the LU-factorization on return. L is lower
*			  triangular, and U is upper triangular with ones on
*		 	  the diagonal.
*               jstat   - Output status:
*                          < 0: Error.
*                          = 0: Ok.
*                          > 0: Warning.
*
* METHOD     : The matrix W is of the following form:
*
*                 1 2 3 4 5 6 7 8 9   ED(.)
*              1  D X X X         Y     1
*              2  X D X X         Y     2
*              3    X D X X       Y     2
*              4        D X X X   Y     1
* 	       5          D X X X Y     1
*	       6          X D X X Y     2
*   	       7          X X D X Y     3
*	       8  Z Z Z Z Z Z Z D Z     -
*	       9  Z Z Z Z Z Z Z Z D     -
*
*
*               w1[1:nur*ik] is used to store the upper nur rows of W,
*		each row containing ik non-zero entries. The positions of
*		W[i,i] ; i=1,2,...,[nur+nlr]  are indicated in ed[1:nur]. In
*               the above example: nur=7, ik=4 and ed[1:7] is indicated to the
*		right. The diagonal elements of W are indicated with D's,
*		and the elements of w1 are indicated with X's.
*
*		w2[1:nur*nrc] is used to store the right nrc columns of W,
* 		each column containing nur entries. In the above example,
*		nrc=1, w2 is indicated with Y's.
*
*		w3[1:nlr*(nur+nlr)] is used to store the lower nlr rows of W.
*		In the above example, nlr=2 and nur+nlr=9 (total number of
*		rows/columns). w3 is indicated with Z's in the above example.
*
* REFERENCES :  Fortran version:
*               E.Aarn[s, CP, 1979-01
*
* CALLS      : s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland
*
*********************************************************************
*/
{
  int kpos = 0;
  int ii, jj;			/* Loop control parameters 		*/
  int ll;
  int nn;			/* Number of rows/columns in W 		*/
  int nlc;			/* Number of left columns in W 		*/
  int di;			/* Pointer to diagonal element of W   	*/
  int midi;			/* Parameters used in elimination alg.  */
  int midl;
  int korr;			/* midl - midi 				*/
  int mur;			/* Used in calculation of index for w3  */
  double wii;			/* Used to store values from matrix W   */
  double wli;

  *jstat = 0;


  /* Test if legal dimension of interpolation problem */

  if (nur < 1 || (nur >= 1 && ik < 1) || nrc < 0 || nlr < 0)
    goto err160;

  nn = nur + nlr;
  nlc = nn - nrc;
  if (ik > nlc)
    goto err160;


  /* Elimination scheme, jump if band part of W is completed */

  for (ii = 0; ii < nur; ii++)
    {
      di = ed[ii];
      wii = w1[(di - 1) * nur + ii];


      /* Test for errors */

      if (ii >= nlc)
	goto err163;
      if ((di < 1) || (ik < di) || (wii == 0.0))
	goto err162;


      /* Jump if W(ii,jj) is trivially zero, jj = ii+1,ii+2,...,nlr */

      if (di < ik)
	{
	  for (jj = di; jj < ik; jj++)
	    w1[jj * nur + ii] /= wii;


	  /* Perform elimination row by row */

	  midi = ii - di;
	  for (ll = ii + 1;; ll++)
	    {
	      /* Jump if ii-th element of rows of band-part has been
               * eliminated */

	      if (ll >= nur)
		break;
	      midl = ll - ed[ll];


	      /* Jump if W(ii,jj) is trivially zero, jj = ll,ll+1,...,nur */

	      if (midl >= ii)
		break;
	      korr = midl - midi;
	      wli = w1[(di - korr - 1) * nur + ll];
	      for (jj = di; jj < ik; jj++)
		w1[(jj - korr) * nur + ll] += -w1[jj * nur + ii] * wli;
	    }

	  /*  Eliminate ii-th column of w3 using ii-th row from w1 */

	  if (nlr > 0)
	    for (ll = 0; ll < nlr; ll++)
	      {
		wli = w3[ii * nlr + ll];
		for (jj = di; jj < ik; jj++)
		  w3[(jj + midi + 1) * nlr + ll] -= w1[jj * nur + ii] * wli;
	      }
	}
    }

  /* Apply the above elimination scheme on w2 */

  if (nrc > 0)
    {
      /* Jump if band part of W is completed or if system error
       * occures, i.e. if w2 contains some diagonal elements of W) */

      for (ii = 0; ii < nur; ii++)
	{
	  /* Test for error */

	  if (ii > nlc)
	    goto err163;

	  di = ed[ii];
	  wii = w1[(di - 1) * nur + ii];
	  for (jj = 0; jj < nrc; jj++)
	    w2[jj * nur + ii] /= wii;


	  /* Perform elimination row by row */

	  midi = ii - di;
	  for (ll = ii + 1;; ll++)
	    {

	      /*  Jump if ii-th element of rows of band-part has been
	          eliminated                                          */

	      if (ll >= nur)
		break;
	      midl = ll - ed[ll];


	      /* Jump if W(ii,jj) is trivially zero, jj = ll,ll+1,...,nur */

	      if (midl >= ii)
		break;
	      korr = midl - midi;
	      wli = w1[(di - korr - 1) * nur + ll];
	      for (jj = 0; jj < nrc; jj++)
		w2[jj * nur + ll] -= w2[jj * nur + ii] * wli;
	    }

	  /*  Eliminate ii-th column of w3 using ii-th row from w2 */

	  for (ll = 0; ll < nlr; ll++)
	    {
	      wli = w3[ii * nlr + ll];
	      for (jj = nlc; jj < nn; jj++)
		w3[jj * nlr + ll] -= w2[(jj - nlc) * nur + ii] * wli;
	    }
	}
    }

  /* Eliminate w3-part of W */

  if (ii >= nn)
    goto out;
  for (; ii < nn; ii++)
    {

      /* 1 <= ii <= nn */

      mur = ii - nur;
      wii = w3[ii * nlr + mur];
      if (wii == (double) 0.0)
	goto err162;


      /*  1 <= ii < nn */

      for (jj = ii + 1; jj < nn; jj++)
	w3[jj * nlr + mur] /= wii;
      for (ll = mur + 1; ll < nlr; ll++)
	{
	  wli = w3[ii * nlr + ll];
	  for (jj = ii + 1; jj < nn; jj++)
	    w3[jj * nlr + ll] -= w3[jj * nlr + mur] * wli;
	}
    }

  goto out;


  /* W may be non-invertible */

err162:
  *jstat = -162;
  s6err ("s1926", *jstat, kpos);
  goto out;

  /* Error in dimension in interpolation problem */

err160:
  *jstat = -160;
  s6err ("s1926", *jstat, kpos);
  goto out;

  /* w2 contains diagonal elements */

err163:
  *jstat = -163;
  s6err ("s1926", *jstat, kpos);
  goto out;

out:
  return;
}
