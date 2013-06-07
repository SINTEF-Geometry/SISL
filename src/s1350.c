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
