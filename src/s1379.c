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
 * $Id: s1379.c,v 1.3 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1379

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
  s1379(double ep[],double ev[],double epar[],int im,int idim,
	SISLCurve **rcurve,int *jstat)
#else
void s1379(ep,ev,epar,im,idim,rcurve,jstat)
     double ep[];
     double ev[];
     double epar[];
     int    im;
     int    idim;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose: To compute the cubic Hermit interpolant to the data given
*          by the points ep, the derivatives ev and paramterization epar.
*          The curve is represented as a B-spline curve.If the first and
*          last points are exactly equal (down to the last bit) the a periodic
*          basis with the first 1 single knot, then a 3 tupple knot is made
*          at the start and 3 tupple knot followed by a single knot at the
*          end. If also the derivatives at the start and end are equal then
*          all knots will be double.
*
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..)
*          ev     - Array containing the derivatives in sequence
*                   (x,y,..,x,y,..)
*          epar   - Parametrization array. The array should be increasing
*                   in value.
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
* Output:
*          rcurve - Pointer to the curve produced
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*     The knot vector will have 4-tupple, 3-tupple or 2-tupple  knots at
*     epar[0] and epar[im-1]. This is decided by checking if the input
*     data is cyclic:
*        if first point != last point  4-tupple knots
*        if first == last and there derivatives different 3-tupple
*        if both position and derivatives equal 2-tupple
*     In the case not 4-tupple knots the new knot intervals are made cyclic.
*
*     Suppose we have reached data point no. j which corresponds to the
*     parameter value z=epar(j), i.e., the knot vector et has a double
*     knot at z, and these knots must be et(2*j+1) and et(2*j+2).
*     then there are only two B-splines which are nonzero at z, namely
*     B(2*j-1) and B(2*j) (remember everything is cubic), and these are
*     also the only B-splines with nonzero derivative at z. we can
*     therefore determine the coefficients of these two B-splines,
*     ec(2*j-1) and ec(2*j), by solving a 2x2 linear system OF
*     equations, and this can be done directly.
*     Suppose the distance from z to the previous knot is h1 (measured
*     in parameter interval, so h1= et(2*j+2)-et(2*j), and the distance
*     to the next knot is h2 = et(2*j+3) -et(2*j+1). then the two
*     B-splines and their derivatives have the following values at z:
*
*                B(2*j-1,z)=h2/(h1+h2),    B(2*j,z)=h1/(h1+h2)
*               dB(2*j-1,z)=-3/(h1+h2),   dB(2*j,z)=3/(h1+h2).
*
*     solving the linear system for ec(2*j-1) and ec(2*j) yields
*
*                    ec(2*j-1)=ep(j) - h1*ev(j)/3,
*                      ec(2*j)=ep(j) + h2*ev(j)/3.
*
*     it can also be easily checked that these equations are valid for
*     the first two and last two coefficients as well, provided one
*     sets h1=0 and h2=0, respectively.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI 1988-11
* Revised by : Tor Dokken, SI 1992-04
* Revised by : Johs. Kaasa,SI 1993-01
*
*********************************************************************
*/
{
  int ki,kj;          /* Loop variables                              */
  int kk;             /* Polynomial order                            */
  int kn;             /* Number of vertices                          */
  int kpoint;         /* Pointer into point and derivative array     */
  int kcoef;          /* Pointer into coefficient array              */
  int kpos=0;         /* Position of error                           */
  int kthis;          /* Current point                               */
  int kstat=0;        /* Status variable                             */
  int kcycpos = 1;    /* Flag telling if first and last points are equal */
  int kcycder = 1;    /* Flag telling if first and last derviatives are equal */
  double *st=SISL_NULL;    /* Knot vector                                 */
  double *scoef=SISL_NULL; /* B-spline vertices                           */
  double th1,th2;     /* Parameter intervals                         */



  /* Check input */

  if (im < 2)   goto err181;
  if (idim < 1) goto err102;

  /* Set the dimension and order of the spline space */

  kn = 2*im;
  kk = 4;

  /* Allocate arrays for temporary storage of knots and vertices */

  st    = newarray(kn+kk,DOUBLE);
  if (st == SISL_NULL) goto err101;
  scoef = newarray(idim*kn,DOUBLE);
  if (scoef == SISL_NULL) goto err101;

  /* Check if the curve is periodic, e.g. if first and last points are
     equal and/or that first and last derivates are equal */

  /*  for (kj=0, kcycpos=1 ; kj<idim && kcycpos == 1 ; kj++)
     if (ep[kj] != ep[idim*(im-1)+kj]) kcycpos =0; */
  for (kj=0, kcycpos=1 ; kj<idim && kcycpos == 1 ; kj++)
     if (DNEQUAL(ep[kj], ep[idim*(im-1)+kj])) kcycpos =0;

  /*  for (kj=0, kcycder=1 ; kj<idim && kcycder == 1 ; kj++)
    if (ev[kj] != ev[idim*(im-1)+kj]) kcycder= 0; */
  for (kj=0, kcycder=1 ; kj<idim && kcycder == 1 ; kj++)
    if (DNEQUAL(ev[kj], ev[idim*(im-1)+kj])) kcycder= 0;

  /* Make the knot vector, first all knots except the two first and the two last */

  for (ki=2,kj=0 ; ki<kn+2 ; ki+=2, kj++)
    st[ki] = st[ki+1] = epar[kj];



  /* Make the two first and two last knots */

  if (kcycder == 1 && kcycpos == 1)
    {
      /* Two first knots to be shifted */

      st[0]= st[1] = epar[0] - (epar[im-1]-epar[im-2]);
      st[kn+2]= st[kn+3] = epar[im-1] + epar[1] - epar[0];
    }
  else if (kcycder ==0 && kcycpos ==1)
    {
      /* First and last knot to be shifted */

      st[0] = epar[0] - (epar[im-1]-epar[im-2]);
      st[1] = st[2];
      st[kn+2] = st[kn];
      st[kn+3] = epar[im-1] + epar[1] - epar[0];
    }
  else
    {
      /* k-regular basis */

      st[0] = st[1] = st[2];
      st[kn+2] = st[kn+3] = st[kn];
    }

  /* Compute knot vector and coefficients as indicated above */

  for (kj=0, kcoef=0, kpoint = 0 ; kj < kn ; kj+=2)
    {
      th1 = st[kj+3] - st[kj+1];
      th2 = st[kj+4] - st[kj+2];

      /*  Compute coefficient no kj */

      kthis = kpoint;
      for (ki=0;ki<idim;ki++,kpoint++)
        {
	  scoef[kcoef++] = ep[kpoint] - th1*ONE_THIRD*ev[kpoint];
        }

      /*  Compute coefficient no kj+1 */

      kpoint = kthis;
      for (ki=0;ki<idim;ki++,kpoint++)
        {
	  scoef[kcoef++] = ep[kpoint] + th2*ONE_THIRD*ev[kpoint];
        }
    }

  /* Make new curve object */

  *rcurve = newCurve(kn,kk,st,scoef,1,idim,1);
  if (*rcurve == SISL_NULL) goto err101;

  /* Remove unneccesarry knots */

  s6crvcheck(*rcurve,&kstat);
  if (kstat<0) goto error;

  /* Periodicity flag */
  if (kcycpos)
    {
       test_cyclic_knots((*rcurve)->et,(*rcurve)->in,(*rcurve)->ik,&kstat);
       if (kstat<0) goto error;
       if (kstat == 2) (*rcurve)->cuopen = SISL_CRV_PERIODIC;
    }

  /* Calculation completed */

  *jstat = 0;
  goto out;


  /* Error in space allocation. Return zero. */


  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1379",*jstat,kpos);
  goto out;


  /* Dimension less than 1*/
 err102: *jstat = -102;
  s6err("s1379",*jstat,kpos);
  goto out;

  /* Too few interpolation conditions */

 err181: *jstat = -181;
  s6err("s1379",*jstat,kpos);
  goto out;

 error:  *jstat =kstat;
  s6err("s1379",*jstat,kpos);
  goto out;

 out:
  if (st != SISL_NULL) freearray(st);
  if (scoef != SISL_NULL) freearray(scoef);

  return;
}
