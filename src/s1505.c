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
 * $Id:
 *
 */
#define S1505

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1505(const SISLSurf *ps1,int ider,int m1, int m2,double *ebder1,double
      *ebder2, int *ileft1,int *ileft2,double eder[],double norm[],int *jstat)
#else
void s1505(ps1,ider,m1,m2,ebder1,ebder2,ileft1,ileft2,eder,norm,jstat)
     SISLSurf   *ps1;
     int    ider;
     int m1;
     int m2;
     double *ebder1;
     double *ebder2;
     int    *ileft1;
     int    *ileft2;
     double eder[];
     double norm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the surface pointed at by ps1 over an m1 * m2
*              grid, assuming that the B-splines have been
*              pre-evaluated, by s1504, over that grid.
*              The knots et1 and et2 and grid points (x[i],y[j]) are not needed.
*              Compute ider derivatives.
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*              m1     - Number of grid points in first direction.
*              m2     - Number of grid points in first direction.
*              ileft1 - Array of indexes to the intervals in the knotvector
*                       in the first parameter direction where each subsequence
*                       of k1 non-zero B-splines are located.
*              ileft2 - Array of indexes to the intervals in the knotvector
*                       in the second parameter direction where each subsequence
*                       of k2 non-zero B-splines are located.
*              ebder1 - Triple array of dimension [(ider+1)*k1*m1] containing
*                       values of the k1 nonzero B-splines and their
*                       derivatives at the points x[0],...,x[m1-1]
*                       (i.e. pre-evaluated B-splines).
*                       These numbers are stored in the following order:
*                       First the (ider+1) derivatives of the first nonzero
*                       B-spline at x[0]. Then the (ider+1) derivatives of
*                       the second nonzero B-spline at x[0], etc.
*                       Later we repeat for x[1],... etc.
*
*              ebder2 - Triple array of dimension [(ider+1)*k2*m2] containing
*                       values of the k2 nonzero B-splines and their
*                       derivatives at the points y[0],...,y[m2-1]
*
* OUTPUT     : eder   - Array where the derivatives of the surface
*                       are placed, dimension
*                         idim * ((ider+1)(ider+2) / 2) * m1 * m2.
*                       The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2)
*                       derivative, etc. at point (x[0],y[0]),
*                       followed by the same information at (x[1],y[0]),
*                       etc.
*              norm   - Normals of surface. Is calculated if ider >= 1.
*                       Dimension is idim*m1*m2.
*                       The normals are not normalized.
*              jstat  - status messages
*                          = 2      : Surface is degenerate
*                                     at some point, normal
*                                     has zero length.
*                          = 1      : Surface is close to
*                                     degenerate at some point
*                                     Angle between tangents,
*                                     less than angular tolerance.
*                          = 0      : ok
*                          < 0      : error
*
* METHOD     : The code and method is similar to that of s1421 except that
*              the B-splines and their derivatives have already been
*              calculated (and are given in ebder1 and ebder2) and
*              we evaluate over a whole grid rather than at one point.
*              The method is to find the control points and control derivatives
*              of each isocurve in x (fixed y value). We then
*              evaluate at each x value along the isocurve.
*
* CALLS      : s6err     - Error handling routine
*              s6strider - Make derivative of rational expression
*
* WRITTEN BY : Michael Floater, SINTEF, May 1998.
*********************************************************************
*/
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of error.                          */
  int i1,i2;          /* Loop variables. */
  int i1pos,i2pos;    /* Offset indexes. */
  int kn1,kn2;        /* The number of B-splines accociated with the knot
			 vectors st1 and st2.                            */
  int kk1,kk2;        /* The polynomial order of the surface in the two
			 directions.                                     */
  int kdim;           /* The dimension of the space in which the surface
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int kleft1,kleft2;  /* Local versions of knot intervals.            */
  int ki,kx,kjh;      /* Control variables in for loops and for stepping
			 through arrays.                                 */
  int kih2;           /* Index for stepping through ebder2. */
  int ky,kl,kl1,kl2;  /* Control variables in for loops and for stepping
			 through arrays.                                 */
  double *scoef;      /* The B-spline coefficients of the surface.
			 This is an array of dimension [kn2*kn1*kdim].   */
  double tt;          /* Dummy variable used for holding an array element
			 in a for loop.                                  */
  double *ew=SISL_NULL;    /* Pointer to an array of dimension [kn1*(ider+1)*kdim]
			 which will be used to store the result of the first
			 matrix multiplication. */
  double *sder=SISL_NULL;  /* Pointer to array used for storage of points, if
			 rational has room for homogenous coordinates. */
  double *enorm=SISL_NULL; /* Array for surface normal. */
  int size;           /* Space occupied by points and derivs at one eval. */
  int sizeh;          /* Space occupied by homogeneous points and derivs . */
  int size1,size2;    /* Useful variables. */
  int ederpos;       /* Index of position in eder. */
  int normpos;       /* Index of position in norm. */

  int knumb2;         /* Necessary size of ew */

  int tot,temp;       /* Temporary variables. */

  kn1 = ps1 -> in1;
  kn2 = ps1 -> in2;
  kk1 = ps1 -> ik1;
  kk2 = ps1 -> ik2;
  kdim = ps1 -> idim;
  /* Check the input. */

  if (kdim < 1) goto err102;
  if (kk1 < 1) goto err115;
  if (kn1 < kk1 || kn2 < kk2) goto err116;
  if (ider < 0) goto err178;

  *jstat = 0;

  if (ps1->ikind == 2 || ps1->ikind == 4)
  {
    scoef = ps1 -> rcoef;
    kdim +=1;
  }
  else
  {
    scoef = ps1 -> ecoef;
  }
  sizeh = kdim*(ider+1)*(ider+2)/2;
  if((sder=newarray(sizeh,DOUBLE)) == SISL_NULL) goto err101;
  if((enorm=newarray(ps1->idim,DOUBLE)) == SISL_NULL) goto err101;

  size = ps1->idim*(ider+1)*(ider+2)/2;
  size1 = (ider+1)*kk1;
  size2 = (ider+1)*kk2;

  /* Allocate space for B-spline values and derivatives and one work array. */

  knumb2 = kn1*(ider+1)*kdim;
  if((ew=newarray(knumb2,DOUBLE)) == SISL_NULL) goto err101;

  ederpos = 0;
  normpos = 0;

  /* Run through grid points in the y direction. */
  for(i2=0, i2pos=0; i2<m2; i2++, i2pos += size2)
  {
    kleft2 = ileft2[i2];

    /* Compute the control points (and control derivatives
       of the v = x[i2] isocurve. */

    /* Set all the elements of ew to 0. */
    for(ki=0; ki<knumb2; ki++) ew[ki] = DZERO;

    /* ki steps through the appropriate kk2 rows of B-spline coefficients
       while kih2 steps through the B-spline value and derivatives for the
       B-spline given by ki.            */

    kih2 = i2pos;
    for (ki=kleft2-kk2+1; ki<=kleft2; ki++)
    {
      /* kx counts through the ider+1 derivatives to be computed.
         kjh steps through ew once for each ki to accumulate the
         contribution from the different B-splines.
         kl1 points to the first component of the first B-spline coefficient
         in row no. ki of the B-spline coefficient matrix.
      */

      kjh = 0; kl1 = kdim*ki*kn1;
      for (kx=0; kx<=ider; kx++)
      {
        /* The value of the B-spline derivative is stored in tt while
           kl2 steps through the kdim components of all the B-spline
           coefficients that multiply nonzero B-splines along st1.
        */

        tt = ebder2[kih2++]; kl2 = kl1;
        for (kl=0; kl<kdim*kn1; kl++,kjh++,kl2++)
        {
           ew[kjh] += scoef[kl2]*tt;
        }
      }
    }

    /* Run through grid points in the x direction
       evaluating along the iso-curve defined by y = y[i2]. */
    for(i1=0, i1pos=0; i1<m1; i1++, i1pos += size1)
    {
      kleft1 = ileft1[i1];

      /* Set all the elements of sder to 0. */
      for(ki=0; ki<sizeh; ki++) sder[ki] = DZERO;

      for(ky=0; ky<=ider; ky++)
      {
	kl1 = kdim * (ky * kn1 + kleft1 - kk1 + 1);
	for(kx=0; kx<=ider-ky; kx++)
	{
          tot = kx + ky;
          temp = ((tot * (tot+1)) >> 1) + ky;
	  kjh = temp * kdim;

	  for(ki=0; ki<kk1; ki++)
          {
	    tt = ebder1[i1pos + (ider+1) * ki + kx];
	    for(kl=0; kl<kdim; kl++)
	    {
	      sder[kjh+kl] += ew[kl1 + kdim * ki + kl] * tt;
	    }
          }
        }
      }

      /* If rational surface calculate the derivatives based on
         derivatives in homogenous coordinates */

      if (ps1->ikind == 2 || ps1->ikind == 4)
      {
        s6strider(sder,ps1->idim,ider,eder+ederpos,&kstat);
        if (kstat<0) goto error;
      }
      else
      {
        for(ki=0; ki<sizeh; ki++) eder[ederpos+ki] = sder[ki];
      }

      /* Calculate normal if idim==3 and ider>0. */

      if (ider>0 && ps1->idim ==3)
        {

          s6crss(eder+ederpos+ps1->idim,eder+ederpos+2*ps1->idim,enorm);

          s6norm(enorm,3,norm+normpos,&kstat);
         }

      ederpos += size;
      normpos += 3;
    }
  }

  goto out;

  /* Not enough memory. */
 err101: *jstat = -101;
  s6err("s1505",*jstat,kpos);
  goto out;

  /* kdim less than 1. */
 err102: *jstat = -102;
  s6err("s1505",*jstat,kpos);
  goto out;

  /* Polynomial order less than 1. */
 err115: *jstat = -115;
  s6err("s1505",*jstat,kpos);
  goto out;

  /* Fewer B-splines than the order. */
 err116: *jstat = -116;
  s6err("s1505",*jstat,kpos);
  goto out;

  /* Illegal derivative requested. */
 err178: *jstat = -178;
  s6err("s1221",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("s1505",*jstat,kpos);
  goto out;

 out:
  /* Free memory. */
  if(sder != SISL_NULL) freearray(sder);
  if(enorm != SISL_NULL) freearray(enorm);
  if (ew != SISL_NULL) freearray(ew);

    return;
}
