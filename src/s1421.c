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
 * $Id: s1421.c,v 1.5 2001-03-19 15:58:49 afr Exp $
 *
 */
#define S1421

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1421(SISLSurf *ps1,int ider,double epar[],int *ilfs,int *ilft,
	   double eder[],double enorm[],int *jstat)
#else
void s1421(ps1,ider,epar,ilfs,ilft,eder,enorm,jstat)
     SISLSurf   *ps1;
     int    ider;
     double epar[];
     int    *ilfs;
     int    *ilft;
     double eder[];
     double enorm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the surface pointed at by ps1 at the parameter
*              value epar. Compute ider derivatives.
*
*
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.    
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*              epar   - Parameter-value at which to calculate. Dimension
*                       of epar is 2.
*
*                
*
* INPUT/OUTPUT : ilfs  - Pointer to the interval in the knotvector
*                        in first parameter direction where epar[0] 
*                        is found. The relation
*                          
*                          et1[ilfs] <= epar[0] < et1[ilfs+1]
* 
*                        where et1 is the knotvektor should hold.
*                        ilfs is set equal to zero at the first call
*                        to the routine. 
*                ilft  - Corresponding to ilfs in the second parameter
*                        direction.
*             
*
*
*
* OUTPUT     : eder   - Array where the derivative of the curve in
*                       apar is placed. The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative, etc. Dimension of eder is 
*                       idim*(1+2+...+(ider+1)).
*              enorm  - Normal of surface. Is calculated if ider >= 1.
*                       Dimension is idim. The normal is not normalized.
*              jstat  - status messages  
*                                         = 2      : Surface is degenerate
*                                                    at the point, normal
*                                                    has zero length
*                                         = 1      : Surface is close to
*                                                    degenerate at the point
*                                                    Angle between tangents,
*                                                    less than angular tolerance
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : Suppose that the given surface is of the form
*
*                 s(u,v) = sum(i,j) c(i,j)*B(i,k1,t1)*B(j,k2,t2)
*
*              where c is the matrix of B-spline coefficients (each c(i,j)
*              is a vector with idim components),
*              B(i,k1,t1) the B-splines accociated with the knot vector t1,
*              and B(j,k2,t2) the B-splines accociated with the knot vector t2.
*              This may be expressed in matrix form as
*
*                           s(u,v) = tran(B2(v)) * C * B1(u),     (1)
*
*              where
*
*                tran(B1(u))=(B(1,k1,t1)(u),B(2,k1,t1)(u),...,B(n1,k1,t1(u)))
*
*              is the vector of B-spline values at u, and tran(a) denotes
*              the transpose of the vector a.
*              It is known that for a given value of u, there are at most
*              k1 (the order of the splines associated with t1) nonzero
*              B-splines. If ilfs has the correct value, these B-splines
*              will be
*
*              B(ilfs-k1+1,k1),B(ilfs+k1+2,k1),...,B(ilfs,k1),
*
*              and similarly in the second parameter direction.
*              This means that in Equation 1 above the matrix C can be
*              reduced to a k2xk1 matrix and the vectors of B-spline values,
*              B1(u) and B2(v), can be reduced to vectors of length
*              k1 and k2 respectively.
*
*              This notation is also valid for derivatives. The D(i,j)
*              derivative of S is given by
*
*                    D(i,j)S(u,v) = tran(D(j)B2(v) * C * D(i)B1(u),
*
*              where D(i)B1(u) denotes the vector of the i'th derivatives
*              of the B-splines accociated with t1, at the point u
*              and similarly in the second parameter direction.
*              Therefore, if in (1) the vector B1(u) is replaced with
*              the matrix DB1 with D(i)B1(u) as the i+1'st column
*              for i=0,1,...ider, and similarly for B2(v),
*              then all the required derivatives DS(u,v) are given by
*              the matrix product
*
*                        DS(u,v) = tran(DB2) * C * DB1.    (2)
*
*              Here DS(u,v) is an iderxider matrix. This is the basis
*              for the algorithm: First the matrix DB2 is computed,
*              then tran(DB2) is multiplied with C and the result stored
*              in the local array ew.
*              Finally DB1 is computed
*              and multiplied with ew and the result stored in eder.
*                               --- Knut Moerken.
*
*              Note that only the elements of DS in the TOP LEFT
*              HAND TRIANGLE are computed. This saves time and
*              space, c.f. s1424   --- Michael Floater.
*
*
* CALLS      : s1220 - Computes B-spline values and derivatives at
*                      a given point.
*              s1219 - Determines ilfs.
*              s6err    - Error handling routine 
*              s6strider - Make derivative of rational expression
*
*
* WRITTEN BY : Michael Floater, SI, 1.9.92. The old version of
*              s1421 called s1424 which wasted (a little) time
*              and space. For example, if ider=2, the derivatives
*              puuv, puvv, and puuvv were calculated (in s1424)
*              and then discarded.
*                 This new version is itself a version of s1424 but
*              calculates no unnecessary derivatives.
*              After testing with clock() this new version of s1421
*              appears to be on average about 10% faster for B-splines
*              and 20% faster for NURBS.
*              Note that for NURBS we call s6strider, a triangular
*              version of s6sratder.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              SISL_NULL tests included
* Revised by : Johannes Kaasa, SINTEF Oslo, Nov. 1995,
*              Made local copies of leftknot.
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of error.                          */
  int kn1,kn2;        /* The number of B-splines accociated with the knot
			 vectors st1 and st2.                            */
  int kk1,kk2;        /* The polynomial order of the surface in the two
			 directions.                                     */
  int kdim;           /* The dimension of the space in which the surface
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int kleft2,kleft1;  /* Local versions of ilfs and ilft which are
			 used in order to avoid the pointers.            */
  int ki,kj,kih,kjh;  /* Control variables in for loops and for stepping
			 through arrays.                                 */
  int kh,kl,kl1,kl2;  /* Control variables in for loops and for stepping
			 through arrays.                                 */
  double *st1,*st2;   /* The knot vectors of the surface. These have
			 length [kn1+kk1] and [kn2+kk2],
			 respectively.                                   */
  double *scoef;      /* The B-spline coefficients of the surface.
			 This is an array of dimension [kn2*kn1*kdim].   */
  double tt;          /* Dummy variable used for holding an array element
			 in a for loop.                                  */
  double *ebder=SISL_NULL; /* Pointer to an array of dimension
			 [max(kk1*(ider+1),kk2*(ider+1))] which will
			 contain the values and ider first derivatives of
			 the kk1 (kk2) nonzero B-splines at epar[0] (epar[1]).
			 These are stored in the following order:
			 First the value, 1. derivative etc. of the
			 first nonzero B-spline, then the same for the
			 second nonzero B-spline and so on.              */
  
  double *ew=SISL_NULL;    /* Pointer to an array of dimension [kk1*(ider+1)*kdim]
			 which will be used to store the result of the first
			 matrix multiplication in (2) above. This array is
			 initialized to all zeros.                       */
  double *sder=SISL_NULL;  /* Pointer to array used for storage of points, if
			 non rational sder points to eder, if rational sder
			 has to be allocated to make room for the homogenous
			 coordinate */
  
  double sdum1[49];   /* Arraye used for ebder */
  double sdum2[147];  /* Array used for ew */
  int knumb1;         /* Necessary size of ebder */   
  int knumb2;         /* Necessary size of ew */   

  int tot,temp;       /* Temporary variables. */
  int kinc;           /* For controlling kih.                */
  
  kleft1 = *ilfs;
  kleft2 = *ilft;
  
  /* Copy surface to local parameters.  */
  
  kn1 = ps1 -> in1;
  kn2 = ps1 -> in2;                                         
  kk1 = ps1 -> ik1;
  kk2 = ps1 -> ik2;
  st1 = ps1 -> et1;
  st2 = ps1 -> et2;
  kdim = ps1 -> idim;
  if (ps1->ikind == 2 || ps1->ikind == 4)
    {
      scoef = ps1 -> rcoef;
      kdim +=1;
      if((sder=newarray(kdim*(ider+1)*(ider+2)/2,DOUBLE)) == SISL_NULL)
         goto err101;
    }
  else
    {
      scoef = ps1 -> ecoef;
      sder = eder;  
    }
  
  /* Check the input. */
  
  if (kdim < 1) goto err102;
  if (kk1 < 1) goto err115;
  if (kn1 < kk1 || kn2 < kk2) goto err116;
  if (ider < 0) goto err178;
  if (st1[kk1-1] == st1[kk1] || st1[kn1-1] == st1[kn1]) goto err117;
  if (st2[kk2-1] == st2[kk2] || st2[kn2-1] == st2[kn2]) goto err117;  
  
  /* Allocate space for B-spline values and derivatives and one work array. */
  
  knumb1 = max(kk1*(ider+1),kk2*(ider+1));
  
  /* ONly allocate ebder if sdum1 too small */
  
  if (knumb1>49)
    { 
      if((ebder = newarray(knumb1,double)) == SISL_NULL) goto err101;
    }
  else
    {
      ebder = &sdum1[0];
      for (ki=0;ki<knumb1;ki++)
	ebder[ki] = DZERO;
    }
  
  if (ebder == SISL_NULL) goto err101;
  
  /* Only allocate ew if sdum2 too small */
  
  knumb2 = (kk1*(ider+1)*kdim);
  if (knumb2>147)
    {
      if((ew=new0array(knumb2,double)) == SISL_NULL) goto err101;
    }
  else
    { 
      ew = &sdum2[0];
      for (ki=0;ki<knumb2;ki++)
	sdum2[ki] = DZERO;
    }
  
  if (ew == SISL_NULL) goto err101;
  
  /* Set all the elements of sder to 0. */
  
  for (ki=0; ki<kdim*(ider+1)*(ider+2)/2; ki++) sder[ki] = DZERO;
  
  /* Compute the values and derivatives of the nonzero B-splines in the
     second parameter direction.                                        */
  
  s1220(st2,kk2,kn2,&kleft2,epar[1],ider,ebder,&kstat);
  
  if (kstat < 0) goto error;
  
  /* Update ilfs (ilft was updated above, in s1220). */
  
  s1219(st1,kk1,kn1,&kleft1,epar[0],&kstat);
  
  if (kstat < 0) goto error;
  
  /* Compute the first matrix product in (2) above. */
  
  /* ki steps through the appropriate kk2 rows of B-spline coefficients
     while kih steps through the B-spline value and derivatives for the
     B-spline given by ki.                                              */
  
  kih = 0;
  for (ki=kleft2-kk2+1; ki<=kleft2; ki++)
    {      
      /* kj counts through the ider+1 derivatives to be computed.
	 kjh steps through ew once for each ki to accumulate the contribution
	 from the different B-splines.
	 kl1 points to the first component of the first B-spline coefficient
	 in row no. ki of the B-spline coefficient matrix that multiplies
	 a nonzero B-spline in the first parameter direction.
	 */
      
      kjh = 0; kl1 = ki*kdim*kn1 + kdim*(kleft1-kk1+1);
      for (kj=0; kj<=ider; kj++)
	{
	  
	  /* The value of the B-spline derivative is stored in tt while
	     kl2 steps through the kdim components of all the B-spline
	     coefficients that multiplies nonzero B-splines along st1. 
	     */
	  
	  tt = ebder[kih++]; kl2 = kl1;
	  for (kl=0; kl<kdim*kk1; kl++,kjh++,kl2++)
	    {
	      ew[kjh] += scoef[kl2]*tt;
	    }
	}
    }
  
  /* Compute the values and derivatives of the nonzero B-splines in the
     first parameter direction.                                        */
  
  s1220(st1,kk1,kn1,&kleft1,epar[0],ider,ebder,&kstat);         
  
  if (kstat < 0) goto error;
  
  /* Compute the remaining matrix product. */
  
  /* kh steps through the ider+1 derivatives in the first parameter direction
     (the rows of ew if we image it as a kk1x(ider+1) matrix with each element
     a kdim dimensional vector) while kl1 steps through the elements of ew
     (again considering each element to have kdim components).                   
     */
  
  kl1 = 0;
  for (kh=0; kh<=ider; kh++)
    {      
      /* ki steps through the kk1 columns of ew (corresponding to the columns
	 of scoef that multiply nonzero B-splines along st1), while kih
	 steps through the B-spline values and derivatives for the nonzero
	 B-splines along st1 (stored in ebder).
	 */

      kinc = 0;
      for (ki=0; ki<kk1; ki++,kinc+=(ider+1))
	{
	  
	  kih = kinc;

	  /* kj counts through the ider+1 derivatives in the first
	     parameter direction (corresponding to the columns of sder).
	     kjh points to the row of sder corresponding to derivatives of
	     order kh in the second parameter direction (if sder is
	     considered a matrix with elements consisting of vectors with
	     kdim components).
	     */
	  
	  for (kj=0; kj<=ider-kh; kj++)
	    {
	      /* Find index for sder (a triangular matrix). */

              tot = kj + kh;

	      temp = ((tot * (tot+1)) >> 1) + kh;

	      kjh = temp * kdim;

	      /* Pick out the current element of ebder.
		 kl2 steps through the kdim components of the (kh,ki)
		 element of ew.
		 */
	      
	      tt = ebder[kih++];
	      kl2 = kl1;
	      for (kl=0; kl<kdim; kl++,kjh++,kl2++)
		{
		  sder[kjh] += ew[kl2]*tt;
		}
	    }
	  kl1 += kdim;
	}
    }
  
  /* Free memory. */
  
  /* If rational surface calculate the derivatives based on derivatives in
     homogenous coordinates */
  
  if (ps1->ikind == 2 || ps1->ikind == 4)
    {
      s6strider(sder,ps1->idim,ider,eder,&kstat);
      if (kstat<0) goto error;
      if(sder != SISL_NULL) freearray(sder);
    }
  
  /* Only free ew and ebder if the were allocated by newarray */
  
  if (knumb1 > 49 && ebder != SISL_NULL)
     freearray(ebder);
  if (knumb2 > 147 && ew != SISL_NULL)
     freearray(ew);
  
  /* Make cross products of tangents, if idim==3 and derivative >0 */
  
  if (ider>0 && ps1->idim ==3)
    {
      double tlen1,tlen2,tnorm,tang=(double)0.0;
      
      s6crss(eder+ps1->idim,eder+2*ps1->idim,enorm);
      
      /*  Make length of tangents and normal */
      
      tlen1 = s6length(eder+ps1->idim,ps1->idim,&kstat);
      tlen2 = s6length(eder+2*ps1->idim,ps1->idim,&kstat);
      tnorm = s6length(enorm,ps1->idim,&kstat);
      
      /*  Calculate angle between tangents */
      
      if (tlen1 != DZERO && tlen2 != DZERO && tnorm != DZERO)
        tang = tnorm/(tlen1*tlen2);
      
      if (tang == DZERO) *jstat = 2;
      else if (tang <= ANGULAR_TOLERANCE) *jstat = 1;   
      else *jstat = 0;
      goto out;
    }
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  /* Not enough memory. */
 err101: *jstat = -101;
  s6err("s1421",*jstat,kpos);
  goto out;
  
  /* kdim less than 1. */
 err102: *jstat = -102;
  s6err("s1421",*jstat,kpos);
  goto out;
  
  /* Polynomial order less than 1. */
 err115: *jstat = -115;
  s6err("s1421",*jstat,kpos);
  goto out;
  
  /* Fewer B-splines than the order. */
 err116: *jstat = -116;
  s6err("s1421",*jstat,kpos);
  goto out;
  
  /* Error in knot vector.
     (The first or last interval of one of the knot vectors is empty.) */
 err117: *jstat = -117;
  s6err("s1421",*jstat,kpos);
  goto out;
  
  /* Illegal derivative requested. */
 err178: *jstat = -178;
  s6err("s1221",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
 error:  *jstat = kstat;
  s6err("s1421",*jstat,kpos);
  goto out;
  
 out:
    *ilfs = kleft1;
    *ilft = kleft2;
    return;
}
