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
 * $Id: ev_cv_off.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define EV_CV_OFF

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   ev_cv_off(SISLCurve *pc1,int ider,double ax,int *ileft,
	       double offset,double eder[],int *jstat)
#else
void ev_cv_off(pc1,ider,ax,ileft,offset,eder,jstat)
     SISLCurve *pc1;
     int ider;
     double ax;
     int *ileft;
     double offset;
     double eder[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and ider first derivatives of the
*              offset of the 2D B-spline curve pointed to by pc1, at 
*              the point with parameter value ax. At the moment we only
*              calculate the first derivative.
*
*
*
* INPUT      : pc1    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              ider   - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              ax     - The parameter value at which to compute
*                       position and derivatives.
*              offset - The offset distance.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*                          
*                          et[ileft] <= ax < et[ileft+1]
* 
*                        should hold. (If ax == et[in] then ileft should
*                        be in-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : eder   - Double array of dimension [(ider+1)*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*                       (The C declaration of eder as a two dimensional array
*                       would therefore be eder[ider+1,idim].)
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : We evaluate the original curve with one derivative higher,
*              and evaluate the expression for the offset.
*
* REFERENCES :
*
*-
* CALLS      : s1220    - Computes B-spline values and derivatives at
*                         a given point.
*              s6err    - Error handling routine 
*              s6ratder - Make derivative of rational expression
*
* WRITTEN BY : Johannes Kaasa, SI, March 1992. Based on s1221. 
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of the error.                      */
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdim;           /* The dimension of the space in which the curve
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int kleft;          /* Local version of ileft which is used in order to
			 avoid the pointer.                              */
  int kder;           /* Local version of ider+1. Since derivatives of 
			 order higher than kk-1 are all zero, we set
			 kder = min(kk-1,ider+1).                        */
  int kind;           /* Type of curve                                   */
  int ki,kj,kih,kjh;  /* Control variables in for loops and for stepping
			 through arrays.                                 */
  int kl,kl1,kl2;     /* Control variables in for loops and for stepping
			 through arrays.                                 */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  double *scoef;      /* Pointer to the first element of the curve's
			 B-spline coefficients. This is assumed to be an
			 array with [kn*kdim] elements stored in the
			 following order:
			 First the kdim components of the first B-spline
			 coefficient, then the kdim components of the
			 second B-spline coefficient and so on.          */
  double tt;          /* Dummy variable used for holding an array element
			 in a for loop.                                  */
  double *ebder=SISL_NULL; /* Pointer to an array of dimension [kk*(ider+1)]
		       which will contain the values and ider first derivatives
			 of the kk nonzero B-splines at ax.
			 These are stored in the following order:
			 First the value, 1. derivative etc. of the
			 first nonzero B-spline, then the same for the
			 second nonzero B-spline and so on.              */
  double *sder=SISL_NULL;  /* Pointer to array used for storage of points, if
			 non rational sder points to eder, if rational sder
			 has to be allocated to make room for the homogenous
			 coordinate */
  double hlpnorm;     /* Length of the cross normal to the curve.        */
  
  /* Copy curve attributes to local parameters.  */
  
  kn = pc1 -> in;
  kk = pc1 -> ik;
  st = pc1 -> et;
  kdim = pc1 -> idim;
  kind = pc1 ->ikind;
  
  if (kind == 2 || kind == 4)
    {
      scoef = pc1 -> rcoef;
      kdim +=1;
    }
  else
      scoef = pc1 -> ecoef;

  sder = newarray(kdim*(ider+2),DOUBLE); 
  if (sder==SISL_NULL) goto err101;
  
  /* Check the input. */
  
  if (pc1->idim != 2) goto err102;
  
  if (kk < 1) goto err110;
  
  if (kn < kk) goto err111;
  
  if (st[kk-1] == st[kk] || st[kn-1] == st[kn]) goto err112;
  
  if (ider < 0 || ider > 2) goto err178;
  
  kder = min(kk-1,ider+1);
  
  /* Allocate space for B-spline values and derivatives. */
  
  ebder = newarray(kk*(kder+1),double);
  if (ebder == SISL_NULL) goto err101;
  
  /* Set all the elements of sder to 0. */
  
  for (ki=0; ki<(ider+2)*kdim; ki++) sder[ki] = DZERO;
  
  /* Compute the values and derivatives of the nonzero B-splines and
     update ileft if necessary.                                      */
  
  s1220(st,kk,kn,ileft,ax,kder,ebder,&kstat);
  
  if (kstat < 0) goto error;
  
  kleft = *ileft;
  
  /* Evaluate the original curve with one derivative higher. */
  
  /* ki steps through the appropriate kk B-spline coefficients while kih steps
     through the B-spline value and derivatives for the B-spline given by ki.*/
  
  kih = 0;
  for (ki=kleft-kk+1; ki<=kleft; ki++)
    {
      
      /* kj counts through the kder+1 derivatives to be computed.
	 kjh steps through sder once for each ki to accumulate the contribution
	 from the different B-splines.
	 kl1 points to the first component of B-spline coefficient no. ki. */
      
      kjh = 0; kl1 = kdim*ki;
      for (kj=0; kj<=kder; kj++)
	{
	  
	  /* The value of the B-spline derivative is stored in tt while
	     kl2 steps through the idim components of this B-spline
	     coefficient. */
	  
	  tt = ebder[kih++]; kl2 = kl1;
	  for (kl=0; kl<kdim; kl++,kjh++,kl2++)
	    {
	      sder[kjh] += scoef[kl2]*tt;
	    }
	}
    }
  
  /* If rational curve calculate the derivatives based on derivatives in
     homogenous coordinates */
  
  if (kind == 2 || kind == 4)
    {
      s6ratder(sder,pc1->idim,ider,sder,&kstat);
      if (kstat<0) goto error;
    }
  
  /* Evaluate the offset expression. */
  
  hlpnorm = sqrt(sder[2]*sder[2] + sder[3]*sder[3]);
  if (hlpnorm < REL_COMP_RES) goto err113;
  
  eder[0] = sder[0] + offset*(-sder[3]/hlpnorm);
  eder[1] = sder[1] + offset*(sder[2]/hlpnorm);
  if (ider > 0)
    {
       eder[2] = sder[2] + offset*(sder[2]*sder[3]*sder[4] - sder[2]
				   *sder[2]*sder[5])/pow(hlpnorm,3); 
       eder[3] = sder[3] + offset*(sder[3]*sder[3]*sder[4] - sder[2]
				   *sder[3]*sder[5])/pow(hlpnorm,3);
    }
  if (ider > 1)
    {
       eder[4] = sder[4] + offset*((pow(sder[2],3)*sder[3] 
		 + sder[2]*pow(sder[3],3))*sder[6] - (pow(sder[2],4) 
		 + (sder[2]*sder[2])*(sder[3]*sder[3]))*sder[7] 
		 + (pow(sder[3],3) - 2*sder[3]*(sder[2]*sder[2]))*
		 (sder[4]* sder[4]) + 2*(pow(sder[2],3) 
		 - 2*sder[2]*(sder[3]*sder[3]))*(sder[4]*
                 sder[5]) + 3*(sder[2]*sder[2])*sder[3]*sder[5]*sder[5])/
	         pow(hlpnorm,5);
       eder[5] = sder[5] + offset*(((sder[2]*sder[2])*(sder[3]*sder[3]) 
		 + pow(sder[3],4))*sder[6] - (pow(sder[2],3)*sder[3] 
		 + sder[2]*pow(sder[3],3))*sder[7] -3*sder[2]
		 *(sder[3]*sder[3])*(sder[4]*sder[4]) + 2*(2*(sder[2]
		 *sder[2]*sder[3]) - pow(sder[3],3))*(sder[4]*sder[5]) + 
		 (2*sder[2]*(sder[3]*sder[3]) - pow(sder[2],3))
		  *(sder[5]*sder[5]))/pow(hlpnorm,5);
    }
  
  /* Free memory. */
  
  freearray(ebder);
  freearray(sder);
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;

/* Not enough memory. */
 err101: *jstat = -101;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* kdim less than 1. */
 err102: *jstat = -102;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Polynomial order less than 1. */
 err110: *jstat = -110;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Fewer B-splines than the order. */
 err111: *jstat = -111;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Error in knot vector.
   (The first or last interval of the knot vector is empty.) */
 err112: *jstat = -112;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Singularity in the original evaluation point.
   (The original curve tangent is zero in this point.) */
 err113: *jstat = -113;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Illegal derivative requested. */
 err178: *jstat = -178;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

/* Error in lower level routine.  */

 error:  *jstat = kstat;
         s6err("ev_cv_off",*jstat,kpos);
         goto out;

out: return;
}
