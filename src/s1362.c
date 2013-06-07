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
 * $Id: s1362.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1362

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1362(SISLCurve *pc1,double aoffset,double enorm[],int idim,
	   int ider,double ax,int *ileft,double eder[],int *jstat)
#else
void s1362(pc1,aoffset,enorm,idim,ider,ax,ileft,eder,jstat)
     SISLCurve *pc1;
     double       aoffset;
     double       enorm[];
     int          idim;
     int          ider;
     double       ax;
     int          *ileft;
     double       eder[];
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and ider<=2 first derivatives of
*              the offset curve to the B-spline curve pointed to by pc1,
*              at the point with parameter value ax.
*
*
*
* INPUT      : pc1    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              aoffset- Offset distance
*              enorm  - Only used when idim==3
*              idim   - Dimension of the space the curve lies in 2 or 3
*              ider   - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              ax     - The parameter value at which to compute
*                       position and derivatives.
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
* METHOD     : The offset curve is made in the following way:
*
*              q(u) = p(u) + aoffset N(u)
*
*              For idim==2: N(u) = (-y'(u),x'(u))
*
*              For idim==3: N(u) = p'(u)xenorm
*
*              The derivatives of these expressions are calculated:
*
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1221, s6length, s6crss, s6scpr, s6err
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, October 1988
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of the error.                      */
  int kder=ider+1;    /* Derivatives necessary to calculate              */
  int ki,kj;          /* Loop variable                                   */
  double sder[12];    /* Storage of derivatives                          */
  double snorm[3];    /* Normal to tangent vector                        */
  double sdnorm[3];   /* Derivative of normal to tangent                 */
  double sddnorm[3];  /* Second derivative to normal of tangent          */
  double tderl;       /* Length of the first derivative                  */
  double tl;          /* The length of the derivative vector             */
  double *spnt;       /* Pointer to derivative                           */
  double *start;      /* Pointer to derivative                           */
  
  if (idim !=2 && idim != 3) goto err105;
  
  if (DEQUAL(aoffset,DZERO))
    {
      s1221(pc1,ider,ax,ileft,eder,&kstat);
      if (kstat<0) goto error;
    }
  else
    {
      s1221(pc1,kder,ax,ileft,sder,&kstat);
      if (kstat<0) goto error;
      
      tderl = s6length(sder+idim,idim,&kstat);
      if (DEQUAL(tderl,DZERO)) tderl = (double)1.0;
      
      /* The tangent length might be very different from 1. scale it and
	 higher order derivatives to give tangent length one, this may give
	 over or underflow in the calculations */
      
      
      for (ki=1,start=sder+idim ; ki<=kder ; ki++,start+=idim)
        {
	  for (kj=0,spnt=start ; kj<idim*(kder+1-ki) ; kj++,spnt++)
            {
	      *spnt /= tderl;
            }
        }
      
      
      if (idim==2)
        {
	  snorm[0]   = -sder[3];
	  snorm[1]   =  sder[2];
	  if(ider>=1)
	    {                                                           
	      sdnorm[0]  = -sder[5];
	      sdnorm[1]  =  sder[4];
	      if(ider>=2)
		{
		  sddnorm[0] = -sder[7];
		  sddnorm[1] =  sder[6];
		}
	    }
        }
      else
        {
	  s6crss(sder+idim,enorm,snorm);
	  if (ider>=1)
	    {
	      s6crss(sder+2*idim,enorm,sdnorm);
	      if (ider>=2)
		s6crss(sder+3*idim,enorm,sddnorm);
	    }
        }
      
      tl = s6length(snorm,idim,&kstat);
      if (DEQUAL(tl,DZERO)) tl = (double)1.0;
      
      /* Calculate position on offset curve */
      
      for (ki=0;ki<idim;ki++)
        eder[ki] = sder[ki] + aoffset*snorm[ki]/tl;
      
      if (ider>=1)
        {
	  double tpndpn;         /* The scalar product: snorm sdnorm */
	  double tl3=tl*tl*tl;   /* tl**3 */
	  
	  tpndpn = s6scpr(snorm,sdnorm,idim);
	  
	  
	  /*  Calculate derivative of offset curve */
	  for (ki=0;ki<idim;ki++)
            {
	      eder[idim+ki] = sder[idim+ki] + aoffset*(sdnorm[ki]/tl
			    - snorm[ki]*tpndpn/tl3);
            }
	  if (ider>=2)
            {
	      double tl5= tl3*tl*tl; /* tl**5 */
	      double tpnddpn;         /* The scalar product: snorm sddnorm */
	      double tdpndpn;         /* The scalar product: sdnorm sdnorm */
	      
	      tpnddpn = s6scpr(snorm,sddnorm,idim);
	      tdpndpn = s6scpr(sdnorm,sdnorm,idim); 
	      
	      
	      /*      Calculate derivative of offset curve */
	      for (ki=0;ki<idim;ki++)
                {
		  eder[2*idim+ki] = sder[2*idim+ki] + aoffset*(sddnorm[ki]/tl
				  - (double)2.0*sdnorm[ki]*tpndpn/tl3
				  - snorm[ki]*(tdpndpn+tpnddpn)/tl3
				  + (double)3.0*snorm[ki]*tpndpn*tpndpn/tl5);
                }
            }
        }
      
      
      /* Scale the derivatives to match the original parametrization,
	 the kder=ider+1 derivative is not given as output  */
      
      for (ki=1,start=eder+idim ; ki<kder ; ki++,start+=idim)
        {
	  for (kj=0,spnt=start ; kj<idim*(kder-ki) ; kj++,spnt++)
            {
	      *spnt *= tderl;
            }
        }
    }
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  
  /* idim not 2 0r 3 */
 err105: *jstat = -105;
  s6err("s1362",*jstat,kpos);
  goto out;
  
  
  /* Error in lower level routine.  */
  
 error:  *jstat = kstat;
  s6err("s1362",*jstat,kpos);
  goto out;
  
 out: return;
}
