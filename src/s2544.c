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
 * $Id: s2544.c,v 1.6 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S2544

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2544(SISLSurf *surf, int ider, int iside1, int iside2, double parvalue[],
      int *leftknot1,int *leftknot2, double norcurv[], int *jstat)
#else
 void s2544(surf, ider, iside1, iside2, parvalue, leftknot1, leftknot2, norcurv,
	    jstat)
      SISLSurf *surf;
      int    ider;
      int    iside1;
      int    iside2;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double norcurv[];
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the Normal curvature of a Surface for given
*                  values (u,v) = (parvalue[0],parvalue[1]), in the
*                  direction (parvalue[2],parvalue[3])
*                  where:
*
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                  
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0 and 1.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*          iside1   - Indicator telling if the derivatives in the first
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*          iside2   - Indicator telling if the derivatives in the second
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*      parvalue     - Parameter-value at which to evaluate pluss the direction
*                     Dimension of parvalue is 4.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     first parameter direction where parvalue[0] is found,
*                     that is:
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the routine.
*
*     leftknot2     - Pointer to the interval in the knot vector in the
*                     second parameter direction where parvalue[1] is found,
*                     that is:
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the routine.
*
*  OUTPUT       :
*     norcurv      - Normal curvature and derivatives of normal curvature
*                    of the surface in (u,v) =(parvalue[0],parvalue[1])
*                    in the direction (parvalue[2],parvalue[3]).
*        jstat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The Normal curvature is given by
*
*                      <Xuu,N>d1*d1 + 2<Xuv,N>d1*d2 + <Xvv,N>d2*d2
*                      --------------------------------------------
*                      <Xu,Xu>d1*d1 + 2<Xu,Xv>d1*d2 + <Xv,Xv>d2*d2
*                  
*                  
*		   Where (d1,d2) is the normalised parameter direction
*                  (parvalue[2],parvalue[3]) and the other numbers are
*                  the factors for the first and second fundamental forms.
*
*                  The routine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate. 
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1422() and s2513().
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Normal Curvature.
*                    The routine returns jstat = 2.
*               (ii) If the surface is close to degenerate, the calculations
*                    can be numerical instable. The routine returns
*                    jstat = 1.
*              (iii) If the surface is Cr the curvature calculated is C(r-2).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.
*
*
* WRITTEN BY :    Ulf J Krystad, SINTEF, Oslo, Norway.            Date: 1995-1
* REVISED BY :    Johannes Kaasa, SINTEF, Oslo, Norway.           Date: 1995-9
* REWRITTEN BY :  Johannes Kaasa, SINTEF, Oslo, Norway.           Date: 1995-11
*                 (Added first derivative of normal curvature).
******************************************************************************
*/
{
  int kwarn = 0;      	 /* Local staus variable(warning).                  */
  double derive[30];     /* Array containing the computed derivatives.      */
  double normal[3];      /* Array containing the computed normalvektor.     */
  double fundform[10];   /* The coefficients of the fundamental forms.
			    The sequence is: E, F, G, e, f, g, P, Q, S, T.  */
  double d1, d2;         /* The normalised parameter direction
			    (parvalue[2],parvalue[3]) and their length      */
  double temp1, temp2;   /* Temporary values				    */
  double length;         /* Square of normal length.                        */
  double a, b, c, d;     /* Utility coefficients.                           */
  double sigma;          /* Utility coefficient.                            */
  double alpha;          /* Utility coefficient.                            */
  double beta;           /* Utility coefficient.                            */
  double gamma;          /* Utility coefficient.                            */
  double delta;          /* Utility coefficient.                            */
  double epsilon;        /* Utility coefficient.                            */
  double mu;             /* Utility coefficient.                            */
  double P, Q, S, T;     /* Third order coefficients.                       */
  double D;              /* Utility coefficient.                            */
  double H;              /* The mean curvature.                             */
  double K;              /* The Gaussian curvature.                         */
  double k1;             /* Max. principal curvature.                       */
  double k2;             /* Min. principal curvature.                       */
  double theta;          /* Angle in the tangent plane.                     */
  double phi;            /* Angle in the tangent plane.                     */
  double psi;            /* Angle in the tangent plane.                     */
  double tanglenA;       /* Length of tangentA.                             */
  double tanglenB;       /* Length of tangentB.                             */
  double asin_result;    /* Result of asin.                                 */
  double sin_contrib;    /* Contribution from sinus.                        */
  double cos_contrib;    /* Contribution form cosinus.                      */
  double max_dir[2];     /* Max. direction of the principal curvature k1.   */
  double min_dir[2];     /* Min. direction of the principal curvature k2.   */
  double tangentA[3];    /* Tangent to the surface.                         */
  double tangentB[3];    /* Tangent to the surface.                         */
  double cross[3];       /* Cross product vector.                           */
  /* ______________________________________________________________________ */
  

  if (ider < 0 || ider > 1) goto err178;
  
  length = sqrt(parvalue[2]*parvalue[2] + parvalue[3]*parvalue[3]);
  if (length < REL_PAR_RES) goto err174;
  d1 = parvalue[2]/length;
  d2 = parvalue[3]/length;


  if (surf == SISL_NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

    s1422(surf, ider+2, iside1, iside2, parvalue, leftknot1, leftknot2, derive,
	  normal, jstat);
    if (*jstat > 0) kwarn = *jstat;

    if (*jstat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (*jstat != 2) /* The surface is not degenerate */
    {
       /* Find factors in fundamental form */
       s2513(surf, 0, ider+2, 1, derive, normal, fundform, jstat);
       
       if (*jstat < 0)
	  goto error;
       
       temp1    = fundform[0]*d1*d1 + 2*fundform[1]*d1*d2 + fundform[2]*d2*d2;
       if (temp1 < REL_PAR_RES) goto err174;
       temp2    = fundform[3]*d1*d1 + 2*fundform[4]*d1*d2 + fundform[5]*d2*d2;
       
       norcurv[0] = temp2/temp1;
       
       if (ider > 0)
       {
	  
	  /* Calculate the derivative of the normal curvature. */
	  
	  if (surf->idim == 3)
	     length = normal[0]*normal[0] + normal[1]*normal[1] +
		normal[2]*normal[2];
	  else if (surf->idim == 1)
	     length = 1. + derive[1]*derive[1] + derive[2]*derive[2];
	  sigma = sqrt(length);
	  
	  /* Find the coefficients in the expression. */
	  
	  alpha   = s6scpr(&derive[3*(surf->idim)], &derive[surf->idim],
			   surf->idim)/length;
	  beta    = s6scpr(&derive[3*(surf->idim)], &derive[2*(surf->idim)],
			   surf->idim)/length;
	  gamma   = s6scpr(&derive[4*(surf->idim)], &derive[surf->idim],
			   surf->idim)/length;
	  delta   = s6scpr(&derive[4*(surf->idim)], &derive[2*(surf->idim)],
			   surf->idim)/length;
	  epsilon = s6scpr(&derive[5*(surf->idim)], &derive[surf->idim],
			   surf->idim)/length;
	  mu      = s6scpr(&derive[5*(surf->idim)], &derive[2*(surf->idim)],
			   surf->idim)/length;
	  
	  a = fundform[1]*fundform[4] - fundform[2]*fundform[3];
	  b = fundform[1]*fundform[3] - fundform[0]*fundform[4];
	  c = fundform[1]*fundform[5] - fundform[2]*fundform[4];
	  d = fundform[1]*fundform[4] - fundform[0]*fundform[5];
	  
	  P = fundform[6] + 3*(a*alpha + b*beta);
	  Q = fundform[7] + c*alpha + d*beta + 2*a*gamma + 2*b*delta;
	  S = fundform[8] + 2*c*gamma + 2*d*delta + a*epsilon + b*mu;
	  T = fundform[9] + 3*(c*epsilon + d*mu);
	  
	  /* Calculate the principal curvatures. */
	  
	  s2543(surf, 0, derive, normal, &k1, &k2, max_dir, min_dir, jstat);
	  
	  if (*jstat < 0)
	     goto error;
	  
	  H = (k1 + k2)/2.;
	  K = k1*k2;
	  
	  D = sqrt(H*H - K);
	  if (fabs(D) < REL_PAR_RES)
	  {
	     /* Umbilical point, derivative is zero. */
	     
	     norcurv[1] = 0.;
	     goto out;
	  }
	  
	  /* Calculate the angles theta, phi and psi. */
	  
	  if (surf->idim == 3)
	  {
	     tangentA[0] = max_dir[0]*derive[3] + max_dir[1]*derive[6];
	     tangentA[1] = max_dir[0]*derive[4] + max_dir[1]*derive[7];
	     tangentA[2] = max_dir[0]*derive[5] + max_dir[1]*derive[8];
	     tangentB[0] = parvalue[2]*derive[3] + parvalue[3]*derive[6];
	     tangentB[1] = parvalue[2]*derive[4] + parvalue[3]*derive[7];
	     tangentB[2] = parvalue[2]*derive[5] + parvalue[3]*derive[8];
	  }
	  else
	  {
	     tangentA[0] = max_dir[0];
	     tangentA[1] = max_dir[1];
	     tangentA[2] = max_dir[0]*derive[1] + max_dir[1]*derive[2];
	     tangentB[0] = parvalue[2];
	     tangentB[1] = parvalue[3];
	     tangentB[2] = parvalue[2]*derive[1] + parvalue[3]*derive[2];
	  }
	  
	  tanglenA = s6length(tangentA, 3, jstat);
	  if (*jstat < 0)
	     goto error;
	  tanglenB = s6length(tangentB, 3, jstat);
	  if (*jstat < 0)
	     goto error;

	  theta = acos(max(-1., min(1., s6scpr(tangentA, tangentB, 3)/
				     (tanglenA*tanglenB))));
	  
	  /* Check rotational direction. */
	  
	  s6crss(normal, tangentA, cross);
	  if (s6scpr(cross, tangentB, 3) < 0.)
	     theta = TWOPI - theta;
	  
	  /* Calculate phi. */

	  phi = acos(max(-1., min(1., (fundform[5] - H*fundform[2])/
				   (fundform[2]*D))));
	  asin_result = asin(max(-1.,min(1., (fundform[4]*fundform[2] - 
			fundform[5]*fundform[1])/(sigma*fundform[2]*D))));
	  if (asin_result < 0)
	     phi = TWOPI - phi;
	  if (fabs(phi - asin_result) > ANGULAR_TOLERANCE &&
	      fabs(PI - phi - asin_result) > ANGULAR_TOLERANCE &&
	      fabs(phi - TWOPI - asin_result) > ANGULAR_TOLERANCE)
	     goto err180;

	  phi /= 2.;
	  
	  /* Calculate psi. */

          psi = acos(max(-1., min(1., sigma/sqrt(fundform[0]*fundform[2]))));
          asin_result = asin(max(-1.,min(1., fundform[1]/
					 sqrt(fundform[0]*fundform[2]))));
	  if (asin_result < 0)
	     psi = TWOPI - psi;
	  if (fabs(psi - asin_result) > ANGULAR_TOLERANCE &&
	      fabs(PI - psi - asin_result) > ANGULAR_TOLERANCE &&
	      fabs(psi - TWOPI - asin_result) > ANGULAR_TOLERANCE)
	     goto err180;
	  
	  /* Calculate the derivative. */
	  
	  sin_contrib = sin(theta + phi);
	  cos_contrib = cos(theta + phi + psi);
	  
	  norcurv[1] = (P*sqrt(fundform[2]*fundform[2]*fundform[2])*
			sin_contrib*sin_contrib*sin_contrib
			+ 3*Q*fundform[2]*sqrt(fundform[0])*cos_contrib*
			sin_contrib*sin_contrib
			+ 3*S*sqrt(fundform[2])*fundform[0]*cos_contrib*
			cos_contrib*sin_contrib
			+ T*sqrt(fundform[0]*fundform[0]*fundform[0])*
			cos_contrib*cos_contrib*cos_contrib)/(sigma*length);
       }
       
    }
    else if (*jstat == 2) /* The surface is degenerated. */
    {
      norcurv[0] = 0.0;
      goto war002;
    }

  }


  /* Successful computations  */

  *jstat = kwarn;
  goto out;


  /* ____________________________________________________________________ */
  /*                           ERROR EXIT				  */
  /* ____________________________________________________________________ */
  
   /* The surface is degenerated at (u,v) */
war002:
  *jstat = 2;
  goto out;

  /* Error. Input (surface) pointer is SISL_NULL. */
err150:
  *jstat = -150;
  s6err("s2544", *jstat, 0);
  goto out;

  /* Degenerate condition. */
err174:
  *jstat = -174;
  s6err("s2544",*jstat,0);
  goto out;
  
  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2544",*jstat,0);
  goto out;

  /* Problems in angle calculation. */
err180:
  *jstat = -180;
  s6err("s2544",*jstat,0);
  goto out;
  
  /* Error in lower level routine.  */
error:
  s6err("s2544",*jstat,0);
  goto out;


  /* ____________________________________________________________________ */
  /*                        THE ONE AND ONLY EXIT			  */
  /* ____________________________________________________________________ */
out:

  return;

}
