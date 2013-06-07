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
 * $Id: s9iterate.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */
#define S9ITERATE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9iterate(double epoint[],double epnt1[],double epnt2[],double epar1[],
	       double epar2[],SISLSurf *psurf1,SISLSurf *psurf2,double astep,
	       double aepsge,double gpnt1[],double gpnt2[],double gpar1[],
	       double gpar2[],int *jstat)
#else
void s9iterate(epoint,epnt1,epnt2,epar1,epar2,psurf1,psurf2,astep,aepsge,
               gpnt1,gpnt2,gpar1,gpar2,jstat)
     double epoint[];
     double epnt1[];
     double epnt2[];
     double epar1[];
     double epar2[];
     SISLSurf   *psurf1;
     SISLSurf   *psurf2;
     double astep;
     double aepsge;
     double gpnt1[];
     double gpnt2[];
     double gpar1[];
     double gpar2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To iterate to an intersection point between two surfaces
*              and a plane.
*
*
*
* INPUT      : epoint - Array containing parts of plane description.
*                       epoint[0:2] contains a position value.
*                       epoint[3:5] contains the normal to the plane
*                       A point in the plane is defined by
*                       epoint[0:2] + astep*epoint[3:5]
*              epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in first surface
*              epnt2  - 0-2 Derivatives + normal of start point for
*                       iteration in second surface
*              epar1  - Parameter pair of start point in first surface
*              epar2  - Parameter pair of start point in second surface
*              psurf1 - Description of first surface
*              psurf2 - Description of second surface
*              astep  - Step length
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in first surface
*              gpnt2  - 0-2 Derivatives + normal of result of iteration
*                       in second surface
*              gpar1  - Parameter pair of result of iteration in first surface
*              gpar2  - Parameter pair of result of iteration in second
*                       surface
*              jstat  - status messages  
*                       = 2      : Iteration diverged or to many iterations
*                       = 1      : iteration converged, singular point found
*                       = 0      : ok, iteration converged
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, June-1988
* Revised by : Tor Dokken, SI, OSLO, Norway, 24-Feb-1989
*              Prepared for degenerate points
* Revised by : Tor Dokken, SI, Oslo, Norway, 3-April-1989
*              Correct handling of small determinats
*
*********************************************************************
*/
{
  int ki;                 /* Variable used in loop                          */
  int kcont;              /* Indicator telling if iteration is not finished */
  int kder = 2;           /* Derivative indicator                           */
  int klfu=0;             /* Pointer into knot vector                       */
  int klfv=0;             /* Pointer into knot vector                       */
  int klfs=0;             /* Pointer into knot vector                       */
  int klft=0;             /* Pointer into knot vector                       */
  int kstat;              /* Status variable                                */
  int knbit;              /* Counter for number of iterations               */
  int kdim = 3;           /* Set dimension to 3                             */
  int kmaxit = 100;       /* Maximal number of iterations allowed           */
  int kpos=1;             /* Position indicator ofr errors                  */
  double spoint[3];       /* SISLPoint in intersection plane                    */
  double *snorm;          /* Pointer to normal vector of intersection plane */
  double *sp,*spu,*spv,*spn; /* Pointers into gpnt1                         */
  double *sq,*sqs,*sqt,*sqn; /* Pointers into gpnt2                         */
  double ta11,ta12,ta21;  /* Variables used in equation systems             */
  double ta22,tb1,tb2;    /* Variables used in equation systems             */
  double sdiff[3];        /* Difference between two vectors                 */
  double tdum1,tdum2;     /* Dummy variables                                */
  double tdum3,tdum;      /* Dummy variables                                */
  double tdist;           /* Distance betweentwo points in iteration        */
  
  
  
  /* Make description of intersection plane */
  
  for (ki=0;ki<3;ki++)
    {
      spoint[ki] = epoint[ki] + astep*epoint[ki+3];
    }
  
  snorm = epoint + 3;
  
  /* Copy input variables to output variables */
  
  memcopy(gpnt1,epnt1,21,DOUBLE); 
  memcopy(gpnt2,epnt2,21,DOUBLE);
  memcopy(gpar1,epar1,2,DOUBLE); 
  memcopy(gpar2,epar2,2,DOUBLE);
  
  /* At the start of the iteration the two point gpnt1 and gpnt2 might be
     very close since we in most cases start from a point on the intersection
     curve. */
  
  /* Set a number of local pointers that ar used often */
  sp  = gpnt1;
  spu = gpnt1 + 3;
  spv = gpnt1 + 6;
  spn = gpnt1 + 18;
  sq  = gpnt2;
  sqs = gpnt2 + 3;
  sqt = gpnt2 + 6;
  sqn = gpnt2 + 18;
  
  kcont = 1;
  knbit = 0;
  
  while (kcont)
    
    {
      
      /* Put a parametric representation of the tangent 
	 plane of surface 1 into
	 the implicit representation of the tangent 
	 plane of surface 2 and also
	 into the implicit representation of 
	 the intersection plane */
      
      ta11 = s6scpr(spu,sqn,kdim);
      ta12 = s6scpr(spv,sqn,kdim);
      ta21 = s6scpr(spu,snorm,kdim);
      ta22 = s6scpr(spv,snorm,kdim);
      
      s6diff(sq,sp,kdim,sdiff);
      tb1  = s6scpr(sdiff,sqn,kdim);
      
      tdum = MAX(fabs(ta11),fabs(ta12));
      tdum = MAX(tdum,fabs(tb1));
      if (tdum == DZERO) tdum = (double)1.0;
      ta11 /= tdum;
      ta12 /= tdum;
      tb1  /= tdum;
      
      s6diff(spoint,sp,kdim,sdiff);
      tb2  = s6scpr(sdiff,snorm,kdim);
      
      tdum = MAX(fabs(ta21),fabs(ta22));
      tdum = MAX(tdum,fabs(tb2));
      if (tdum == DZERO) tdum = (double)1.0;
      ta21 /= tdum;
      ta22 /= tdum;
      tb2  /= tdum;
      
      /* Calculate determinant of equation system */
      
      tdum1 = ta11*ta22 - ta12*ta21;
      tdum  = MAX(fabs(ta11),fabs(ta22));
      tdum  = MAX(fabs(ta12),tdum);
      tdum  = MAX(fabs(ta21),tdum);
      
      if (DEQUAL((tdum+tdum1),tdum)) tdum1 =DZERO;
      
      
      /* If tdum1 = 0.0, then the equation system is singular, 
	 iteration not possible */

      if (DNEQUAL(tdum1,DZERO))
        {
	  gpar1[0] += (tb1*ta22-tb2*ta12)/tdum1;
	  gpar1[1] += (ta11*tb2-ta21*tb1)/tdum1;
        }
      
      /* Put a parametric representation of the 
	 tangent plane of surface 2 into
	 the implicit representation of the 
	 tangent plane of surface 1 and also
	 into the implicit representation 
	 of the intersection plane */
      
      ta11 = s6scpr(sqs,spn,kdim);
      ta12 = s6scpr(sqt,spn,kdim);
      ta21 = s6scpr(sqs,snorm,kdim);
      ta22 = s6scpr(sqt,snorm,kdim);
      
      s6diff(sp,sq,kdim,sdiff);
      tb1  = s6scpr(sdiff,spn,kdim);
      
      s6diff(spoint,sq,kdim,sdiff);
      tb2  = s6scpr(sdiff,snorm,kdim);
      
      /*Calculate determinant of equation system */

      tdum2 = ta11*ta22 - ta12*ta21;
      
      tdum2 = ta11*ta22 - ta12*ta21;
      tdum  = MAX(fabs(ta11),fabs(ta22));
      tdum  = MAX(fabs(ta12),tdum);
      tdum  = MAX(fabs(ta21),tdum);
      
      if (DEQUAL((tdum+tdum2),tdum)) tdum2 =DZERO;
      
      /* If tdum2 = 0.0, then the equation system is singular, 
	 iteration not possible */

      if (DNEQUAL(tdum2,DZERO))
        {
	  gpar2[0] += (tb1*ta22-tb2*ta12)/tdum2;
	  gpar2[1] += (ta11*tb2-ta21*tb1)/tdum2;
        }
      
      /* Calculate values of new points */
      
      s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+18,&kstat); 
      if (kstat<0) goto error;
      
      /* If the surface normal has zero length no use in continuing */
      
      if (kstat == 2) goto war02;
      
      s1421(psurf2,kder,gpar2,&klfs,&klft,gpnt2,gpnt2+18,&kstat); 
      if (kstat<0) goto error;
      
      /* If the surface normal has zero length no use in continuing */
      
      if (kstat == 2) goto war02;
      
      /* Make difference between the two points, 
	 and calculate length of difference */
      s6diff(gpnt1,gpnt2,kdim,sdiff);
      tdum3 = s6length(sdiff,kdim,&kstat);
      if (kstat==0) 
        {
	  /* Length is zero iteration has converged */

	  kcont = 0;
        }
      
      if (knbit==0)
        {
	  /* First iteration intitate distance variable, if the equation
	     systems were not singular */

	  if (DEQUAL(tdum1,DZERO) || DEQUAL(tdum2,DZERO)) goto war02;
	  tdist = tdum3;
	  knbit = 1;
        }
      else
        {
	  /* More than one iteration done, stop if distance is not decreasing.
	     Then decide if we converge distance between the points is within
	     the tolerance and the last step had singular or none singular
	     equation systems. */

	  knbit = knbit + 1;
	  if (tdum3>=tdist)
            {
	      /* Distance is not decreasing */
	      if (tdist <= aepsge)
                {
		  /* Distance within tolerance */
		  if (DEQUAL(tdum1,DZERO) || DEQUAL(tdum2,DZERO))
                    {
		      /* Singular equation system */
		      goto war01;
                    }
		  else
                    {
		      /* Nonsingular equation system */
		      goto war00;
                    }
                }
	      else
                {
		  /* Distance is not within tolerance, divergence */
		  goto war02;
                }
            }
	  /*      Distance still decreasing */
	  
	  tdist = tdum3;
        }
      
      /*  Make sure that not to many iteration are being done */
      if (knbit > kmaxit) goto war02;
    }
  
  
  /* Iteration converged */
 war00:
  
  *jstat = 0;
  goto out;
  
  /* Iteration converged, singular point found */
 war01: 
  *jstat = 1;
  goto out;
  
  /* To many iterations or iteration diverging */
 war02: 
  *jstat = 2;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : 
    *jstat = kstat;
  s6err("s9iterate",*jstat,kpos);
  goto out;
  
 out:
  return;
}
