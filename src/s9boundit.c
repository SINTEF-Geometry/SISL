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
 * $Id: s9boundit.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S9BOUNDIT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9boundit(double epnt1[],double epnt2[],double epar1[],double epar2[],
	       SISLSurf *psurf1,SISLSurf *psurf2,double apar,int idir,double aepsge,
	       double gpnt1[],double gpnt2[],double gpar1[],double gpar2[],int *jstat)
#else
void s9boundit(epnt1,epnt2,epar1,epar2,psurf1,psurf2,apar,idir,aepsge,
               gpnt1,gpnt2,gpar1,gpar2,jstat)
     double epnt1[];
     double epnt2[];
     double epar1[];
     double epar2[];
     SISLSurf   *psurf1;
     SISLSurf   *psurf2;
     double apar;
     int    idir;
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
*              and a constant parameter line.
*
*
*
* INPUT      : epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in first surface
*              epnt2  - 0-2 Derivatives + normal of start point for
*                       iteration in second surface
*              epar1  - Parameter pair of start point in first surface
*              epar2  - Parameter pair of start point in second surface
*              psurf1 - Description of first surface
*              psurf2 - Description of second surface
*              apar   - Parameter value
*              idir   - Parameter direction
*                         idir = 1: The line has constant parameter value
*                                   in first direction of first surface
*                         idir = 2: The line has constant parameter value
*                                   in second direction of first surface
*                         idir = 3: The line has constant parameter value
*                                   in first direction of second surface
*                         idir = 4: The line has constant parameter value
*                                   in second direction of second surface
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
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, oct-1988
* Revised by : Tor Dokken, Si, Oslo, Norway, 24-feb-1989
*              Finds degenerate points
* Revised by : Tor Dokken, SI, Oslo, Norway, 1989-April-04
*              Correction of wrong normal vector usage
*
*********************************************************************       
*/
{
  int kcont;              /* Indicator telling if iteration is not finished */
  int kder = 2;           /* Derivative indicator                           */
  int klfu=0;             /* Pointer into knot vector                       */
  int klfv=0;             /* Pointer into knot vector                       */
  int klfs=0;             /* Pointer into knot vector                       */
  int klft=0;             /* Pointer into knot vector                       */
  int kstat;              /* Status variable                                */
  int knbit=0;            /* Counter for number of iterations               */
  int kdim = 3;           /* Set dimension to 3                             */
  int kmaxit = 100;       /* Maximal number of iterations allowed           */
  int kpos=1;             /* Position indicator ofr errors                  */
  double snorm1[3];       /* Normalvector to constant parameter line        */
  double snorm2[3];       /* Normalvector to constant parameter line        */
  double *sp,*spu,*spv,*spn; /* Pointers into gpnt1                         */
  double *sq,*sqs,*sqt,*sqn; /* Pointers into gpnt2                         */
  double ta11,ta12,ta21;  /* Variables used in equation systems             */
  double ta22,tb1,tb2;    /* Variables used in equation systems             */
  double sdiff[3];        /* Difference between two vectors                 */
  double tdum2;           /* Dummy variables                                */
  double tdum3;           /* Dummy variables                                */
  double tdist;           /* Distance betweentwo points in iteration        */
  double tdu,tdv,tds,tdt; /* Increments of parameter values                 */
  
  
  /* Copy input variables to output variables */
  
  memcopy(gpnt1,epnt1,21,DOUBLE); 
  memcopy(gpnt2,epnt2,21,DOUBLE);
  memcopy(gpar1,epar1,2,DOUBLE); 
  memcopy(gpar2,epar2,2,DOUBLE);
  
  /* At the start of the iteration the two point gpnt1 and gpnt2 might be
     very close since we in most cases start from a point on the intersection
     curve. */
  
  /* Set a number of local pointers that are used often */
  sp  = gpnt1;
  spu = gpnt1 + 3;
  spv = gpnt1 + 6;
  spn = gpnt1 + 18;
  sq  = gpnt2;
  sqs = gpnt2 + 3;
  sqt = gpnt2 + 6;
  sqn = gpnt2 + 18;
  
  kcont = 1;
  
  while (kcont)
    
    {
      if (idir==1 || idir==2)
        {
	  /* The constant parameter direction is in the first surface, intersect
	     with implicit representation of tangent plane of second surface.
	     Independent of which parameter direction is constant we want to
	     make an equation:
	     du*ta11 + dv*ta12 = tb1
	     describing the connection between du and dv. Afterwards du or dv can
	     be fixed and dv or du calculated
	   
	     Put a parametric representation of the tangent plane of surface 1 into
	     the implicit representation of the tangent plane of surface 2.
	   */ 
	  
	  ta11 = s6scpr(spu,sqn,kdim);
	  ta12 = s6scpr(spv,sqn,kdim);
	  s6diff(sq,sp,kdim,sdiff);
	  
	  tb1  = s6scpr(sdiff,sqn,kdim);
	  
	  /* Now we can branch on the constant parameter direction */
	  
	  if (idir == 1)
            {
	      /* First parameter is constant  */
	      
	      tdu = apar - gpar1[0];
	      if (DNEQUAL(ta12,DZERO))
		{
		  tdv = (tb1-tdu*ta11)/ta12;
		}
	      else
		{
		  /* spv is normal to normalvector */
		  goto war02;
		}
	      
	      gpar1[0]  = apar;
	      gpar1[1] += tdv;
            }
	  else
            {
	      /* Second parameter direction constant */
	      tdv = apar - gpar1[1];
	      if (DNEQUAL(ta11,DZERO))
		{
		  tdu = (tb1-tdv*ta12)/ta11;
		}
	      else
		{
		  /* spv is normal to normalvector */
		  goto war02;
		}
	      gpar1[0] += tdu;
	      gpar1[1]  = apar;
            }
	  
	  /* Calculate the point found in first surface */
	  
	  s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  /* If the surface has normal of zero length leave the routine */
	  
	  if (kstat == 2) goto war02;
	  
	  /* Make the difference of the found point and sq */
	  
	  s6diff(gpnt1,sq,kdim,sdiff);
	  
	  
	  /* Project the point onto surface 2 along the normal sqn */
	  
	  
	  /* Make two normals to the normal of surface two in last point */
	  
	  s6twonorm(sqn,snorm1,snorm2,&kstat);
	  if (kstat<0) goto error;
	  
	  ta11 = s6scpr(sqs,snorm1,kdim);
	  ta12 = s6scpr(sqt,snorm1,kdim);
	  ta21 = s6scpr(sqs,snorm2,kdim);
	  ta22 = s6scpr(sqt,snorm2,kdim);
	  
	  tb1  = s6scpr(sdiff,snorm1,kdim);
	  
	  tb2  = s6scpr(sdiff,snorm2,kdim);
	  
	  /*      Calculate determinant of equation system */
	  tdum2 = ta11*ta22 - ta12*ta21;
	  
	  /* If tdum2 = 0.0, then the equation system is singular, iteration not
	     possible. */
	  if (DNEQUAL(tdum2,DZERO))
            {
	      gpar2[0] += (tb1*ta22-tb2*ta12)/tdum2;
	      gpar2[1] += (ta11*tb2-ta21*tb1)/tdum2;
            }
	  
	  /* Calculate point in second surface */
	  
	  s1421(psurf2,kder,gpar2,&klfs,&klft,gpnt2,gpnt2+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  /* If the surface has normal of zero length leave the routine */
	  
	  if (kstat == 2) goto war02;
        }
      
      else
        {
	  /* idir==3 or idir==4 */
	  
	  /*  The constant parameter direction is in the second surface, intersect
	      with implicit representation of tangent plane of second surface.
	      Independent of which parameter direction is constant we want to
	      make an equation:
	      ds*ta11 + dt*ta12 = tb1
	      describing the connection between ds and dt. Afterwards ds or dt can
	      be fixed and dt or ds calculated
	   
	      Put a parametric representation of the tangent plane of surface 2 into
	      the implicit representation of the tangent plane of surface 1.
	   */ 
	  
	  ta11 = s6scpr(sqs,spn,kdim);
	  ta12 = s6scpr(sqt,spn,kdim);
	  s6diff(sp,sq,kdim,sdiff);
	  
	  tb1  = s6scpr(sdiff,spn,kdim);
	  
	  /* Now we can branch on the constant parameter direction */
	  
	  if (idir == 3)
            {
	      /* First parameter is constant  */
	      
	      tds = apar - gpar2[0];
	      if (DNEQUAL(ta12,DZERO))
		{
		  tdt = (tb1-tds*ta11)/ta12;
		}
	      else
		{
		  /* sqt is normal to normalvector */
		  goto war02;
		}
	      
	      gpar2[0]  = apar;
	      gpar2[1] += tdt;
            }     
	  else
            {
	      /* Second parameter direction constant */
	      tdt = apar - gpar2[1];
	      if (DNEQUAL(ta11,DZERO))
		{
		  tds = (tb1-tdt*ta12)/ta11;
		}
	      else
		{
		  /* spv is normal to normalvector */
		  goto war02;
		}
	      gpar2[0] += tds;
	      gpar2[1]  = apar;
	      
            }
	  
	  /* Calculate the point found in first surface */
	  
	  s1421(psurf2,kder,gpar2,&klfs,&klft,gpnt2,gpnt2+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  /* If the surface has normal of zero length leave the routine */
	  
	  if (kstat == 2) goto war02;
	  
	  /* Make the difference of the found point and sq */
	  
	  s6diff(gpnt2,sp,kdim,sdiff);
	  
	  
	  /* Project the point onto surface 2 along the normal spn */
	  
	  
	  /* Make two normals to the normal of surface one in last point */
	  
	  s6twonorm(spn,snorm1,snorm2,&kstat);
	  if (kstat<0) goto error;
	  
	  
	  /* Put a parametric representation of the tangent plane of surface 1 into
	     the implicit representation of the tangent planes of the constant
	     parameter line of surface 2 */
	  
	  ta11 = s6scpr(spu,snorm1,kdim);
	  ta12 = s6scpr(spv,snorm1,kdim);
	  ta21 = s6scpr(spu,snorm2,kdim);
	  ta22 = s6scpr(spv,snorm2,kdim);
	  
	  tb1  = s6scpr(sdiff,snorm1,kdim);
	  
	  tb2  = s6scpr(sdiff,snorm2,kdim);
	  
	  /* Calculate determinant of equation system */

	  tdum2 = ta11*ta22 - ta12*ta21;
	  
	  /* If tdum2 = 0.0, then the equation system is singular, iteration not
	     possible. */

	  if (DNEQUAL(tdum2,DZERO))
            {
	      gpar1[0] += (tb1*ta22-tb2*ta12)/tdum2;
	      gpar1[1] += (ta11*tb2-ta21*tb1)/tdum2;
            }
	  
	  /* Calculate point in first surface */
	  
	  s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  /* If the surface has normal of zero length leave the routine */
	  
	  if (kstat == 2) goto war02;
        }
      
      
      /* Make difference between the two points, 
	 and calculate length of difference */

      s6diff(gpnt1,gpnt2,kdim,sdiff);
      tdum3 = s6length(sdiff,kdim,&kstat);
      knbit = knbit + 1;
      
      if (kstat==0) 
        {
	  /* Length is zero iteration has converged   */
	  kcont = 0;
	  goto war00;
        }
      
      if (knbit<=1)
        {
	  /* First iteration intitate distance variable, if the equation
	     systems were not singular */

	  if (DEQUAL(tdum2,DZERO)) goto war02;
	  tdist = tdum3;
        }
      else
        {
	  /* More than one iteration done, stop if distance is not decreasing.
	     Then decide if we converge distance between the points is within
	     the tolerance and the last step had singular or none singular
	     equation systems. */

	  if (tdum3>=tdist)
            {
	      /* Distance is not decreasing */
	      if (tdist <= aepsge)
                {
		  /* Distance within tolerance */
		  if (DEQUAL(tdum2,DZERO))
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
	  /* Distance still decreasing */
	  
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
  
  error : *jstat = kstat;
    s6err("s9boundit",*jstat,kpos);
    goto out;
  
  out:
    return;
}
