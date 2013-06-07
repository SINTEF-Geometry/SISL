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
 * $Id: s9boundimp.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define S9BOUNDIMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9boundimp(double epnt1[],double epar1[],SISLSurf *psurf1,double eimpli[],
		int ideg,double apar,int idir,double aepsge,
		double gpnt1[],double gpar1[],int *jstat)
#else
void s9boundimp(epnt1,epar1,psurf1,eimpli,ideg,apar,idir,aepsge,
               gpnt1,gpar1,jstat)
     double epnt1[];
     double epar1[];
     SISLSurf   *psurf1;
     double eimpli[];
     int    ideg;
     double apar;
     int    idir;
     double aepsge;
     double gpnt1[];
     double gpar1[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To iterate to an intersection point between a B-spline
*              surfaces, an implicit surface and a constant parameter line.
*
*
* INPUT      : epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in B-spline surface
*              epar1  - Parameter pair of start point in B-spline surface
*              psurf1 - Description of B-spline surface              
*              eimpli - Description of implicit surface
*              ideg   - Degree of implicit surface
*                        ideg=1: Plane
*                        ideg=2; Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              apar   - The constant parameter value
*              idir   - The parameter direction of the constant parameter
*                       values:
*                         idir = 1: The line has constant parameter value
*                                   in first direction
*                         idir = 2: The line has constant parameter value
*                                   in second direction
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in B-spline surface
*              gpar1  - Parameter pair of result of iteration in B-spline
*                       surface
*              jstat  - status messages  
*                       = 2      : Iteration diverged or to many iterations
*                       = 1      : iteration converged, singular point found
*                       = 0      : ok, iteration converged
*                       < 0      : error
*
*
* METHOD     : We want to find the intersection point between the three
*              surfaces.
*
*              Ideg=1:
*                P(s,t)
*                AX = 0     The implicit represented plane given by econic
*                s  = apar  If idir = 1
*                t  = apar  If idir = 2
*
*
*              Ideg=2:
*                P(s,t)
*                XAX = 0 The implicit second degree surface
*                s  = apar  If idir = 1
*                t  = apar  If idir = 2
*
*              By making a Newton iteration on the functions we get when
*              P(s,t) is put into the implicit equations we can iterate to
*              an intersection point.
*
* USE        : The function is only working in 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 4-July-1988
* Revised by : Tor Dokken, SI, Oslo, Norway, 24-feb-1989
*              Detects now degenerate points
* Revised by : Tor Dokken, SI, Oslo; Norway, march 1989
*              Corrected sequence of statements
* Revised by : Mike Floater, SI, 1991-01
*                   Add perspective and circular silhouettes (ideg=1004,ideg=1005)
*
*********************************************************************
*/
{
  int kcont;              /* Indicator telling if iteration is not finished */
  int kder = 2;           /* Derivative indicator                           */
  int klfu=0;             /* Pointer into knot vector                       */
  int klfv=0;             /* Pointer into knot vector                       */
  int kstat;              /* Status variable                                */
  int knbit=0;            /* Counter for number of iterations               */
  int kmaxit = 100;       /* Maximal number of iterations allowed           */
  int kpos=0;             /* Position indicator ofr errors                  */
  int ksize;              /* Number of doubles for storage of derivateves
			     and normal vector */
  int ksizem3;            /* ksize - 3                                      */
  double *sp,*spu,*spv,*spn; /* Pointers into gpnt1                         */
  double ta11,ta12,tb1;   /* Variables used in equation systems             */
  double tdu,tdv;         /* Increments of u and v parameter directions     */
  double tdist;           /* Distance between two points in iteration        */
  double tcurdst;         /* Distance between points in both surfaces       */
  double sder[3];         /* Derivatives of comb. of impl. surf and par.surf*/
  double sproj[3];        /* Projection direction                           */
  
  
  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles are to be calculated */
  
  if (ideg==1003 || ideg==1004 || ideg==1005)
    {
      kder = 3;
      ksize = 33;
    }
  else
    {
      ksize = 21;
      kder =2;
    }
  ksizem3 = ksize -3;
  
  /* Copy input variables to output variables */
  
  memcopy(gpnt1,epnt1,21,DOUBLE); 
  memcopy(gpar1,epar1,2,DOUBLE); 
  
  /* At the start of the iteration the point gpnt1 is put into both implicit
     equations */
  
  /* Set a number of local pointers that are used often */
  sp  = gpnt1;
  spu = gpnt1 + 3;
  spv = gpnt1 + 6; 
  spn = gpnt1 + 18;
  
  kcont = 1;
  
  while (kcont)
    
    {
      /*  Independent of which parameter direction is constant we want to
       *   make an equation:
       *           du*ta11 + dv*ta12 = tb1
       *   describing the connection between du and dv. Afterwards du or dv can
       *   be fixed and dv or du calculated
       */
      
      /* Calculate value and derivatives of the parametric surface put into
	 the equation of the implicit surface */
      
      s1331(gpnt1,eimpli,ideg,1,sder,sproj,&kstat);
      
      ta11 = sder[1];
      ta12 = sder[2];
      tb1  = -sder[0];
      
      
      /*  Now we can branch on the constant parameter dircection */
      
      if (idir == 1)
        {
	  /* First parameter is constant  */
	  
	  tdu = apar - gpar1[0];
	  if (DNEQUAL(ta12,DZERO) )
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
      
      
      /*  Calculate value of new point */
      
      s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+18,&kstat); 
      if (kstat<0) goto error;
      
      /*  Stop iteration if degenerate point */
      if (kstat == 2) goto war02;
      
      /*  Find distance between point and point on implicit surface along sproj
       */
      tcurdst = s1309(gpnt1,sproj,eimpli,ideg,&kstat);
      if (kstat < 0) goto error;
      
      tcurdst = fabs(tcurdst);
      
      /*  tcurdst now contains the distance between the point in the parametric
	  surface and the projection along sproj of this point onto the implicit
	  surface if ideg== 1,2 or 1001. In the case ideg==1003,1004,1005 we have a
	  silhouette line and tcurdst contains the angle PI minus the angle 
	  between the view direction and the normal of the surface */
      
      
      /*  We continue iteration so long as the error tcurdst is not decreasing */
      
      knbit = knbit + 1;
      
      if (DEQUAL(tcurdst,DZERO))
        {
	  /* Length is zero iteration has converged   */
	  kcont = 0;
	  goto war00;
        }
      
      if (knbit<=1)
        {
	  /* First iteration intitate distance variable, if the equation
	     systems were not singular */
	  tdist = tcurdst;
        }
      else
        {
	  /*  More than one iteration done, stop if distance or angle is not
	      decreasing. */
	  if (tcurdst>=tdist)
            {
	      /*  Distance or angle is not decreasing */
	      if (  (ideg < 1003 && tdist <= aepsge) ||
                    (  (ideg==1003 || ideg==1004 || ideg==1005) &&
                       tdist <= ANGULAR_TOLERANCE))
                {               
		  /*  Distance within tolerance */
		  goto war00; 
		  
                }
	      else
                {
		  /* Distance is not within tolerance, divergence */
		  goto war02;
                }
            }
	  /* Distance still decreasing */
	  
	  tdist = tcurdst;
        }
      
      /*  Make sure that not to many iteration are being done */
      if (knbit > kmaxit) goto war02;
    }
  
  
  /* Iteration converged */
 war00:
  
  *jstat = 0;
  goto out;
  
  /* To many iterations or iteration diverging */
 war02: *jstat = 2;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s9boundimp",*jstat,kpos);
  goto out;
  
 out:
  return;
}
