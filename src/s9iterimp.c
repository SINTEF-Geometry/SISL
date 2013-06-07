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
 * $Id: s9iterimp.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */
#define S9ITERIMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9iterimp(double epoint[],double epnt1[],double epar1[],SISLSurf *psurf1,
	       double eimpli[],int ideg,double astep,double aepsge,
	       double gpnt1[],double gpar1[],int *jstat)
#else
void s9iterimp(epoint,epnt1,epar1,psurf1,eimpli,ideg,astep,aepsge,
               gpnt1,gpar1,jstat)
     double epoint[];
     double epnt1[];
     double epar1[];
     SISLSurf   *psurf1;
     double eimpli[];
     int    ideg;
     double astep;
     double aepsge;
     double gpnt1[];
     double gpar1[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To iterate to an intersection point between a B-spline
*              surfaces, an implicit surface and a plane.
*
*
* INPUT      : epoint - Array containing parts of plane description.
*                       epoint[0:2] contains a position value.
*                       epoint[3:5] contains the normal to the plane
*                       A point in the plane is defined by
*                       epoint[0:2] + astep*epoint[3:5]
*              epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in B-spline surface
*              epar1  - Parameter pair of start point in B-spline surface
*              psurf1 - Description of B-spline surface              
*              eimpli - Description of implicit surface
*              ideg   - Degree of implicit surface
*                        ideg=1:    Plane              
*                        ideg=2;    Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              astep  - Step length
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
*                AX = 0  The implicit represented plane given by econic
*                BX = 0  The implicit represented plane giving the step
*
*
*              Ideg=2:
*                P(s,t)
*                XAX = 0 The implicit second degree surface
*                BX  = 0 The implicit represented plane giving the step
*
*
*              Ideg=1001; Torus surface
*                P(s,t)
*                Torus described by center, normal, big and small radius
*
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
* Revised by : Tor Dokken, SI, Oslo, Norway, 24-Feb-1989
*              Detects degenerate points
* Revised by : Tor Dokken, SI, Oslo, Norway, March 1989
*              Test for almost zero determinant introduced.
* Revised by : Mike Floater, SI, 1991-01
*                   Add perspective and circular silhouettes (ideg=1004,ideg=1005)
*
*********************************************************************
*/
{
  int ki;                 /* Variable used in loop                          */
  int kcont;              /* Indicator telling if iteration is not finished */
  int kder = 1;           /* Derivative indicator                           */
  int klfu=0;             /* Pointer into knot vector                       */
  int klfv=0;             /* Pointer into knot vector                       */
  int kstat;              /* Status variable                                */
  int knbit;              /* Counter for number of iterations               */
  int kdim = 3;           /* Set dimension to 3                             */
  int kmaxit = 100;       /* Maximal number of iterations allowed           */
  int kpos=1;             /* Position indicator ofr errors                  */
  int ksing;              /* Singularity indicator                          */
  int ksize;              /* Number of doubles for storage of derivateves
			     and normal vector */
  int ksizem3;            /* ksize - 3                                      */
  double spoint[3];       /* SISLPoint in intersection plane                    */
  double *snorm;          /* Pointer to normal vector of intersection plane */
  double sbinorm[3];      /* Vector normal tu curve tangent                 */
  double *sp,*spu,*spv,*spn; /* Pointers into gpnt1                         */
  double sprev[3];        /* Coordinates of previous point in iteration     */
  double ta11,ta12,ta21;  /* Variables used in equation systems             */
  double ta22,tb1,tb2;    /* Variables used in equation systems             */
  double sdiff[3];        /* Difference between two vectors                 */
  double tdum,tdum1;      /* Dummy variables                                */
  double tdist;           /* Error so fare in iteration                     */
  double tcurdst;         /* Error at current step in the iteration         */
  double sder[3];         /* Derivatives of comb. of impl. surf and par.surf*/
  double sproj[3];        /* Projection direction                           */
  double tlnorm;          /* Length of normal vector of step plane          */
  double tdiststep;       /* Distance from step plane                       */
  double titer;           /* Iteration criteria                             */
  
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
  
  /* Make description of intersection plane */
  
  tlnorm = s6length(epoint+3,3,&kstat);
  if (kstat<0) goto error;
  if (DEQUAL(tlnorm,DZERO)) tlnorm = (double)1.0;
  
  for (ki=0;ki<3;ki++)
    {
      spoint[ki] = epoint[ki] + astep*epoint[ki+3]/tlnorm;
    }
  
  snorm = epoint + 3;
  
  /* Copy input variables to output variables */
  
  memcopy(gpnt1,epnt1,ksize,DOUBLE); 
  memcopy(gpar1,epar1,2,DOUBLE); 
  
  /* At the start of the iteration the point gpnt1 is put into both implicit
     equations */
  
  /* Set a number of local pointers that are used often */
  sp  = gpnt1;
  spu = gpnt1 + 3;
  spv = gpnt1 + 6;
  spn = gpnt1 + ksizem3;
  
  kcont = 1;
  knbit = 0;
  
  while (kcont)
    
    {             
      
      /* For all degrees we have to put the B-spline surface into the plane
	 determining the step length before we branch degree on the degree
	 of the implicit surface. */
      
      ta21 = s6scpr(spu,snorm,kdim);
      ta22 = s6scpr(spv,snorm,kdim);
      s6diff(spoint,sp,kdim,sdiff);
      tb2  = s6scpr(sdiff,snorm,kdim);
      tdum = max(fabs(ta21),fabs(ta22));
      tdum = max(tdum,fabs(tb2));
      if (DEQUAL(tdum,DZERO)) tdum = (double)1.0;
      ta21 /= tdum;
      ta22 /= tdum;
      tb2  /= tdum;
      
      
      /* Calculate value and derivatives of the parametric surface put into
	 the equation of the implicit surface */
      
      s1331(gpnt1,eimpli,ideg,1,sder,sproj,&kstat);
      
      ta11 = sder[1];
      ta12 = sder[2];
      tb1  = -sder[0];
      
      tdum = max(fabs(ta11),fabs(ta12));
      tdum = max(tdum,fabs(tb1));
      if (DEQUAL(tdum,DZERO)) tdum = (double)1.0;
      ta11 /= tdum;
      ta12 /= tdum;
      tb1  /= tdum;
      
      
      /* Calculate determinant of equation system */
      
      tdum1 = ta11*ta22 - ta12*ta21;
      tdum  = MAX(fabs(ta11),fabs(ta22));
      tdum  = MAX(fabs(ta12),tdum);
      tdum  = MAX(fabs(ta21),tdum);
      
      if (DEQUAL((tdum+tdum1),tdum)) tdum1 =DZERO;
      
      
      /* If tdum1 = 0.0, then the equation system is singular, try an
	 alternative setup of the equation system */
      
      if (tdum1 == DZERO && ideg < 1003)
        {
	  s6crss(sproj,snorm,sbinorm);
	  ta11 = s6scpr(spu,sbinorm,kdim);
	  ta12 = s6scpr(spv,sbinorm,kdim);
	  tb1  = s6scpr(sdiff,sbinorm,kdim);
	  
	  tdum = max(fabs(ta11),fabs(ta12));
	  tdum = max(tdum,fabs(tb1));
	  if (DEQUAL(tdum,DZERO)) tdum = (double)1.0;
	  ta11 /= tdum;
	  ta12 /= tdum;
	  tb1  /= tdum;
	  
	  /* Calculate determinant of equation system */
	  
	  tdum1 = ta11*ta22 - ta12*ta21;
	  tdum  = MAX(fabs(ta11),fabs(ta22));
	  tdum  = MAX(fabs(ta12),tdum);
	  tdum  = MAX(fabs(ta21),tdum);
	  
	  if (DEQUAL((tdum+tdum1),tdum)) tdum1 =DZERO;
        }
      
      if (DNEQUAL(tdum1,DZERO))
        {
	  gpar1[0] += (tb1*ta22-tb2*ta12)/tdum1;
	  gpar1[1] += (ta11*tb2-ta21*tb1)/tdum1;
        }
      
      /* Calculate value of new points */
      
      s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+ksizem3,&kstat); 
      if (kstat<0) goto error;
      
      /* If the surface normal has zero length no use in continuing */
      
      if (kstat == 2) goto war02;
      
      
      tcurdst = s1309(gpnt1,sproj,eimpli,ideg,&kstat);
      if (kstat < 0) goto error;
      if (kstat ==2) goto war02;
      
      tcurdst = fabs(tcurdst);
      
      /* Calculate distance from step plane */
      
      s6diff(spoint,gpnt1,kdim,sdiff);
      tdiststep  = fabs(s6scpr(sdiff,snorm,kdim)/tlnorm);
      
      
      /* tcurdst now contains the distance between the point in the parametric
	 surface and the projection along sproj of this point onto the implicit
	 surface if ideg== 1,2 or 1001. In the case ideg==1003,1004,1005 we have a
	 silhouette line and tcurdst contains the angle PI minus the angle 
	 between the view direction and the normal of the surface */
      
      /* We continue iteration so long as the error titer is not decreasing */
      
      
      if (ideg==1003 || ideg==1004 || ideg==1005)
        {
	  /* tcurdst contains an angle and is compared with ANGULAR_TOLERANCE,
	     while tdiststep contains a distance and is compared with aepsge.
	     To make a measure of these that is consistent they have to have
	     the same unit measure thus we make.   */
	  
	  titer = tcurdst*aepsge + tdiststep*ANGULAR_TOLERANCE;
        }
      else
        titer = tcurdst + tdiststep;
      
      if (DEQUAL(tcurdst,DZERO) && DEQUAL(tdiststep,DZERO))
        {
	  /* Length is zero iteration has converged   */
	  kcont = 0;
        }
      
      if (knbit==0)
        {
	  /* First iteration intitate distance variable, if the equation
	     systems were not singular */

	  if (DEQUAL(tdum1,DZERO)) goto war02;
	  tdist = titer;
	  knbit = 1;
        }
      else
        {
	  /* More than one iteration done, stop if distance is not decreasing.
	     Then decide if we converge distance between the points is within
	     the tolerance and the last step had singular or none singular
	     equation systems. */

	  knbit = knbit + 1;
	  if (titer>=tdist)
            {
	      /* Distance is not decreasing */

	      if (fabs(s6dist(sprev,gpnt1,kdim)) <= aepsge)
                {               
		  /* Distance within tolerance */
		  
		  /* Check if singularity reached. 
		     This is the case if tdum1=0.0
		     or if the relative distance between
		     the input point and
		     the output point is within the 
		     relative computer resolution */

		  ksing = 1;
		  for (ki=0 ; ki<3; ki++)
                    {
		      tdum = MAX(fabs(epnt1[ki]),fabs(gpnt1[ki]));
		      if (DEQUAL(tdum,DZERO)) tdum = (double)1.0;
		      if(fabs(epnt1[ki]-gpnt1[ki])/tdum > REL_COMP_RES)
                        ksing = 0;
                    }
		  if (DEQUAL(tdum1,DZERO) || ksing == 1)
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
	  
	  tdist = titer;
        }
      
      /*  Make sure that not to many iteration are being done */
      if (knbit > kmaxit) goto war02;
      /*  Remember this point */
      memcopy(sprev,gpnt1,3,DOUBLE);
    }
  
  /* Iteration converged, calculate also second derivatives */
 war00:
  
  kder = 2;
  s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+ksizem3,&kstat); 
  if (kstat<0) goto error;
  
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
  s6err("s9iterimp",*jstat,kpos);
  goto out;
  
 out:
  return;
}
