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
 * $Id: crvcrvtang.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */
#define CRV_CRV_TANG

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static
void
   c_c_t_s9corr(double [],double,double,double,double,double,double);
static
void
   c_c_t_s9dir(double *,double *,double *,double [],double [],double [],
		  double [],double [], int);
#else
static void c_c_t_s9corr();
static void c_c_t_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
   crv_crv_tang(SISLCurve *pc1, SISLCurve *pc2, double aepsge,
	   double guess_par[], double iter_par[], int *jstat)
#else
     void crv_crv_tang(pc1, pc2, aepsge, guess_par, iter_par, jstat)
     SISLCurve   *pc1;
     SISLCurve   *pc2;
     double aepsge;
     double guess_par[];
     double iter_par[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration to find a tangent between the two
*              curves pc1 and pc2.
*
*
* INPUT      : pc1       - Pointer to the first curve.
*              pc2       - Pointer to the second curve.
*              aepsge    - Geometry resolution.
*              guess_par - Guess parameter values in pc1 and pc2.
*
*
*
* OUTPUT     : iter_par  - Tangential parameter values in pc1 and pc2.
*              jstat   - status messages  
*                                < 0   : error.
*
*
* METHOD     : Newton iteration on the surface ((pc2x - pc1x,pc2y - pc1y)
*              *(-(dpc1/ds)y,(dpc1/ds)x), (pc2x - pc1x,pc2y - pc1y)*
*              (-(dpc2/dt)y,(dpc2/dt)x)), towards the point (0,0).
*              Guess_par and iter_par must not be separated by a tangential
*              discontinuity.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SI, March 1992 (based on s1773).
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int kleft1=0;             /* Variables used in the evaluator.            */
  int kleft2=0;             /* Variables used in the evaluator.            */
  int kder=1;               /* Order of derivatives to be calulated        */
  int kdim;                 /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  double *sval =SISL_NULL;       /* Value ,first and second derivatiev of surf. */ 
  double *sdiff;            /* Difference between the point and the surf.  */
  double enext[2];          /* Parameter values                            */
  double snext[2];          /* Parameter values                            */
  
  double estart[2];         /* Lower borders in the parameter plane.       */
  double eend[2];           /* Higher borders in the parameter plane.      */
  double zero_pnt[2];       /* The zero point.                             */
  
  /* Test input.  */
  
  if (pc1->idim != 2 || pc2->idim != 2) goto err106;
  
  kdim = 2;
  
  /* Initiation. */
  
  estart[0] = pc1->et[pc1->ik - 1];
  estart[1] = pc2->et[pc2->ik - 1];
  eend[0] = pc1->et[pc1->in];
  eend[1] = pc2->et[pc2->in];
  zero_pnt[0] = (double)0.;
  zero_pnt[1] = (double)0.;
  enext[0] = guess_par[0];
  enext[1] = guess_par[1];
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  tdelta[0] = eend[0] - estart[0];
  tdelta[1] = eend[1] - estart[1];
  
  /* Allocate local used memory */
  
  sval = newarray(4*kdim,double);
  if (sval == SISL_NULL) goto err101;
  
  sdiff = sval + 3*kdim;
  
  /* Initiate variables.  */
  
  tprev = (double)HUGE;
  
  /* Evaluate 0-1.st derivatives of surface */
  
  eval_2_crv(pc1, pc2, kder, enext, &kleft1, &kleft2,
	     sval, &kstat);
  if (kstat < 0) goto error;
  
  /* Compute the distanse vector and value and the new step. */
  
  c_c_t_s9dir(&tdist,td,td+1,sdiff,zero_pnt,sval,sval+kdim,sval+kdim+kdim,kdim);
  
  /* Correct if we are not inside the parameter intervall. */
  
  t1[0] = td[0];
  t1[1] = td[1];
  c_c_t_s9corr(t1,enext[0],enext[1],estart[0],eend[0],estart[1],eend[1]);
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 50; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = enext[0] + t1[0];
      snext[1] = enext[1] + t1[1];
      
      eval_2_crv(pc1, pc2, kder, snext, &kleft1, &kleft2,
    	         sval, &kstat);
      if (kstat < 0) goto error;
      
      /* Compute the distanse vector and value and the new step. */
      
      c_c_t_s9dir(&tdist,tdn,tdn+1,sdiff,zero_pnt,
	    sval,sval+kdim,sval+kdim+kdim,kdim);
      
      /* Check if the direction of the step have change. */
      
      kdir = (s6scpr(td,tdn,2) >= DZERO);     /* 0 if changed. */
      
      /* Ordinary converging. */
      
      if (tdist < tprev/(double)2 || kdir)
	{
          enext[0] += t1[0];
          enext[1] += t1[1];
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  c_c_t_s9corr(t1,enext[0],enext[1],estart[0],eend[0],estart[1],eend[1]);
	  
          if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	  
          tprev = tdist;
	}
      
      /* Not converging, corrigate and try again. */
      
      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          knbit--;
	}
    }
  
  /* Iteration stopped, test if point founds found is within resolution */
  
  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;
  
  iter_par[0] = enext[0];
  iter_par[1] = enext[1];
  
  /* Iteration completed.  */
  
  goto out;
  
  /* Error in allocation */
  
 err101: *jstat = -101;
  s6err("crv_crv_tang",*jstat,kpos);
  goto out;                  
  
  /* Error in input. Conflicting dimensions.  */
  
 err106: *jstat = -106;
  s6err("crv_crv_tang",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("crv_crv_tang",*jstat,kpos);
  goto out;                  
  
 out:    if (sval != SISL_NULL) freearray(sval);
}

#if defined(SISLNEEDPROTOTYPES)
static
void
   c_c_t_s9corr(double gd[],double acoef1,double acoef2,
		   double astart1,double aend1,double astart2,double aend2)
#else
static void c_c_t_s9corr(gd,acoef1,acoef2,astart1,aend1,astart2,aend2)
     double gd[];
     double acoef1;
     double acoef2;
     double astart1;
     double aend1;
     double astart2;
     double aend2;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To be sure that we are inside the boorder of the
*              parameter plan. If we are outside clipping is used
*	       to corrigate the step value.
*
*
* INPUT      : acoef1  - Coeffisient in the first direction.
*              acoef2  - Coeffisient in the second direction.
*              astart1 - The lower boorder in first direction.
*              aend1   - The higher boorder in first direction.
*              estart2 - The lower boorder in second direction.
*              eend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gd    - Old and new step value.
*
*
* METHOD     : We are cutting a line inside a rectangle.
*	       In this case we always know that the startpoint of
*	       the line is inside the rectangel, and we may therfor
*	       use a simple kind of clipping.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{
  if (acoef1 + gd[0] < astart1)  gd[0] = astart1 - acoef1;
  else if (acoef1 + gd[0] > aend1) gd[0] = aend1 - acoef1;
  
  if (acoef2 + gd[1] < astart2)  gd[1] = astart2 - acoef2;
  else if (acoef2 + gd[1] > aend2) gd[1] = aend2 - acoef2;
}

#if defined(SISLNEEDPROTOTYPES)
static
void
   c_c_t_s9dir(double *cdist,double *cdiff1,double *cdiff2,
		  double gdiff[],double eval1[],double eval2[],
		  double eder1[],double eder2[], int idim)
#else
static void c_c_t_s9dir(cdist,cdiff1,cdiff2,gdiff,eval1,eval2,eder1,eder2,idim)
     double *cdist;
     double *cdiff1;
     double *cdiff2;
     double gdiff[];
     double eval1[];
     double eval2[];
     double eder1[];
     double eder2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance vector and value beetween
*	       a point and a point on a surface.
*	       And to compute a next step on both parameter direction
*	       This is equivalent to the nearest way to the
*	       parameter plan in the tangent plan from a point in the
*	       distance surface between a point and a surface.
*
*
* INPUT      : eval1 - Value in point.
*              eval2 - Value in point on surface.
*              eder1 - Derevative in first parameter direction.
*              eder2 - Derevative in second parameter direction.
*	       idim  - Dimension of space the surface lie in.
*
*
* OUTPUT     : gdiff   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection 
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection 
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to the point.
*	       The differens vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the point and
*	       point on the surface is "diff" and the two derevativ vectors
*	       are "der1" and "der2", and the two wanted parameter distance
*	       are "dt1" and "dt2", then we get the following system of 
*	       equations:
*		 <-dist+dt1*der1+dt2*der2,der1> = 0
*		 <-dist+dt1*der1+dt2*der2,der2> = 0
*	       This is futher:
*
*		 | <der1,der1>   <der1,der2> |  | dt1 |   | <dist,der1> |
*		 |                           |  |     | = |     	|
*		 | <der1,der2>   <der2,der2> |  | dt2 |   | <dist,der2> |
*
*	       The solution of this matrix equation is the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{                        
  int kstat;		          /* Local status variable. */
  register double tdet;		  /* Determinant */
  register double t1,t2,t3,t4,t5; /* Variables in equation system */
  
  /* Computing the different vector */
  
  s6diff(eval1,eval2,idim,gdiff);
    
  /* Computing the length of the different vector. */
  
  *cdist = s6length(gdiff,idim,&kstat);
  
  t1 =  s6scpr(eder1,eder1,idim);
  t2 =  s6scpr(eder1,eder2,idim);
  t3 =  s6scpr(eder2,eder2,idim);
  t4 =  s6scpr(gdiff,eder1,idim);
  t5 =  s6scpr(gdiff,eder2,idim);
  
  /* Computing the determinant. */
  
  tdet = t1*t3 - t2*t2;
  
  if (DEQUAL(tdet,DZERO))
    {
      *cdiff1 = DZERO;
      *cdiff2 = DZERO;
    }
  else 
    {
      /* Using Cramer's rule to find the solution of the system. */
      
      *cdiff1 =  (t4*t3-t5*t2)/tdet;
      *cdiff2 =  (t1*t5-t2*t4)/tdet;
    }
}
