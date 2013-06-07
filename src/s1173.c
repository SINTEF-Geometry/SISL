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
 * $Id: s1173.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1173

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void s1173_s9dir(double *,double *,double *,double [],double [],
			double [],double);
static void s1173_s9corr(double [],double,double,double,double,double,double);
static double s1173_s9del(double *,double *,double *,int);
#else
static void s1173_s9dir();
static void s1173_s9corr();
static double s1173_s9del();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1173(SISLPoint *ppoint, SISLSurf *psurf, double aepsge,double estart[],
     double eend[], double enext[], double gpos[],int *jstat)
#else
void s1173(ppoint,psurf,aepsge,estart,eend,enext,gpos,jstat)
     SISLPoint *ppoint;
     SISLSurf  *psurf;
     double       aepsge;
     double       estart[];
     double       eend[];
     double       enext[];
     double       gpos[];
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a onedimensional surface and a level value (point).
*              The function finds a closest point(maximum) or an
*              intersection point.
*
*
* INPUT      : ppoint   - Pointer to the point(level value) in.
*              psurface - Pointer to the surface.
*              aepsge   - Geometry resolution.
*              estart   - Start values of parameter intervalls.
*              eend     - End value of parameter intervalls.
*
*
*
* OUTPUT     : gpos    - Parameter values of the surface in intersection
*                        or closest point( maximum).
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, JUNE 1989
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int kleft1=0;             /* Variables used in the evaluator.            */
  int kleft2=0;             /* Variables used in the evaluator.            */
  int kder=2;               /* Order of derivatives to be calulated        */
  int kdim=1;               /* Dimension of space the surface lies in      */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  double *sval =SISL_NULL;       /* Value ,first and second derivatiev of surf. */ 
  double *sdiff;            /* Difference between the point and the surf.  */
  double *snorm;            /* Normal vector of the surface, dummy.        */
  double snext[2];          /* Parameter values                            */
  
  /* Test input.  */
  
  if (ppoint->idim != psurf->idim) goto err106;
  if (ppoint->idim != kdim) goto err106;
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  tdelta[0] = psurf->et1[psurf->in1] - psurf->et1[psurf->ik1 - 1];
  tdelta[1] = psurf->et2[psurf->in2] - psurf->et2[psurf->ik2 - 1];
  
  
  /* Allocate local used memory */
  
  sval = newarray(8*kdim,double);
  if (sval == SISL_NULL) goto err101;
  
  sdiff = sval + 6*kdim;
  snorm = sdiff + kdim;
  
  /* Initiate variables.  */
  
  tprev = (double)HUGE;
  
  
  /* Evaluate 0-1.st derivatives of surface */
  
  s1421(psurf,kder,enext,&kleft1,&kleft2,sval,snorm,&kstat);
  if (kstat < 0) goto error;
  
  /* Compute the distanse vector and value and the new step. */
  
  s1173_s9dir(&tdist,td,td+1,sdiff,ppoint->ecoef,sval,aepsge);
  
  
  /* Correct if we are not inside the parameter intervall. */
  
  
  t1[0] = td[0];
  t1[1] = td[1];
  s1173_s9corr(t1,enext[0],enext[1],estart[0],eend[0],estart[1],eend[1]);
  
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 50; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = enext[0] + t1[0];
      snext[1] = enext[1] + t1[1];
      
      s1421(psurf,kder,snext,&kleft1,&kleft2,sval,snorm,&kstat);
      if (kstat < 0) goto error;
      
      
      /* Compute the distanse vector and value and the new step. */
      
      s1173_s9dir(&tdist,tdn,tdn+1,sdiff,ppoint->ecoef,sval,aepsge);
      
      
      /* Check if the direction of the step have change. */
      
      kdir = (s6scpr(td,tdn,2) >= DZERO);     /* 0 if changed. */
      
      
      /* Ordinary converging. */
      
      if (tdist <= tprev || kdir)
	{
          enext[0] += t1[0];
          enext[1] += t1[1];
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  s1173_s9corr(t1,enext[0],enext[1],estart[0],eend[0],estart[1],eend[1]);
	  
	  
          if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	  
          tprev = tdist;
	}
      
      /* Not converging, corrigate and try again. */
      
      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
	}
    }
  
  /* Iteration stopped, test if point is within resolution */
  
  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;
  
  /* Test if the iteration is close to a knot */
  if (DEQUAL(enext[0],psurf->et1[kleft1]))
    gpos[0] = psurf->et1[kleft1];
  else if (DEQUAL(enext[0],psurf->et1[kleft1+1]))
    gpos[0] = psurf->et1[kleft1+1];
  else
    gpos[0] = enext[0];
  
  if (DEQUAL(enext[1],psurf->et2[kleft2]))
    gpos[1] = psurf->et2[kleft2];
  else if (DEQUAL(enext[1],psurf->et2[kleft2+1]))
    gpos[1] = psurf->et2[kleft2+1];
  else
    gpos[1] = enext[1];
  
  
  /* Iteration completed.  */
  
  
  goto out;
  
  
  /* Error in allocation */
  
 err101: *jstat = -101;
  s6err("s1173",*jstat,kpos);
  goto out;                  
  
  /* Error in input. Conflicting dimensions.  */
  
 err106: *jstat = -106;
  s6err("s1173",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1173",*jstat,kpos);
  goto out;                  
  
 out:    if (sval != SISL_NULL) freearray(sval);
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1173_s9corr(double gd[], double acoef1,double acoef2,double astart1,
		   double aend1,double astart2, double aend2)
#else
static void s1173_s9corr(gd,acoef1,acoef2,astart1,aend1,astart2,aend2)
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
static void
s1173_s9dir(double *cdist, double *cdiff1, double *cdiff2,
		  double gdiff[], double evalp[], double evals[], double aepsge)
#else
static void s1173_s9dir(cdist,cdiff1,cdiff2,gdiff,evalp,evals,aepsge)
     double *cdist;
     double *cdiff1;
     double *cdiff2;
     double gdiff[];
     double evalp[];
     double evals[];
     double aepsge;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance vector and value beetween
*	       a point and a point on a surface.
*	       And to compute a next step on both parameter direction
*
*
* INPUT      : evalp - Value in point.
*              evals - Value and derivatives in point on surface.
*              aepsge  - The geometry resolution.
*
* OUTPUT     : gdiff   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection 
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection 
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*
*
* METHOD     : This is a one dimensional case. We want to minimize
*
*                   min <(point - S(x,y),point - S(x,y)>
*                    x,y
*
*              This is equivalent to     
*
*                    x,y : (point-S(x,y))Sx(x,y) = 0      
*                          (point-S(x,y))Sy(x,y) = 0      
*
*              Using Taylor we get:
*                  
*                    x,y : (point-S-DxSx-DySy)(Sx+DxSxx+DySxy) =0
*                          (point-S-DxSx-DySy)(Sy+DySyy+DxSxy) =0
*
*              If we neglect the Dx**2,Dy**2 and Dx*Dy parts, we get:
*
*   x,y: [(point-S)Sxx - Sx**2]Dx + [(point-S)Sxy - Sx*Sy]Dy = - (point-S)Sx 
*        [(point-S)Syy - Sy**2]Dy + [(point-S)Sxy - Sx*Sy]Dx = - (point-S)Sy 
*
*	       The solution of this matrix equation is the
*	       following function.
*              If, however, the matrix is singular, 
*              we estimate Dx and Dy separately
*              by Newton iteration in one variable. 
*              We combine these two increments
*              by convexcombination to get the minimum distance to move(!=0)
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, June 1989
*
*********************************************************************
*/                       
{                        
  int kstat=0;		      /* Local status variable.                    */
  double tdiv;		      /* Determinant                               */
  double ta11,ta12,ta21,ta22; /* The matrix                  		   */
  double tmax;                /* The largest value in matrix               */
  double tb1,tb2;             /* The right hand side.                      */
  double tval,tderx,tderxx;   /* Function and deriv. 
				 values in one-dimentional case */
  double tdery,tderyy;
  double tderxy;
  double tdeltax,tdeltay;   /* Locals for the step value to be determined. */
  double ttemp;             /* Temporary value. */
  
  if (aepsge < 0) kstat=1;
  
  /* Computing the different vector */
  s6diff(evalp,evals,1,gdiff);
  
  /* Computing the length of the different vector. */
  *cdist = s6length(gdiff,1,&kstat);
  
  /* Init */
  tval   = evals[0];
  tderx  = evals[1];
  tdery  = evals[2];
  tderxx = evals[3];
  tderxy = evals[4];
  tderyy = evals[5];
  tdeltax = DZERO;
  tdeltay = DZERO;
  *cdiff1  = DZERO;
  *cdiff2  = DZERO;
  
  
  /* Building the matrix. */
  
  ta11 = (gdiff[0]*tderxx - tderx*tderx);
  ta12 = (gdiff[0]*tderxy - tderx*tdery);
  ta21 = (gdiff[0]*tderxy - tderx*tdery);
  ta22 = (gdiff[0]*tderyy - tdery*tdery);
  tb1  = -gdiff[0]*tderx;
  tb2  = -gdiff[0]*tdery;
  
  if (DEQUAL(tb1,DZERO) && DEQUAL(tb2,DZERO))
    {
      /* Finished, we have found a max. */
    }
  else
    {
      tdiv    = ta11*ta22 - ta21*ta12;
      tmax = max(fabs(ta11),max(fabs(ta12),max(fabs(ta21),fabs(ta22))));
      
      if (fabs(tdiv) > tmax*REL_COMP_RES)
	{
	  /* The matrix is ok, solve the system using Cramers rule. */
	  tdeltax = tb1*ta22 - tb2*ta12;    
	  tdeltay = ta11*tb2 - ta21*tb1;
	  tdeltax /= tdiv;
	  tdeltay /= tdiv;
	}
      else
	{
	  /* The matrix is nearly singular, 
	     use Newton on each parameter direction*/
	  tdeltax = s1173_s9del(gdiff,&tderx,&tderxx,1);
	  tdeltay = s1173_s9del(gdiff,&tdery,&tderyy,1);
	  
	  
	  if (fabs(tdeltax) < REL_COMP_RES || fabs(tdeltay) < REL_COMP_RES )
	    /* If one is very small, we use them as they are. */
	    ;
	  else
	    {
	      /* Use the shortest step; min (1-k)Dx + kDy */
	      ttemp   = tdeltay*tdeltax/(tdeltax*tdeltax + tdeltay*tdeltay);
	      tdeltax = tdeltay*ttemp;
	      tdeltay = tdeltax*ttemp;
	      
	    }
	  
	} 
    }  
  
  *cdiff1  = tdeltax;
  *cdiff2  = tdeltay;
  
}

#if defined(SISLNEEDPROTOTYPES)
static double
s1173_s9del(double *eco, double *eco1, double *eco2, int idim)
#else
static double s1173_s9del(eco,eco1,eco2,idim)
     double *eco;
     double *eco1;
     double *eco2;
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance on the parameter line to a point
*            on the curve on which the tangent is orthogonal to the
*            difference vector from this point on the curve to the
*            point in the space.
*
*
* INPUT      : eco   - The differens vector beetween the point and the
*                      current posision on the curve.
*              eco1  - The first derevative vector in the  current posision
*                      on the curve.
*              eco2  - The second derevative vector in the  current posision
*                      on the curve.
*              idim  - Dimension of space the vectors lie in.
*
*
* OUTPUT     : s1173_s9del - The computed parameter distance.
*
*
* METHOD     : We have to find the parameter distance "dt" from
*              the equation:
*                <ecoef-dt*ecoef1,ecoef1+dt*ecoef2> = 0.
*              This may be written:
*                  <ecoef,ecoef1> + <ecoef,ecoef2>*dt
*                - <ecoef1,ecoef1>*dt + <ecoef1,ecoef2>*dt*dt = 0
*              The following function is the solution of this second
*              degree equation. We are always using the solution
*              with the smallest absolute value.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Mar 1989
*
*********************************************************************
*/                       
{
  double t1,t2,t3,t4,t5,t6;   /* Constants in equation.                 */
  
  t1 =  s6scpr(eco,eco1,idim);
  t3 =  s6scpr(eco1,eco1,idim);
  t2 =  t3 - s6scpr(eco,eco2,idim);
  t4 =  -(double)2 * s6scpr(eco1,eco2,idim);
  
  
  
  if (DEQUAL(t4,DZERO))    /* The second degree part is degenerated. */
    {
      if (DEQUAL(t2,DZERO)) 
	{
          if (DEQUAL(t3,DZERO))            return DZERO;
          else                             return (t1/t3);
	}
      else                                  return (t1/t2);
    }
  else                /* An ordinary second degree equation.    */
    {
      t5 = t2*t2 - (double)2*t4*t1;
      if (t5 < DZERO)                       return (t1/t3);
      else
	{
          t6 = sqrt(t5);
          t5 = (t2 + t6)/t4;
          t6 = (t2 - t6)/t4;
	  t1 *= t3;
	  
	  
          /* We have two solutions and we want to use the one
	     with the same sign as we get while using an other
	     metode t1/t3. If both solutions have the same
	     sign we use the one with smallest value. */
	  
          if (t1 < DZERO)
	    {
	      if (t5 <= DZERO && t6 <= DZERO)
		{
		  if (t5 > t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 <= DZERO)        return t5;
	      else if (t6 <= DZERO)        return t6;
              else                         return min(t5,t6);
	    }
	  else if (t1 > DZERO)
	    {
	      if (t5 >= DZERO && t6 >= DZERO)
		{
		  if (t5 < t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 >= DZERO)        return t5;
	      else if (t6 >= DZERO)        return t6;
              else                         return max(t5,t6);
	    }
	  else                             return min(fabs(t5),fabs(t6));
	}
    }
}












