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
 * $Id: s1773.c,v 1.3 2001-03-19 15:58:53 afr Exp $
 *
 */
#define S1773

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static
void
s1773_s9corr(double [],double,double,double,double,double,double);
static
void
s1773_s9dir(double *,double *,double *,double [],double [],double [],
		  double, int,int *);
#else
static void s1773_s9corr();
static void s1773_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1773(SISLPoint *ppoint,SISLSurf *psurf,double aepsge,
	   double estart[],double eend[],double enext[],double gpos[],int *jstat)
#else
void s1773(ppoint,psurf,aepsge,estart,eend,enext,gpos,jstat)
     SISLPoint  *ppoint;
     SISLSurf   *psurf;
     double aepsge;
     double estart[];
     double eend[];
     double enext[];
     double gpos[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a surface and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameters is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint   - The point in the closest point problem.
*              psurf    - The surface in the closest point problem.
*              aepsge   - Geometry resolution.
*              estart   - Surface parameters giving the start of the search
*                         area (umin, vmin).
*              eend     - Surface parameters giving the end of the search
*                         area (umax, vmax).
*              enext    - Surface guess parameters for the closest point
*                         iteration.
*
*
*
* OUTPUT     : gpos    - Resulting surface parameters from the iteration.
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
* WRITTEN BY : Arne Laksaa, SI, May 1989
* Revised by : Johannes Kaasa, SINTEF Oslo, August 1995.
*              Introduced a local copy of enext, to avoid changes.
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
  int kdeg;                 /* Degenaracy flag.                            */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  double *sval =SISL_NULL;       /* Value ,first and second derivatiev of surf. */ 
  double *sdiff;            /* Difference between the point and the surf.  */
  double *snorm;            /* Normal vector of the surface, dummy.        */
  double snext[2];          /* Parameter values                            */
  double guess[2];          /* Local copy of enext.                        */
  
  guess[0] = enext[0];
  guess[1] = enext[1];
  
  /* Test input.  */
  
  if (ppoint->idim != psurf->idim) goto err106;
  
  kdim = ppoint -> idim;
  
  if (kdim == 1)
    {
      s1173(ppoint,psurf,aepsge,estart,eend,guess,gpos,&kstat);
      if (kstat < 0)
        goto error;
      else
        {
	  if (DNEQUAL(gpos[0],estart[0]) &&
	      DNEQUAL(gpos[0],eend[0]) &&
	      DNEQUAL(gpos[1],estart[1]) &&
	      DNEQUAL(gpos[1],eend[1])) 
	    *jstat = (kstat==1 ? 1:3);
	  else
	    *jstat = 0;
	  goto out;
        }
    }
  
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
  /* printf("\n lin: \n %#20.20g %#20.20g",
     guess[0],guess[1]); */
  
  s1421(psurf,kder,guess,&kleft1,&kleft2,sval,snorm,&kstat);
  if (kstat < 0) goto error;
  
  /* Compute the distanse vector and value and the new step. */
  
  s1773_s9dir(&tdist,td,td+1,sdiff,ppoint->ecoef,sval,
	      aepsge,kdim,&kdeg);
  
  /* Correct if we are not inside the parameter intervall. */
  
  t1[0] = td[0];
  t1[1] = td[1];
  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 30; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = guess[0] + t1[0];
      snext[1] = guess[1] + t1[1];
      
      s1421(psurf,kder,snext,&kleft1,&kleft2,sval,snorm,&kstat);
      if (kstat < 0) goto error;
      
      /* Compute the distanse vector and value and the new step. */
      
      s1773_s9dir(&tdist,tdn,tdn+1,sdiff,ppoint->ecoef,
	    sval,aepsge,kdim,&kdeg);
      
      /* Check if the direction of the step have change. */
      
      kdir = (s6scpr(td,tdn,2) >= DZERO);     /* 0 if changed. */
      
      /* Ordinary converging. */
      
      if (tdist < tprev/(double)2 || kdir)
	{
	   guess[0] += t1[0];
	   guess[1] += t1[1];
  
	  /* printf("\n %#20.20g %#20.20g",
	     guess[0],guess[1]); */
  
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
          tprev = tdist;

	  if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	}
      
      /* Not converging, adjust and try again. */
      
      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          /* knbit--;  */
	}
      if (guess[0]==guess[0]+t1[0] &&
	  guess[1]==guess[1]+t1[1]) break;
    }
  
  /* Iteration stopped, test if point founds found is within resolution */
  
  if (tdist <= aepsge)
  {
     *jstat = 1;
     /* printf("\n SUCCESS!!"); */
     
  }
  else if(kdeg)
     *jstat = 9;
  else
     *jstat = 2;
  
  gpos[0] = guess[0];
  gpos[1] = guess[1];
  
  /* Iteration completed.  */
  
  goto out;
  
  /* Error in allocation */
  
 err101: *jstat = -101;
  s6err("s1773",*jstat,kpos);
  goto out;                  
  
  /* Error in input. Conflicting dimensions.  */
  
 err106: *jstat = -106;
  s6err("s1773",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1773",*jstat,kpos);
  goto out;                  
  
 out:    if (sval != SISL_NULL) freearray(sval);
}

#if defined(SISLNEEDPROTOTYPES)
static
void
s1773_s9corr(double gd[],double acoef1,double acoef2,
		   double astart1,double aend1,double astart2,double aend2)
#else
static void s1773_s9corr(gd,acoef1,acoef2,astart1,aend1,astart2,aend2)
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
s1773_s9dir(double *cdist,double *cdiff1,double *cdiff2,
		  double PS[],double eval1[],double eval2[],
		  double aepsge, int idim,int *jstat)
#else
static void s1773_s9dir(cdist,cdiff1,cdiff2,PS,
			eval1,eval2,aepsge,idim,jstat)
     double *cdist;
     double *cdiff1;
     double *cdiff2;
     double PS[];
     double eval1[];
     double eval2[];
     double aepsge;
     int    idim;
     int    *jstat;
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
*              eval2 - Value +1 and 2. der in surface.
*	       aepsge- Geometry tolerance.
*	       idim  - Dimension of space the surface lie in.
*
*
* OUTPUT     : PS   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection 
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection 
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*              jstat   - 0 OK, new No degeneracy.
*                        1 Degeneracy.
*
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to the point.
*	       The difference vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the point and
*	       point on the surface is "diff" and the two derivativ vectors
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
  int kstat=0;		          /* Local status variable.       */
  register double tdet;		  /* Determinant                  */
  register double t1,t2,t3,t4,t5; /* Variables in equation system */
  register double *S, *Su, *Sv;
  /* register double *Suv, *Suu, *Svv; */
                                  /* Pointers to surf values      */
  register double ref, ang;       /* Referance value, angle       */
  register double l1, l2;         /* Vector norm                  */
  register double min_ang=10e-11; /* Min angle                    */
  /* ____________________________________________________________ */
  
  /* Init */
  *jstat = 0;
  *cdiff1 = DZERO;
  *cdiff2 = DZERO;
  
  /* Set pointers */
  S   = eval2;
  Su  = S   + idim;
  Sv  = Su  + idim;
  /* Suu = Sv  + idim;
  Suv = Suu + idim;
  Svv = Suv + idim; */

  /* Degenerate if Su=0 v Sv=0 v Su||Sv */
  l1 = s6length(Su,idim,&kstat);
  l2 = s6length(Sv,idim,&kstat);
  ang = s6ang(Su,Sv,idim);
  if (min(l1,l2) < aepsge || ang < min_ang) *jstat = 1;

  /* Computing difference vector and lenght */
  s6diff(eval1,S,idim,PS);
  *cdist = s6length(PS,idim,&kstat);
  
  if (*jstat == 1)
  {
     if (l1 < aepsge)
     {
	if (l2 > aepsge)
	   /* Su = 0 */
	   *cdiff2 = s6scpr(PS,Sv,idim)/l2*l2;
     }
     else if (l2 < aepsge)
	   /* Sv = 0 */
	   *cdiff1 = s6scpr(PS,Su,idim)/(l1*l1);
     else /* Su,Sv || */
     {
	/* Best strategy? */
	*cdiff1 = s6scpr(PS,Su,idim)/(l1*l1);
      }
	
  }
  else /* *jstat == 0 */
     
  {
     
     t1 =  s6scpr(Su,Su,idim) ; /* - s6scpr(PS,Suu,idim);*/
     t2 =  s6scpr(Su,Sv,idim) ; /* - s6scpr(PS,Suv,idim);*/
     t3 =  s6scpr(Sv,Sv,idim) ; /* - s6scpr(PS,Svv,idim);*/
     t4 =  s6scpr(PS,Su,idim);
     t5 =  s6scpr(PS,Sv,idim);
     
     ref = max(fabs(t1),fabs(t2));
     ref = max(ref,fabs(t3));
     /* Computing the determinant. */
     
     tdet = t1*t3 - t2*t2;
     
     if (DEQUAL(ref+fabs(tdet),ref))
     {
	*jstat = 1;
     }
     else 
     {
	/* Using Cramer's rule to find the solution of the system. */
	
	*cdiff1 =  (t4*t3-t5*t2)/tdet;
	*cdiff2 =  (t1*t5-t2*t4)/tdet;
     }
  }
}
