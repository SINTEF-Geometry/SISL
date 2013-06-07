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
 * $Id: shsing.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */

#define SHSING

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void shsing_s9corr(double [], double [], double[]);
static void shsing_s9dir(double [], double [], double[]);
#else
static void shsing_s9corr();
static void shsing_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  shsing(SISLSurf *psurf1,SISLSurf *psurf2,double limit[],
	 double enext[], double gpos[],int *jstat)
#else
void shsing(psurf1,psurf2,limit,enext,gpos,jstat)
     SISLSurf         *psurf1;
     SISLSurf         *psurf2;
     double       limit[];
     double       enext[];
     double       gpos[];
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To find one point in each surface where the two normals
*              are parallel to each other and to the difference vector 
*              between the two points.
*
* INPUT      : psurf1 - Pointer to first surface
*              psurf2 - Pointer to second surface
*              limit  - Parameter boarders of both surfaces.
*                       limit[0] - limit[1] Parameter interval first
*                                  surface first direction.
*                       limit[2] - limit[3] Parameter interval first 
*                                  surface second direction.
*                       limit[4] - limit[5] Parameter interval second 
*                                  surface first direction.
*                       limit[6] - limit[7] Parameter interval second 
*                                  surface second direction.
*              enext    - Parameter start value for iteration(4 values).
*
*
* OUTPUT     : gpos    - Parameter values of the found singularity.(4 values)
*              jstat   - status messages  
*                                = 1   : Extremum found.
*                                = 0   : Extremum NOT found.
*                                < 0   : error.
*
*
* METHOD     :  - Start with a guess value (u,v) in domain of 
*                 surface 1 (S(u,v))
*           (a) - Find domain value (r,t) of closest point (to S(u,v) 
*                 in surface 2 (Q(r,t))
*               - If vf1(u,v) = <Su,Normal(Q> and vf2(u,v)= <Sv,Normal(Q> 
*                 is small enough  stop.  (<,> means scalar prod.)
*               - Find du and dv by taylorizing vf1 and vf2.
*                 This include finding the derivatives of the closest point 
*                 function (r(u,v),t(u,v)) with respect to u and v. 
*                 (called h(u,v) in article, see coments in shsing_s9dir)
*               - u:= u+du v:= v+dv, goto (a)
*
*
* REFERENCES : Solutions of tangential surface and curve intersections.
*              R P Markot and R L Magedson
*              Computer-Aided Design; vol. 21, no 7 sept. 1989, page 421-429
*
*
* WRITTEN BY : Ulf J. Krystad, SI, JANUARY 1991
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int ki;                   /* Loop control                                */
  int kleftt=0;             /* Variables used in the evaluator.            */
  int klefts=0;             /* Variables used in the evaluator.            */
  int kleftu=0;             /* Variables used in the evaluator.            */
  int kleftv=0;             /* Variables used in the evaluator.            */
  int kder=2;               /* Order of derivatives to be calulated        */
  int kdim=3;               /* Dimension of space the surface lies in      */
  int knbit;                /* Number of iterations                        */
  double tdelta[4];         /* Length of parameter intervals.              */
  double tdist;             /* The current norm of the cross product       */
                            /* between the two normals                     */
  double tprev;             /* The current norm of the cross product       */
                            /* between the two normals                     */
  double td[4],t1[4],tdn[4];/* Distances between old and new parameter     */
			    /* value in the four parameter directions.     */
  double sval1[21];         /* Value ,first and second derivatiev of surf. */ 
  double *snorm1=sval1+18;  /* Normal vector of the surface                */
  double sval2[21];         /* Value ,first and second derivatiev of surf. */ 
  double *snorm2=sval2+18;  /* Normal vector of the surface                */
  double snext[4];          /* Parameter values                            */
  double temp[3];           /* Temp vector storing cross products.         */
  double start[2];          /* Parameters limit of second surface, used in */
                            /* call to closest point                       */
  double end[2];            /* Parameters limit of second surface, used in */
                            /* call to closest point                       */
  double guess[2];          /* Start point for closest point iteration     */
  double tol = (double)10000.0*REL_COMP_RES; /* equality tol. in par.space */
  SISLPoint *ppoint=SISL_NULL;   /* Contains the current position in first      */ 
                            /* surface used in closest point iteration     */
  int max_iter=20;          /* Maximal number of iteration allowed         */

  /* --------------------------------------------------------------------- */
  
  /* Test input.  */
  if (psurf1->idim != kdim) goto err106;
  if (psurf2->idim != kdim) goto err106;
  
  /* Fetch referance numbers from the serach intervals for the surfaces.  */
  tdelta[0] = limit[1] - limit[0];
  tdelta[1] = limit[3] - limit[2];
  tdelta[2] = limit[5] - limit[4];
  tdelta[3] = limit[7] - limit[6];

  /* Set limit values, used in closest point iteration */
  start[0] = limit[4];
  start[1] = limit[6];
  end[0]   = limit[5];
  end[1]   = limit[7];
  
  /* Create point, used in closest point iteration */
  ppoint = newPoint(sval1,3,0);
  
  /* Collapsed ? */
  for (ki=0;ki<4;ki++) if (tdelta[ki] < tol) goto errsmall;  
  
  /* Initiate output variables.  */
  for (ki=0;ki<4;ki++)     gpos[ki] = enext[ki];

  /* Evaluate 0.-2. derivatives of first surface */
  s1421(psurf1,kder,gpos,&kleftt,&klefts,sval1,snorm1,&kstat);
  if (kstat < 0) goto error;

  /* Get closest point in second surface. */
  guess[0] = gpos[2];
  guess[1] = gpos[3];
  s1773(ppoint,psurf2,REL_COMP_RES,start,end,guess,gpos+2,&kstat);
  if (kstat < 0) goto error;
  
  /* Evaluate 0.-2. derivatives of second surface */
  s1421(psurf2,kder,gpos+2,&kleftu,&kleftv,sval2,snorm2,&kstat);
  if (kstat < 0) goto error;

  /* Get length of normal cross product */
  s6crss(snorm1,snorm2,temp);
  tprev = s6length(temp,kdim,&kstat);
  
  /* Compute the Newton stepdistance vector in first surface. */
  shsing_s9dir(td,sval1,sval2);
  
  /* Adjust if we are not inside the parameter intervall. */
  for (ki=0;ki<4;ki++)    t1[ki] = td[ki];

  shsing_s9corr(t1,gpos,limit);
  
  /* Iteratation loop.  */
  
  for (knbit = 0; knbit < max_iter; knbit++)
    {
      
      for (ki=0;ki<2;ki++)    snext[ki] = gpos[ki] + t1[ki];
   
      /* Evaluate 0.-2. derivatives of first surface */
      s1421(psurf1,kder,snext,&kleftt,&klefts,sval1,snorm1,&kstat);
      if (kstat < 0) goto error;
      
      /* Get closest point in second surface. */
      guess[0] = gpos[2];
      guess[1] = gpos[3];
      s1773(ppoint,psurf2,REL_COMP_RES,start,end,guess,snext+2,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate 0.-2. derivatives of second surface */
      s1421(psurf2,kder,snext+2,&kleftu,&kleftv,sval2,snorm2,&kstat);
      if (kstat < 0) goto error;

      /* Get length of normal cross product */
      s6crss(snorm1,snorm2,temp);
      tdist = s6length(temp,kdim,&kstat);
  
      /* Compute the Newton stepdistance vector. */
      shsing_s9dir(tdn,sval1,sval2);
      
      if (tdist <= tprev)
	{
	  /* Ordinary converging. */
	  
	  for (ki=0;ki<4;ki++)
	    {
	      gpos[ki] = snext[ki];
	      td[ki] = t1[ki] = tdn[ki];
	    }
	  
	  /* Adjust if we are not inside the parameter intervall. */
	  shsing_s9corr(t1,gpos,limit);
		  
          if ((fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES) &&
	      (fabs(t1[2]/tdelta[2]) <= REL_COMP_RES) &&	   
	      (fabs(t1[3]/tdelta[3]) <= REL_COMP_RES))
	      {
		for (ki=0;ki<2;ki++) gpos[ki] += t1[ki];
		/* Evaluate 0.-2. derivatives of first surface */
		s1421(psurf1,kder,gpos,&kleftt,&klefts,sval1,snorm1,&kstat);
		if (kstat < 0) goto error;
		
		/* Get closest point in second surface. */
		guess[0] = gpos[2];
		guess[1] = gpos[3];
		s1773(ppoint,psurf2,REL_COMP_RES,start,end,guess,gpos+2,&kstat);
		if (kstat < 0) goto error;
		break;
	      }
	  tprev = tdist;
	}
      
      else
	{
	  /* Not converging, half step length try again. */
	  
	  for (ki=0;ki<4;ki++) t1[ki] /= (double)2;
	}
    }
  
  /* Iteration stopped, test if point is extremum */
  /* Unsure about what i right here , angle between normals and difference vector ?? */
  if (tdist <= tol)
    *jstat = 1;
  else
    *jstat = 0;

 
  /* Test if the iteration is close to a knot */
  if (fabs(gpos[0] - psurf1->et1[kleftt])/tdelta[0] < tol)
    gpos[0] = psurf1->et1[kleftt];
  else if (fabs(gpos[0] - psurf1->et1[kleftt+1])/tdelta[0] < tol)
    gpos[0] = psurf1->et1[kleftt+1];
  
  if (fabs(gpos[1] - psurf1->et2[klefts])/tdelta[1] < tol)
    gpos[1] = psurf1->et2[klefts];
  else if (fabs(gpos[1] - psurf1->et2[klefts+1])/tdelta[1] < tol)
    gpos[1] = psurf1->et2[klefts+1];

  if (fabs(gpos[2] - psurf2->et1[kleftu])/tdelta[2] < tol)
    gpos[2] = psurf2->et1[kleftu];
  else if (fabs(gpos[2] - psurf2->et1[kleftu+1])/tdelta[2] < tol)
    gpos[2] = psurf2->et1[kleftu+1];
  
  if (fabs(gpos[3] - psurf2->et2[kleftv])/tdelta[3] < tol)
    gpos[3] = psurf2->et2[kleftv];
  else if (fabs(gpos[3] - psurf2->et2[kleftv+1])/tdelta[3] < tol)
    gpos[3] = psurf2->et2[kleftv+1];
  
  /* Iteration completed.  */
  goto out;
  
  /* --------------------------------------------------------------------- */ 
  /* Error in input. Dimension not equal to 3 */
  err106: *jstat = -106;
  s6err("shsing",*jstat,kpos);
  goto out;                  

  /* Error in input. One parameter interval colapsed. */
  errsmall: *jstat = -200;
  s6err("shsing",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("shsing",*jstat,kpos);
  goto out;                  
  
 out:if(ppoint) freePoint(ppoint);
}

#if defined(SISLNEEDPROTOTYPES)
static void
  shsing_s9corr(double gd[], double coef[],double limit[])
#else
static void shsing_s9corr(gd,coef,limit)
     double gd[];
     double coef[];
     double limit[];
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To be sure that we are inside the boarder of the
*              parameter plan. If we are outside clipping is used
*	       to adjust the step value.
*
*
* INPUT      : coef    - Current position.
*              limit   - Parameter boarders of both surfaces.
*
*
*
* INPUT/OUTPUT : gd    - Proposed delta values.
*
*
* METHOD     : Cutting the line towards the parameter box.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, NOVEMBER 1990
*
*********************************************************************
*/                       
{
  int ki;

  for (ki=0;ki<4;ki++)
    if (coef[ki] + gd[ki] < limit[2*ki])        gd[ki] = limit[2*ki]    - coef[ki];
    else if (coef[ki] + gd[ki] > limit[2*ki+1]) gd[ki] = limit[2*ki +1] - coef[ki];
  

}

#if defined(SISLNEEDPROTOTYPES)
static void
  shsing_s9dir(double cdiff[],double evals[],double evalq[])
#else
static void shsing_s9dir(cdiff,evals,evalq)
     double cdiff[];
     double evals[];
     double evalq[];
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the increments in the first surface domain.
*              
*
*
* INPUT      : evals - Value and derivatives  on first surface.
*              evalq - Value and derivatives  on second surface.
*              
*
* OUTPUT     : cdiff1  - Parameter increments in two directions.
*            
*
*
* METHOD     : See comments in main header.
*              The only thing missing in the article is the derivation of the
*              function h(u,v). Calling this function (r(u,v), t(u,v))
*              we know that
*              For all (u,v) (<,> meaning scalar product)
*              <S(u,v)-Q(r(u,v),t(u,v)),Qt(r(u,v),t(u,v))> = 0
*              <S(u,v)-Q(r(u,v),t(u,v)),Qr(r(u,v),t(u,v))> = 0
*              This means that (derivation by u and v)
*              <Su-[Qt*tu+Qr*ru],Qt> + <S-Q,Qtt*tu + Qtr*ru> = 0
*              <Su-[Qt*tu+Qr*ru],Qr> + <S-Q,Qrt*tu + Qrr*ru> = 0
*              <Sv-[Qt*tv+Qr*rv],Qt> + <S-Q,Qtt*tv + Qtr*rv> = 0
*              <Sv-[Qt*tv+Qr*rv],Qr> + <S-Q,Qrt*tv + Qrr*rv> = 0
*              Solving these four equations gives us ru,rv,tu,tv. 
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, JANUARY 1991
*
*********************************************************************
*/                       
{                        


  int ki;                             /* Loop control.                               */
  int kdim = 3;                       /* Dim of object space.                        */
  double *sval;                       /* Pointer to first surface value              */
  double *s_u,*s_v,*s_uu,*s_uv,*s_vv; /* Pointers to first surface derivatives       */
  double *ns;                         /* Pointer to first surface normal             */
  double *qval;                       /* Pointer to second surface value             */
  double *q_t,*q_r,*q_tt,*q_tr,*q_rr; /* Pointer to second surface derivatives       */
  double *nq;                         /* Pointer to second surface normal            */
  double nq_u[3], nq_v[3];            /* Derivatives of second surface normal (with u and v !) */
  double help1[3], help2[3];          /* Help arrays                                 */
  double help3[3], help4[3];          /* Help arrays                                 */
  double matr[4];                     /* Matrix in linear equation to be solved      */
  int    piv[2];                      /* Pivotation array                            */
  double sq[3];                       /* The difference cevtor S-Q                   */
  double h_u[2];                      /* The partial derivative of h() by u          */
  double h_v[2];                      /* The partial derivative of h() by v          */
  int kstat;                          /* Local status                                */
  
  /* ------------------------------------------------------------------------------- */
  
  cdiff[0] = DZERO;
  cdiff[1] = DZERO;
  cdiff[2] = DZERO;
  cdiff[3] = DZERO;

  /* Init, Set pointers to input values */
  sval = evals;
  qval = evalq;
  
  s_u   = sval + kdim;
  s_v   = s_u   + kdim;
  s_uu  = s_v   + kdim;
  s_uv  = s_uu  + kdim;
  s_vv  = s_uv  + kdim;
  ns    = s_vv  + kdim;

  q_t   = qval + kdim;
  q_r   = q_t   + kdim;
  q_tt  = q_r   + kdim;
  q_tr  = q_tt  + kdim;
  q_rr  = q_tr  + kdim;
  nq    = q_rr  + kdim;

  /* Get the difference vector S-Q */
  s6diff(sval,qval,kdim,sq);
  
  /* Find the derivatives of the h() function by solving 2 2x2 systems (same matrix) */
  matr[0] = s6scpr(q_tt,sq,kdim) - s6scpr(q_t,q_t,kdim);
  matr[1] = s6scpr(q_tr,sq,kdim) - s6scpr(q_t,q_r,kdim);
  matr[2] = matr[1];
  matr[3] = s6scpr(q_rr,sq,kdim) - s6scpr(q_r,q_r,kdim);  

  h_u[0] = -s6scpr(s_u,q_t,kdim);
  h_u[1] = -s6scpr(s_u,q_r,kdim);

  h_v[0] = -s6scpr(s_v,q_t,kdim);
  h_v[1] = -s6scpr(s_v,q_r,kdim);
  

  /* Factorize matrix */
  s6lufacp(matr,piv,2,&kstat);
  if (kstat != 0) goto out;
  
  /* Solve */
  s6lusolp(matr,h_u,piv,2,&kstat);
  if (kstat != 0) goto out;

  /* Solve */
  s6lusolp(matr,h_v,piv,2,&kstat);
  if (kstat != 0) goto out;

  /* Construct matrix for finding du and dv */
  for (ki=0;ki<kdim;ki++) 
    {
      help1[ki] = q_tt[ki]*h_u[0] + q_tr[ki]*h_u[1];
      help2[ki] = q_tr[ki]*h_u[0] + q_rr[ki]*h_u[1];
    }
  s6crss(help1,q_r,help3);
  s6crss(q_t,help2,help4);
  
  for (ki=0;ki<3;ki++) nq_u[ki] = help3[ki] + help4[ki];

  for (ki=0;ki<kdim;ki++) 
    {
      help1[ki] = q_tt[ki]*h_v[0] + q_tr[ki]*h_v[1];
      help2[ki] = q_tr[ki]*h_v[0] + q_rr[ki]*h_v[1];
    }
  s6crss(help1,q_r,help3);
  s6crss(q_t,help2,help4);
  
  for (ki=0;ki<3;ki++) nq_v[ki] = help3[ki] + help4[ki];

  for (ki=0;ki<4;ki++) matr[ki] = DZERO;
  
  for (ki=0;ki<3;ki++) 
    {
      matr[0] += s_uu[ki]*nq[ki] + s_u[ki]*nq_u[ki];
      matr[1] += s_uv[ki]*nq[ki] + s_u[ki]*nq_v[ki];
      matr[2] += s_uv[ki]*nq[ki] + s_v[ki]*nq_u[ki];
      matr[3] += s_vv[ki]*nq[ki] + s_v[ki]*nq_v[ki];
    }
  
  /* solve the linear 2x2 system */

  s6lufacp(matr,piv,2,&kstat);
  if (kstat != 0) 
    {
      if( DNEQUAL(matr[0],DZERO)) cdiff[0] = - s6scpr(s_u,nq,kdim)/matr[0];
      else if( DNEQUAL(matr[1],DZERO)) cdiff[1] = - s6scpr(s_u,nq,kdim)/matr[1];
      else if( DNEQUAL(matr[2],DZERO)) cdiff[0] = - s6scpr(s_v,nq,kdim)/matr[2];
      else if( DNEQUAL(matr[3],DZERO)) cdiff[1] = - s6scpr(s_v,nq,kdim)/matr[3];

    }
  else 
    {
      cdiff[0] = - s6scpr(s_u,nq,kdim);
      cdiff[1] = - s6scpr(s_v,nq,kdim);
      s6lusolp(matr,cdiff,piv,2,&kstat);
    }
  out:;
  
}
