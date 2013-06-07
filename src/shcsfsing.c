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
 * $Id: shcsfsing.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */

#define SHCSFSING

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void shcsfsing_s9corr(double [], double [], double[]);
static void shcsfsing_s9dir(double [], double [], double[]);
#else
static void shcsfsing_s9corr();
static void shcsfsing_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  shcsfsing(SISLCurve *pcurve,SISLSurf *psurf,double limit[],
	    double enext[], double gpos[],int *jstat)
#else
void shcsfsing(pcurve,psurf,limit,enext,gpos,jstat)
     SISLCurve        *pcurve;
     SISLSurf         *psurf;
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
* PURPOSE    : To find one magical points.
*
* INPUT      : pcurve - Pointer to curve
*              psurf  - Pointer to surface
*              limit  - Parameter boarders of both geometries
*                       limit[0] - limit[1] Parameter interval curve.
*                       limit[2] - limit[3] Parameter interval surface 1.dir
*                       limit[4] - limit[5] Parameter interval surface 2.dir
*              enext    - Parameter start value for iteration(3 values).
*
*
* OUTPUT     : gpos    - Parameter values of the found singularity.(3 values)
*              jstat   - status messages  
*                                = 1   : Extremum found.
*                                = 0   : Extremum NOT found.
*                                < 0   : error.
*
*
* METHOD     :  - Start with a guess value (u,v) in domain of curve (S(w))
*           (a) - Find domain value (r,t) of closest point (to S(w) 
*                 in surface(Q(r,t))
*               - If vf(u,v) = <Sw,Normal(Q)> is small enough  stop 
*                      (<,> means scalar prod.)
*               - Find dwa by taylorizing vf.
*                 This include finding the derivatives of the closest 
*                 point function (r(w),t(w))
*                 with respect to w. (called h(w) in article, see coments 
*                 in shcsfsing_s9dir)
*               - w:= w+dw, goto (a)
*
*
* REFERENCES : Solutions of tangential surface and curve intersections.
*              R P Markot and R L Magedson
*              Computer-Aided Design; vol. 21, no 7 sept. 1989, page 421-429
*
*
* WRITTEN BY : Ulf J. Krystad, SI, SEPTEMBER 1992
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                       */
  int kpos = 0;             /* Position of error.                           */
  int ki;                   /* Loop control                                 */
  int kleft=0;             /* Variables used in the evaluator.              */
  int kleftu=0;             /* Variables used in the evaluator.             */
  int kleftv=0;             /* Variables used in the evaluator.             */
  int kder=2;               /* Order of derivatives to be calulated         */
  int kdim=3;               /* Dimension of space the surface lies in       */
  int knbit;                /* Number of iterations                         */
  double tdelta[3];         /* Length of parameter intervals.               */
  double tdist;             /* The current scalar product.                  */
  double tprev;             /* The previous scalar product.                 */
  double td[3],t1[3],tdn[3];/* Distances between old and new parameter      */
			    /* value in the three parameter directions.     */
  double crv_val[21];         /* Value ,first and second derivatiev of surf.*/ 
  double *crv_tang=crv_val+3;  /* Curve tangent                             */
  double surf_val[21];         /* Value ,first and second deriv. of surf.   */ 
  double *surf_norm=surf_val+18;  /* Normal vector of the surface           */
  double snext[3];          /* Parameter values                             */
  double start[2];          /* Parameters limit of surface, used in         */
                            /* call to closest point                        */
  double end[2];            /* Parameters limit of surface, used in         */
                            /* call to closest point                        */
  double guess[2];          /* Start point for closest point iteration      */
  double tol = (double)10000.0*REL_COMP_RES; /* equality tol. in par.space  */
  SISLPoint *ppoint=SISL_NULL;   /* Contains the current position of the curve   */ 
                            /* used in closest point iteration              */
  int max_iter=20;          /* Maximal number of iteration allowed          */

  /* ---------------------------------------------------------------------- */
  
  /* Test input.  */
  if (pcurve->idim != kdim) goto err106;
  if (psurf->idim != kdim) goto err106;
  
  /* Fetch referance numbers from the serach intervals for the surfaces.  */
  tdelta[0] = limit[1] - limit[0];
  tdelta[1] = limit[3] - limit[2];
  tdelta[2] = limit[5] - limit[4];

  /* Set limit values, used in closest point iteration */
  start[0] = limit[2];
  start[1] = limit[4];
  end[0]   = limit[3];
  end[1]   = limit[5];
  
  /* Create point, used in closest point iteration */
  ppoint = newPoint(crv_val,3,0);
  
  /* Collapsed ? */
  for (ki=0;ki<3;ki++) if (tdelta[ki] < tol) goto errsmall;  
  
  /* Initiate output variables.  */
  for (ki=0;ki<3;ki++)     gpos[ki] = enext[ki];

  /* Evaluate 0.-2. derivatives ofthe curve */
  s1221(pcurve,kder,gpos[0],&kleft,crv_val,&kstat);
  if (kstat < 0) goto error;

  /* Get closest point in second surface. */
  guess[0] = gpos[1];
  guess[1] = gpos[2];
  s1773(ppoint,psurf,REL_COMP_RES,start,end,guess,gpos+1,&kstat);
  if (kstat < 0) goto error;
  
  /* Evaluate 0.-2. derivatives of surface */
  s1421(psurf,kder,gpos+1,&kleftu,&kleftv,surf_val,surf_norm,&kstat);
  if (kstat < 0) goto error;

  /* Get length of dot product */
  tprev = fabs(s6scpr(crv_tang,surf_norm,kdim));
  
  /* Compute the Newton stepdistance vector in first surface. */
  shcsfsing_s9dir(td,crv_val,surf_val);
  
  /* Adjust if we are not inside the parameter intervall. */
  for (ki=0;ki<3;ki++)    t1[ki] = td[ki];

  shcsfsing_s9corr(t1,gpos,limit);
  
  /* Iteratation loop.  */
  
  for (knbit = 0; knbit < max_iter; knbit++)
    {
      
      snext[0] = gpos[0] + t1[0];
   
      /* Evaluate 0.-2. derivatives of curve */
      s1221(pcurve,kder,snext[0],&kleft,crv_val,&kstat);
      if (kstat < 0) goto error;
      
      /* Get closest point in surface. */
      guess[0] = gpos[1];
      guess[1] = gpos[2];
      s1773(ppoint,psurf,REL_COMP_RES,start,end,guess,snext+1,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate 0.-2. derivatives of surface */
      s1421(psurf,kder,snext+1,&kleftu,&kleftv,surf_val,surf_norm,&kstat);
      if (kstat < 0) goto error;

      /* Get length of dot product */
      tdist = fabs(s6scpr(crv_tang,surf_norm,kdim));
 
      /* Compute the Newton stepdistance vector. */
      shcsfsing_s9dir(tdn,crv_val,surf_val);
      
      if (tdist <= tprev)
	{
	  /* Ordinary converging. */
	  
	  for (ki=0;ki<3;ki++) gpos[ki] = snext[ki];
	  td[0] = t1[0] = tdn[0];

	  /* Adjust if we are not inside the parameter intervall. */
	  shcsfsing_s9corr(t1,gpos,limit);
		  
          if ((fabs(t1[0]/tdelta[0]) <= REL_COMP_RES))
	      {
		gpos[0] += t1[0];
		/* Evaluate 0.-2. derivatives of curve */
		s1221(pcurve,kder,gpos[0],&kleft,crv_val,&kstat);
		if (kstat < 0) goto error;
		
		/* Get closest point in surface. */
		guess[0] = gpos[1];
		guess[1] = gpos[2];
		s1773(ppoint,psurf,REL_COMP_RES,start,end,guess,gpos+1,&kstat);
		if (kstat < 0) goto error;
		break;
	      }
	  tprev = tdist;
	}
      
      else
	{
	  /* Not converging, half step length try again. */
	   
	   t1[0] /= (double)2;
	}
    }
  
  /* Iteration stopped, test if point is extremum */
  /* Unsure about what is right here */
  if (tdist <= tol)
    *jstat = 1;
  else
    *jstat = 0;

 
  /* Test if the iteration is close to a knot */
  if (fabs(gpos[0] - pcurve->et[kleft])/tdelta[0] < tol)
    gpos[0] = pcurve->et[kleft];
  else if (fabs(gpos[0] - pcurve->et[kleft+1])/tdelta[0] < tol)
    gpos[0] = pcurve->et[kleft+1];

  if (fabs(gpos[1] - psurf->et1[kleftu])/tdelta[1] < tol)
    gpos[1] = psurf->et1[kleftu];
  else if (fabs(gpos[1] - psurf->et1[kleftu+1])/tdelta[1] < tol)
    gpos[1] = psurf->et1[kleftu+1];
  
  if (fabs(gpos[2] - psurf->et2[kleftv])/tdelta[2] < tol)
    gpos[2] = psurf->et2[kleftv];
  else if (fabs(gpos[3] - psurf->et2[kleftv+1])/tdelta[2] < tol)
    gpos[2] = psurf->et2[kleftv+1];
  
  /* Iteration completed.  */
  goto out;
  
  /* --------------------------------------------------------------------- */ 
  /* Error in input. Dimension not equal to 3 */
  err106: *jstat = -106;
  s6err("shcsfsing",*jstat,kpos);
  goto out;                  

  /* Error in input. One parameter interval colapsed. */
  errsmall: *jstat = -200;
  s6err("shcsfsing",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("shcsfsing",*jstat,kpos);
  goto out;                  
  
 out:if(ppoint) freePoint(ppoint);
}

#if defined(SISLNEEDPROTOTYPES)
static void
  shcsfsing_s9corr(double gd[], double coef[],double limit[])
#else
static void shcsfsing_s9corr(gd,coef,limit)
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
* WRITTEN BY : Ulf J. Krystad, SI, SEPTEMBER 1992
*
*********************************************************************
*/                       
{
  int ki;

  for (ki=0;ki<3;ki++)
    if (coef[ki] + gd[ki] < limit[2*ki])        gd[ki] = limit[2*ki]    - coef[ki];
    else if (coef[ki] + gd[ki] > limit[2*ki+1]) gd[ki] = limit[2*ki +1] - coef[ki];
  

}

#if defined(SISLNEEDPROTOTYPES)
static void
  shcsfsing_s9dir(double cdiff[],double evals[],double evalq[])
#else
static void shcsfsing_s9dir(cdiff,evals,evalq)
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
* INPUT      : evals - Value and derivatives  on curve
*              evalq - Value and derivatives  on surface.
*              
*
* OUTPUT     : cdiff1  - Parameter increments in one direction.
*            
*
*
* METHOD     : See comments in main header.
*              The only thing missing in the article is the derivation of the
*              function h(w). Calling this function (r(w), t(w))
*              we know that
*              For all w (<,> meaning scalar product)
*              <S(w)-Q(r(w),t(w)),Qt(r(w),t(w))> = 0
*              <S(w)-Q(r(w),t(w)),Qr(r(w),t(w))> = 0
*              This means that (derivation by w)
*              <Sw-[Qt*tw+Qr*rw],Qt> + <S-Q,Qtt*tw + Qtr*rw> = 0
*              <Sw-[Qt*tw+Qr*rw],Qr> + <S-Q,Qrt*tw + Qrr*rw> = 0
*              Solving these two equations gives us rw,tw. 
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, SEPTEMBER 1992
*
*********************************************************************
*/                       
{                        


  int ki;                          /* Loop control.                               */
  int kdim = 3;                       /* Dim of object space.                        */
  double *sval;                       /* Pointer to first surface value              */
  double *s_w,*s_ww;                  /* Pointers to curve derivatives       */
  double *qval;                       /* Pointer to surface value             */
  double *q_t,*q_r,*q_tt,*q_tr,*q_rr; /* Pointer to surface derivatives       */
  double *nq;                         /* Pointer to surface normal            */
  double nq_w[3];            	      /* Derivatives of surface normal (with w!) */
  double help1[3], help2[3];          /* Help arrays                                 */
  double help3[3], help4[3];          /* Help arrays                                 */
  double matr[4];                     /* Matrix in linear equation to be solved      */
  int    piv[2];                      /* Pivotation array                            */
  double sq[3];                       /* The difference cevtor S-Q                   */
  double h_w[2];                      /* The derivative of h() by w          */
  double f_val, f_deriv;              /* Value and derivative of the 1 D function. */
   int kstat;                          /* Local status                                */
  
  /* ------------------------------------------------------------------------------- */
  
  cdiff[0] = DZERO;
  cdiff[1] = DZERO;
  cdiff[2] = DZERO;

  /* Init, Set pointers to input values */
  sval = evals;
  qval = evalq;
  
  s_w   = sval + kdim;
  s_ww  = s_w   + kdim;

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

  h_w[0] = -s6scpr(s_w,q_t,kdim);
  h_w[1] = -s6scpr(s_w,q_r,kdim);
  

  /* Factorize matrix */
  s6lufacp(matr,piv,2,&kstat);
  if (kstat != 0) goto out;
  
  /* Solve */
  s6lusolp(matr,h_w,piv,2,&kstat);
  if (kstat != 0) goto out;


  /* Construct matrix for finding dw */
  s6crss(q_tt,q_r,help3);
  s6crss(q_t,q_tr,help4);
  
  for (ki=0;ki<3;ki++) help1[ki] = (help3[ki] + help4[ki])*h_w[0];

  s6crss(q_tr,q_r,help3);
  s6crss(q_t,q_rr,help4);
  
  for (ki=0;ki<3;ki++) help2[ki] = (help3[ki] + help4[ki])*h_w[1];

  for (ki=0;ki<3;ki++) nq_w[ki] = (help1[ki] + help2[ki]);
  
  f_val = s6scpr(s_w, nq, kdim);
  f_deriv = s6scpr(s_ww, nq, kdim) + s6scpr(s_w, nq_w, kdim);
  
  if (DNEQUAL(f_val + fabs(f_deriv), f_val))
      cdiff[0] = - f_val/f_deriv;
     
  out:;
  
}
