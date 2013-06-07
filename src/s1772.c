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
 * $Id: s1772.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */
#define S1772

#include "sislP.h"

#define SINGULAR 1.0e-16
#define copy2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]=(b)[ki]
#define copy3(a,b,c,d) for (ki=0;ki<(d);ki++) (a)[ki]=(b)[ki]=(c)[ki]
#define incr2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]+=(b)[ki]
#define decr2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]-=(b)[ki]
#define set_order(a) if((a)==1) {s_v=s_uu;order=0;} else {s_v=s_v1;order=1;}
/* extern int debug_flag; */
/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static
void
s1772_s9corr(double[],double[],double,double,double[],double[],int*);
static
void
s1772_s9dir(double*,double[],double[],double[],double[],double[],double[],
	    double[],double[],double[],double[],double[],int,int,int*);
static
void
s1772_s6sekant1(SISLCurve *,SISLSurf *,double[],double,double *,double,double,
		double[],double,double[],double[],double[],double[],int *);
static
int
s1772_s6local_pretop(double,double[],double[],double[],double[],double[],
		     double[],double[],double[],double[],double[],double[],
		     int,int*);
#else
static void s1772_s9corr();
static void s1772_s9dir();
static void s1772_s6sekant1();
static int s1772_s6local_pretop();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1772(SISLCurve *pcurve,SISLSurf *psurf,double aepsge,
	   double astart1,double estart2[],double aend1,double eend2[],
	   double anext1,double enext2[],double *cpos1,double gpos2[],
	   int *jstat)
#else
void s1772(pcurve,psurf,aepsge,astart1,estart2,aend1,eend2,
	   anext1,enext2,cpos1,gpos2,jstat)
     SISLCurve  *pcurve;
     SISLSurf   *psurf;
     double aepsge;
     double astart1;
     double estart2[];
     double aend1;
     double eend2[];
     double anext1;
     double enext2[];
     double *cpos1;
     double gpos2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a surface to find a closest point
*              or an intersection point.
*
*
* INPUT      : pcurve    - Pointer to the curve in the intersection.
*              psurf     - Pointer to the surface in the intersection.
*              aepsge    - Geometry resolution.
*              astart1   - Start parameter value of the curve.
*              estart2[] - Start parameter value of surface.
*              aend1     - End parameter value of the curve.
*              eend2[]   - End 1.st parameter value of the surface.
*              anext1    - Start parameter value of the iteration on
*                          the curve.
*              enext2[]  - Start parameter value of the iteration on
*                          the surface.
*              jstat     -	=1   : A quick version is to be used.  
*				else : Importent to find a solution.
*
*
* OUTPUT     : cpos1   - Parameter value of the curve in intersection
*                        point.
*              gpos2[] - Parameter value of the surface in
*                        intersection point.
*              jstat   - status messages  
*                                = 3   : A minimum distanse found.
*                                = 2   : Nothing found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in tree parameter directions.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Nov 1992
*
*********************************************************************
*/                       
{
  int ki;		    /* Counter.					   */
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int left[3];              /* Variables used in the evaluator.            */
  int dim;                  /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int p_dir;                /* Changing direction in par-space.            */
  int g_up,ng_up,g_dir;     /* Changing direction in geometric space.      */
  int order;		    /* Order of methode.			   */
  int sing = 0;		    /* Mark that singularity has ocured.	   */	
  double *c0=SISL_NULL;          /* Value  of curve.				   */ 
  double *c_t;		    /* First derivatiev of curve.		   */ 
  double *c_tt;		    /* Second derivatiev of curve.		   */ 
  double *s0;               /* Value of surf. 				   */
  double *s_u;		    /* First derivatiev in first dir of surf. 	   */
  double *s_v;		    /* First derivatiev in second dir of surf.	   */
  double *s_uu;		    /* Second derivatiev in first dir of surf. 	   */
  double *s_vv;		    /* Second derivatiev in second dir of surf.	   */
  double *s_uv;		    /* Cross derivatiev of surf.	   	   */
  double *s_v1;		    /* First derivatiev in second dir of surf.	   */
  double *norm;		    /* Normal to the surface.			   */
  double *diff;             /* Difference between the curve and the surf.  */
  double *prev_diff;        /* Previous difference.			   */
  double delta[3];          /* Parameter interval of the curve and surface.*/
  double d[3];		    /* Clipped distances between old and new par.
			       value in the tree parameter directions.     */
  double c_d[3];	    /* Computed distances ....			   */
  double nc_d[3];	    /* New computed distances ....		   */
  double dist;              /* Distance between position and origo.        */
  double prev_dist;         /* Previous difference between the curves.     */
  double par_val[3];        /* Parameter values                            */
  double local[45];
  int corr = 0, div2 = 0;


  
  
  /* Test input.  */
  
  if (pcurve->idim != psurf->idim) goto err106;  
  dim = pcurve->idim;
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  delta[0] = psurf->et1[psurf->in1] - psurf->et1[psurf->ik1 - 1];
  delta[1] = psurf->et2[psurf->in2] - psurf->et2[psurf->ik2 - 1];
  delta[2] = pcurve->et[pcurve->in] - pcurve->et[pcurve->ik - 1];
  
  /* Allocate local used memory and set value pointers.*/

  if (dim > 3)
  {
     c0 = newarray((15)*dim,double);
     if (c0 == SISL_NULL) goto err101;
  }
  else
     c0 = local;
  
  s0 = c0 + 3*dim;
  diff = s0 + 10*dim;
  prev_diff = diff+dim;
  c_t = c0+dim;
  c_tt = c_t+dim;
  s_u = s0+dim;
  s_uu = s_u+dim;
  s_v1 = s_uu+dim;
  s_uv  = s_v1+dim;
  s_vv = s_uv+dim+dim;
  norm = s_vv+dim;
    
  /* Initiate variables.  */

  copy2(par_val,enext2,2);
  par_val[2] = anext1;
  left[0]=left[1]=left[2]=0;  
  
  for (ki=1; ki<3; ki++)
  {
     set_order(ki);
     
     /* Evaluate 0-2.st derivatives of curve */
     
     if (par_val[2] == aend1)
	s1227(pcurve,1+order,par_val[2],left+2,c0,&kstat);
     else
	s1221(pcurve,1+order,par_val[2],left+2,c0,&kstat);
     if (kstat < 0) goto error;
     
     /* Evaluate 0-2.st derivatives of surface */
     
     s1424(psurf,1+order,1+order,par_val,left,left+1,s0,&kstat);
     if (kstat < 0) goto error;
     
     /* Compute the distanse vector and value and the new step. */
     
     s1772_s9dir(&dist,diff,c_d, c0,c_t,c_tt,
		 s0,s_u,s_v,s_uu,s_uv,s_vv, dim,order,&kstat);
     if (kstat < 0) goto error;
     if (kstat == 1) 		/* Singular matrix. */
     {
	if (order == 1) goto singular;
     }
     else break;
  }
  
  /* Correct if we are not inside the parameter intervall. */
  
  s6crss(s_u,s_v,norm);
  g_up = ((s6scpr(diff,norm,dim) >= DZERO) ? 1 : -1);
  copy2(d,c_d,3);
  s1772_s9corr(d,par_val, astart1,aend1,estart2,eend2,&corr);
  prev_dist = dist;
  copy2(prev_diff,diff,dim);      
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 30; knbit++)
  {
     incr2(par_val,d,3);
     
     while (1)
     {
	/* Evaluate 0-2.st derivatives of curve */
	
	if (par_val[2] == aend1)
	   s1227(pcurve,1+order,par_val[2],left+2,c0,&kstat);
	else
	   s1221(pcurve,1+order,par_val[2],left+2,c0,&kstat);
	if (kstat < 0) goto error;
	
	/* Evaluate 0-2.st derivatives of surface */
	
	s1424(psurf,1+order,1+order,par_val,left,left+1,s0,&kstat);
	if (kstat < 0) goto error;
	
	/* Compute the distanse vector and value and the new step. */
	
	
	s1772_s9dir(&dist,diff,nc_d, c0,c_t,c_tt,
		    s0,s_u,s_v,s_uu,s_uv,s_vv, dim,order,&kstat);
	if (kstat < 0) goto error;      
	if (kstat == 1) 		/* Singular matrix. */
	{
	   sing++;
	   if (order == 1) goto singular;
	   else	 set_order(2);		/* Change to order 2. */
	}
	else
	{
	   s6crss(s_u,s_v,norm);
	   ng_up = ((s6scpr(diff,norm,dim) >= DZERO) ? 1 : -1);
	   
	   g_dir = (ng_up+g_up != 0);			/* 0 if changed. */
	   p_dir = (s6scpr(c_d,nc_d,3) >= DZERO);	/* 0 if changed. */
	   
	   if (!order && g_dir && (!p_dir || dist > 0.3*prev_dist))
	   {
	      if (div2) div2 = 0;
	      set_order(2);
	      /*  if (debug_flag) printf("\n order-2 ");*/
	   }
	   else if (order && !g_dir)
	   {
	      if (sing) goto singular;
	      if (div2) div2 = 0;
	      set_order(1);
	      /*  if (debug_flag) printf("\n  order-1 "); */
	   }
	   else
	   {
	      if (sing) sing = 0;
	      break;
	   }
	}
     }
     
     if (corr)
	if (!(p_dir && g_dir)) corr = 0;

     if (dist < prev_dist)
     {
	if (div2) div2 = 0;
	
	/* Corrigate if we are not inside the parameter intervall. */
	
	g_up = ng_up;
	copy3(d,c_d,nc_d,3);
	s1772_s9corr(d,par_val, astart1,aend1,estart2,eend2,&corr);
	prev_dist = dist;
	copy2(prev_diff,diff,dim);
	
	/* Testing */
	/*	if (quick && corr > 3) break; */
	if (corr > 3) break;
     }    
     else if ( corr > 3 ||
	     ((fabs(d[0]/delta[0]) <= REL_COMP_RES) &&
	      (fabs(d[1]/delta[1]) <= REL_COMP_RES) &&
	      (fabs(d[2]/delta[2]) <= REL_COMP_RES)))     break;
     else
     {
	/* Not converging, corrigate and try again. */
	/*  if (debug_flag) printf(" *h*:%d ",knbit);*/

	if (corr) corr++;
	if (dist > prev_dist && div2) break;
	div2++;
        decr2(par_val,d,3);
	d[0] /= 2; d[1] /= 2; d[2] /= 2;
     }
  }
	
	/* Iteration stopped, test if point found is within resolution */
  
  goto not_singular;
  
singular:
   
   /*  if (!quick && dist > aepsge) */
     if (dist > aepsge)
     {
	ki = s1772_s6local_pretop(dist,diff,norm,c0,c_t,c_tt,
			    s0,s_u,s_v,s_uu,s_uv,s_vv,dim,&kstat);
	if (kstat < 0) goto error;
	if (ki == 0)
	{
	   s1772_s6sekant1(pcurve,psurf,par_val,c_d[2],&dist,aepsge,
			   astart1,estart2,aend1,eend2,c0,s0,norm,&kstat);  
	   if (kstat < 0) goto error;
	}
     }
     
not_singular:	
  if (dist <= aepsge)
  {
     /* if (debug_flag) printf("\n FOUND: %d dist = %g",knbit,dist); */

    *jstat = 1;
  }
  else
  {
     /*if (debug_flag) printf("\n no: %d dist = %g",knbit,dist);*/

     s6crss(s_u,s_v,norm);
     if ((PIHALF-s6ang(c_t,norm,dim)) < ANGULAR_TOLERANCE)
	*jstat = 3;
     else
	*jstat = 2;
  }
  
  /* if (knbit > 25)
     if (debug_flag) printf("\n *****status: %d dist: %f \tknbit: %d",
			    *jstat,dist,knbit); */

  *cpos1 = par_val[2];
  gpos2[0] = par_val[0];
  gpos2[1] = par_val[1];
  
  /* Iteration completed.  */
  
  goto out;
  
  /* Error in allocation */
  
  err101: 
    *jstat = -101;
    s6err("s1772",*jstat,kpos);
    goto out;                  
  
  /* Error in input. Conflicting dimensions.  */
  
  err106: 
    *jstat = -106;
    s6err("s1772",*jstat,kpos);
    goto out;                  
  
  /* Error in lower level routine.  */
  
  error : 
    *jstat = kstat;
    s6err("s1772",*jstat,kpos);
    goto out;                  
  
  out:
    if (c0 != local && c0 != SISL_NULL) freearray(c0);
}

#if defined(SISLNEEDPROTOTYPES)
static
void
s1772_s9corr(double gd[],double acoef[],double astart1,double aend1,
	     double astart2[],double aend2[],int *corr)
#else
static void s1772_s9corr(gd,acoef,astart1,aend1,astart2,aend2,corr)
     double gd[];
     double acoef[];
     double astart1;
     double aend1;
     double astart2[];
     double aend2[];
     int *corr;
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
* INPUT      : acoef[]   - Parameter values.
*              astart1   - The lower boorder in curve.
*              aend1     - The higher boorder in curve.
*              estart2[] - The lower boorder in surface.
*              eend2[]   - The higher boorder in surface.
*
*
*
* INPUT/OUTPUT : gdn   - Old and new step value.
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
* WRITTEN BY : Arne Laksaa, SI, Nov 1992.
*
*********************************************************************
*/                       
{
  int lcorr = 0;
  if (acoef[0] + gd[0] < astart2[0])  
    {
       gd[0] = astart2[0] - acoef[0]; 
       lcorr=1;
    }
  else if (acoef[0] + gd[0] > aend2[0]) 
    {
       gd[0] = aend2[0] - acoef[0]; 
       lcorr=1;
    }
  
  if (acoef[1] + gd[1] < astart2[1])  
    {
       gd[1] = astart2[1] - acoef[1]; 
       lcorr=1;
    }
  else if (acoef[1] + gd[1] > aend2[1]) 
    {
       gd[1] = aend2[1] - acoef[1]; 
       lcorr=1;
    }
  
  if (acoef[2] + gd[2] < astart1)  
    {
       gd[2] = astart1 - acoef[2]; 
       lcorr=1;
    }
  else if (acoef[2] + gd[2] > aend1) 
    {
       gd[2] = aend1 - acoef[2]; 
       lcorr=1;
    }
  
  if (lcorr) 
    (*corr)++;
  else 
    (*corr) = 0;
}
   
#if defined(SISLNEEDPROTOTYPES)
   static
      void
s1772_s9dir(double *dist,double diff[],double delta[],
		  double f[],double f_t[],double f_tt[],
		  double g[],double g_u[],double g_v[],
		  double g_uu[],double g_uv[],double g_vv[],
		  int dim,int second,int* jstat)
#else
static void s1772_s9dir(dist,diff,delta,f,f_t,f_tt,g,g_u,g_v,g_uu,
			g_uv,g_vv,dim,second,jstat)
     double *dist;
     double diff[];
     double delta[];
     double f[];
     double f_t[];
     double f_tt[];
     double g[];
     double g_u[];
     double g_v[];
     double g_uu[];
     double g_uv[];
     double g_vv[];
     int    dim;
     int    second;
     int*   jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance vector and the length of the
*              distance vector.
*	       And to compute a next step on all tree parameter direction
*	       towards an intersection or a closest point.
*
*
* INPUT      : f    - Value in point on curve.
*              f_t  - First derevative of the curve.
*              f_tt - Second derevative of the curve.
*              g    - Value in point on surface.
*	       g_u  - Derevative of the surface in first par-dir.
*	       g_u  - Derevative of the surface in second par-dir.
*	       g_uu - Second derevative of the surface in first par-dir.
*	       g_uv - Cross derevative of the surface.
*	       g_vv - Second derevative of the surface in second par-dir.
*	       dim    - Dimension of space the curve/surface lie in.
*              second - 1 if we have to use second order method 
*			0 if only first order method.
*
* OUTPUT     : dist    - The lengt of the different vector.
*	       diff[]  - The differens vector.
*	       delta[] - Relative step parameter values towards intersection on
*                        the curve delta[2] and the surface (delta[0],delta[1]).
*
*
* METHOD     : We have a point on the curve and a point on the surface.
*	       The distance vector between this two point are going towards
*	       eider zero length or to be orthogonal to the derivative on
*	       the curve and the derivatives on the surface.
*	       The method is to use Newton itarations.
*	       The dot product between the distance vecore and the derivatives
*	       are going to be zero.
*	       We then use Taylor expasion bouth on the position and
*	       the derivatives. Then we just remove the second degree parts
*	       of the eqations. The matrix is then splited into a first order
*	       part A:
*		     | -g_u |	
*		 K = | -g_v |    ,and       A = K*K(transpost)
*		     |  f_t | 	
*
*	         where f_t is the first derivative of the curve,
*	         and  g_u is the first derivative of the surface in first dir,
*	         and  g_v is the first derivative of the surface in second dir,
*
*	       and a second order part B:
*		     | -<d,g_uu>   -<d,g_uv> 		 |
*		 B = | -<d,g_uv>   -<d,g_vv>		 |
*		     |  			<d,f_tt> |
*
*	         where d is the distance vector,
*	         and f_tt is second derivative of the curve,
*	         and  g_uu is second derivative of the surface in first dir,
*	         and  g_vv is second derivative of the surface in second dir,
*	         and  g_uv is cross derivative of the surface.
*
*	       We then has the following possible eqations:
*
*	                 (A+B)*delta = -K*d, or 
*		             A*delta = -K*d, or
*	          K(transpost)*delta = -d.
*
*
*	       The solutions of these matrix equations are the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Nov 1992.
*
*********************************************************************
*/                       
{                        
  int kstat;			/* Local status variable. 		  */
  double a1,a2,a3,a4,a5,a6;	/* The A matrix, diagonal and A12 A13 A23.*/
  double b1,b2,b3,b4;		/* The B matrix, diagonal and B23.	  */
  double A[9],mat[9];		/* Matrix in linear equation to be solved */
  double h[3];			/* Left side in the equation.		  */
  double x[3];			/* Left side in the equation.		  */
  double r[3];			/* Left side in the equation.		  */
  double det;			/* Determinant for matrix.		  */
  long double ss,aa,xx,bb;	/* For use in iterative improvement.      */
  int    piv[3];		/* Pivotation array                       */
  int k,k3,j;			/* Counters.				  */
  
  
  /* Computing the different vector */
  
  s6diff(f,g,dim,diff);
  
  /* Computing the length of the different vector. */
  
  *dist = s6length(diff,dim,&kstat);
  if (kstat<0) goto error;
  
  if (second || dim != 3)
  {
     a1 = s6scpr(f_t,f_t,dim);
     a2 = s6scpr(g_u,g_u,dim);
     a3 = s6scpr(g_v,g_v,dim);
     a4 = s6scpr(f_t,g_u,dim);
     a5 = s6scpr(f_t,g_v,dim);
     a6 = s6scpr(g_u,g_v,dim);
  }
  
  if (second)
  {
     b1 = s6scpr(diff,f_tt,dim);
     b2 = s6scpr(diff,g_uu,dim);
     b3 = s6scpr(diff,g_vv,dim);
     b4 = s6scpr(diff,g_uv,dim);
  }
  else b1=b2=b3=b4=0;

  if (second || dim != 3)
  {  
     mat[0] = a2-b2;	mat[1] = a6-b4;		mat[2] = -a4;
     mat[3] = a6-b4;	mat[4] = a3-b3;		mat[5] = -a5;
     mat[6] = -a4;	mat[7] = -a5;		mat[8] = a1+b1;
     
     h[0] =  s6scpr(diff,g_u,dim);
     h[1] =  s6scpr(diff,g_v,dim);
     h[2] = -s6scpr(diff,f_t,dim);
  }
  else
  {
     mat[0] = g_u[0];	mat[1] = g_v[0];	mat[2] = -f_t[0];
     mat[3] = g_u[1];	mat[4] = g_v[1];	mat[5] = -f_t[1];
     mat[6] = g_u[2]; 	mat[7] = g_v[2];	mat[8] = -f_t[2];
     
     h[0] =  diff[0];
     h[1] =  diff[1];
     h[2] =  diff[2];
  }
  
  for (k=0;k<9;k++) A[k]=mat[k];
  for (k=0;k<3;k++) x[k]=h[k];
  
  det = A[0]*(A[4]*A[8]-A[5]*A[7])
      - A[1]*(A[3]*A[8]-A[5]*A[6])
      + A[2]*(A[3]*A[7]-A[4]*A[6]);
  if (fabs(det) < 1.0e-16)
  {
     *jstat = 1;
     goto out;
  }  
     
  /* solve the linear 3x3 system */

  /*  s1772_s6lufacp(mat,piv,&kstat); */
  s6lufacp(mat,piv,3,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }  
  
  s6lusolp(mat,x,piv,3,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }
  
  for (k=0;k<3;k++) delta[k] = x[k];

  for (k=k3=0; k<3; k++,k3+=3)
  {
     for (ss=0.0,j=0; j<3; j++)
     {
	aa = A[j+k3];
	xx = x[j];
	ss += aa*xx;
     }
     bb = h[k];
     ss = bb-ss;
     r[k] = ss;
  }
  s6lusolp(mat,r,piv,3,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }

  for (k=0;k<3;k++) delta[k] = x[k] + r[k];
  
  /* if (debug_flag) printf("\nITERATIV IMPROVES: r = (%g %g %g) ",
			 delta[0]-x[0],delta[1]-x[1],delta[2]-x[2]); */

  *jstat = 0;
  goto out;

  error : 
    *jstat = kstat;
    s6err("s1772_s9dir",*jstat,0);
    goto out;                  
	       
  out: 
    return;
}


#if defined(SISLNEEDPROTOTYPES)
   static
      void
s1772_s6sekant1(SISLCurve *pcurve,SISLSurf *psurf,
		 double  par_val[], double delta, double *dist, double aepsge,
		 double astart1,double estart2[],double aend1,double eend2[],
		 double c0[], double s0[], double norm[],
		 int *jstat)
#else
static void s1772_s6sekant1(pcurve,psurf,par_val,delta,dist,aepsge,
			   astart1,estart2,aend1,eend2,c0,s0,norm,jstat)
     SISLCurve  *pcurve;
     SISLSurf   *psurf;
     double  par_val[];
     double  delta;
     double  *dist;
     double aepsge;
     double astart1;
     double estart2[];
     double aend1;
     double eend2[];
     double c0[];
     double s0[];
     double norm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Sekant methode iteration on the distance function between
*              a curve and a surface to find a closest point
*              or an intersection point.
*
*
* INPUT      : pcurve    - Pointer to the curve in the intersection.
*              psurf     - Pointer to the surface in the intersection.
*              delta     - Parameter distance on curve beetveen start values.
*              aepsge    - Geometry resolution.
*              c0        - Array for use in evaluation.
*              s0        - Array for use in evaluation.
*              norm      - Array for use in evaluation.
* INPUT/
* OUTPUT     : par_val[] - Parameter value of the surface in
*                          intersection point.
*              dist      - Distance in space.
* OUTPUT     : jstat     - status messages  
*                                = 3   : A minimum distanse found.
*                                = 2   : Nothing found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Sekant mothode in tree parameter directions.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Nov 1992
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*
*********************************************************************
*/                       
{
  int ki,kj;		    /* Counter.					   */
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int dim;                  /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  double cu_val[2];	    /* Parameter values on curve.		   */
  double new_cu_val;	    /* New parameter value on curve.		   */
  double *diff;		    /* Difference vector between curve surface.    */
  double y[2],new_y,delta_y;/* Signed distance.				   */
  SISLPoint *pt=SISL_NULL;	    /* Point for use in closest point point/surface*/
  int cu_left = 0;	    /* Keep left knot information for evaluator.   */
  int s_left1 = 0;	    /* Keep left knot information for evaluator.   */
  int s_left2 = 0;	    /* Keep left knot information for evaluator.   */
  int shift = 0;	    /* Mark that the diriction have been changed.  */

  *jstat = 0;
  
  /* Test input.  */
  
  if (pcurve->idim != psurf->idim) goto err106;  
  dim = pcurve->idim;
  diff = c0 + dim;
   
  if ((pt = newPoint(c0,dim,0)) == SISL_NULL) goto err101;

  if (delta == 0.0) delta =1e-15;
  
  if ((par_val[2] == astart1 && delta < 0.0) ||
      (par_val[2] == aend1   && delta > 0.0))
  {
     delta = -delta;
     shift++;
  }
  
  if (fabs(delta) < (aend1 -astart1)/100.0)
  {
     if (delta < 0.0)
	delta = (astart1 - aend1)/100.0;
     else
	delta = (aend1 - astart1)/100.0;
  }
  else if (fabs(delta) > (aend1 -astart1)/10.0)
  {
     if (delta < 0.0)
	delta = (astart1 - aend1)/10.0;
     else
	delta = (aend1 - astart1)/10.0;
  }


  cu_val[0] = par_val[2];
  s1221(pcurve,0,cu_val[0],&cu_left,pt->ecoef,&kstat);
  if (kstat < 0) goto error;      
  s1773(pt,psurf,aepsge,estart2,eend2,par_val,par_val,&kstat);
  if (kstat < 0) goto error;      
  s1421(psurf,1,par_val,&s_left1,&s_left2,s0,norm,&kstat);
  if (kstat < 0) goto error;
  for(kj=0; kj<dim; kj++) diff[kj] = s0[kj] - pt->ecoef[kj];
  new_y = s6norm(norm,dim,norm,&kstat);
  if (kstat == 0)
  {
     (*dist)=s6length(diff,dim,&kstat);
     new_cu_val = cu_val[0];
     goto out;
  }
  if (((*dist)=s6length(diff,dim,&kstat)) < aepsge)
  {
     new_cu_val = cu_val[0];
     goto out;
  }
  y[0] = s6scpr(norm,diff,dim);
  cu_val[1] = cu_val[0] + delta;
  
  for (ki=0; ki<20; ki++)
  {
     s1221(pcurve,0,cu_val[1],&cu_left,pt->ecoef,&kstat);
     if (kstat < 0) goto error;      
     s1773(pt,psurf,aepsge,estart2,eend2,par_val,par_val,&kstat);
     if (kstat < 0) goto error;      
     s1421(psurf,1,par_val,&s_left1,&s_left2,s0,norm,&kstat);
     if (kstat < 0) goto error;
     for(kj=0; kj<dim; kj++) diff[kj] = s0[kj] - pt->ecoef[kj];
     new_y = s6norm(norm,dim,norm,&kstat);
     if (kstat == 0)
     {
	(*dist)=s6length(diff,dim,&kstat);
	new_cu_val = cu_val[1];
	goto out;
     }
     if (((*dist)=s6length(diff,dim,&kstat)) < aepsge)
     {
	new_cu_val = cu_val[1];
	goto out;
     }
     y[1] = s6scpr(norm,diff,dim);
     new_y = y[1]/y[0];
     if (new_y > 1.0000000000001)
     {
	if (shift)
	{
	   new_cu_val = cu_val[1];
	   goto out;
	}
	delta = -delta;
	/* ALA, UJK, sept 93, update cu_val[1]*/
	cu_val[1] = cu_val[0] + delta;
	shift++;	
     }
     else if (y[0]*y[1] <= 0.0 || fabs(new_y) < 0.5) break;
     else
     {
	if (cu_val[1]+delta <= aend1 && 
	    cu_val[1]+delta >= astart1) cu_val[1] += delta;
	else if (cu_val[1] < aend1)  	cu_val[1] = aend1;
	else if (cu_val[1] > astart1)   cu_val[1] = astart1;
	else 
	{
	   new_cu_val = cu_val[1];
	   goto out;
	}
     }
  }
  
  if (ki == 20)
  {
     *jstat = 2;
     goto out;
  }

  for (knbit=0; knbit < 50; knbit++)
  {
     delta_y = y[0]-y[1];
     if (fabs(delta_y) < REL_COMP_RES) break;
     
     new_cu_val = cu_val[1] + y[1]*(cu_val[1]-cu_val[0])/delta_y;
     if (new_cu_val >= aend1)
     {
	new_cu_val = aend1;
	if (cu_val[0] == aend1 || cu_val[1] == aend1) goto out;
     }
     else if (new_cu_val <= astart1)
     {
	new_cu_val = astart1;
	if (cu_val[0] == astart1 || cu_val[1] == astart1) goto out;
     }

     s1221(pcurve,0,new_cu_val,&cu_left,pt->ecoef,&kstat);
     if (kstat < 0) goto error;      
     s1773(pt,psurf,aepsge,estart2,eend2,par_val,par_val,&kstat);
     if (kstat < 0) goto error;      
     s1421(psurf,1,par_val,&s_left1,&s_left2,s0,norm,&kstat);
     if (kstat < 0) goto error;
     for(kj=0; kj<dim; kj++) diff[kj] = s0[kj] - pt->ecoef[kj];
     new_y = s6norm(norm,dim,norm,&kstat);
     if (kstat == 0)
     {
	(*dist) = s6length(diff,dim,&kstat);
	goto out;
     }
     if (((*dist)=s6length(diff,dim,&kstat)) < aepsge) goto out;
     new_y = s6scpr(norm,diff,dim);
     
     if ((y[0] < 0.0 && y[1] > 0.0) ||
	 (y[0] > 0.0 && y[1] < 0.0))
     {
	if ((new_y > 0.0 && y[0] > 0.0) ||
	    (new_y < 0.0 && y[0] < 0.0))
	{
	   cu_val[0] = new_cu_val;
	   y[0] = new_y;
	}
	else
	{
	   cu_val[1] = new_cu_val;
	   y[1] = new_y;
	}
     }
     else
     {
	if ( y[0] < 0.0 && new_y > 0.0)
	{
	   if (y[0] < y[1])
	   {
	      cu_val[0] = new_cu_val;
	      y[0] = new_y;
	   }
	   else
	   {
	      cu_val[1] = new_cu_val;
	      y[1] = new_y;
	   }
	}
	else if ( y[0] > 0.0 && new_y < 0.0)
	{
	   if (y[0] > y[1])
	   {
	      cu_val[0] = new_cu_val;
	      y[0] = new_y;
	   }
	   else
	   {
	      cu_val[1] = new_cu_val;
	      y[1] = new_y;
	   }
	}
	else if (y[0] > 0.0)
	{
	   if (y[0] > y[1])
	   {
	      if (new_y >=  y[0]) break;
	      cu_val[0] = new_cu_val;
	      y[0] = new_y;
	   }
	   else 
	   {
	      if (new_y >=  y[1]) break;
	      cu_val[1] = new_cu_val;
	      y[1] = new_y;
	   }
	     
	}
	else if (y[0] < 0.0)
	{
	   if (y[0] < y[1])
	   {
	      if (new_y <=  y[0]) break;
	      cu_val[0] = new_cu_val;
	      y[0] = new_y;
	   }
	   else 
	   {
	      if (new_y <=  y[1]) break;
	      cu_val[1] = new_cu_val;
	      y[1] = new_y;
	   }   
	}	   
     }
  }
  
  /* Iteration completed.  */
  
  goto out;
  
  /* Error in allocation */
  
  err101:
    *jstat = -101;
    s6err("s1772_s6sekant1",*jstat,kpos);
    goto out;                  
    
  /* Error in input. Conflicting dimensions.  */
  
  err106: 
    *jstat = -106;
    s6err("s1772_s6sekant1",*jstat,kpos);
    goto out;                  
  
  /* Error in lower level routine.  */
  
  error : 
    *jstat = kstat;
    s6err("s1772_s6sekant1",*jstat,kpos);
    goto out;                  
  
  out:
    par_val[2] = new_cu_val;
    if(pt) freePoint(pt);
}


#if defined(SISLNEEDPROTOTYPES)
   static
      int
s1772_s6local_pretop(double dist,double diff[],double normal[],
		     double f[],double f_t[],double f_tt[],
		     double s[],double s_u[],double s_v[],
		     double s_uu[],double s_uv[],double s_vv[],
		     int dim, int*jstat)
#else
static int s1772_s6local_pretop(dist,diff,normal,f,f_t,f_tt,s,s_u,s_v,s_uu,
			s_uv,s_vv,dim,jstat)
     double dist;
     double diff[];
     double normal[];
     double f[];
     double f_t[];
     double f_tt[];
     double s[];
     double s_u[];
     double s_v[];
     double s_uu[];
     double s_uv[];
     double s_vv[];
     int    dim;
     int    *jstat;
#endif 
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : To find if we have a minimum or a maximum or a turning
*             point situation. This function assume that it is a singular 
*	      situation.
*
*
*
* INPUT      : dist      - The lengt of the different vector.
*	       diff[]    - The differens vector.
*	       normal[]  - The normal vector on surface.
*              f    - Value in point on curve.
*              f_t  - First derevative of the curve.
*              f_tt - Second derevative of the curve.
*              s    - Value in point on surface.
*	       s_u  - Derevative of the surface in first par-dir.
*	       s_u  - Derevative of the surface in second par-dir.
*	       s_uu - Second derevative of the surface in first par-dir.
*	       s_uv - Cross derevative of the surface.
*	       s_vv - Second derevative of the surface in second par-dir.
*	       dim    - Dimension of space the curve/surface lie in.
*   
*   OUTPUT  :  return value - 	= -1: degenerated system.
*				=  0: Max or turning point.
*				=  1: Minimum position.
*
*              jstat        - Status variable.
*                       	< 0 : error.
*                       	= 0 : ok.
*                       	> 0 : warning.
*
*
*   METHOD  : Computing and interpretation of curvatures.
*
*   REFERENCES : 
*-
*   CALLS      :
*
*   WRITTEN BY : Arne Laksaa, SI, Oslo, des. 1992.
*
************************************************************************
*/
{
  int kstat = 0;	/* Status variable.				*/
  int ki;		/* Counter.					*/
  int return_val;	/* For return value.				*/
  double a1,a2,a3,a4;   /* Matrix.					*/
  double *S_u = SISL_NULL;	/* Normalized s_u.				*/
  double *S_v;		/* Normalized s_v.				*/
  double *S_uxS_v;	/* Cross between S_u and S_v.			*/
  double *s_d;		/* Second derevative in diriction f_t.		*/
  double *N;		/* Normalized normal.				*/
  double *d_uv;		/* Normalized direction vector in par-plane.	*/
  double local[17];	/* Local array for allocations.			*/
  
  *jstat = 0;
  
  if (s6ang(diff,normal,dim) > ANGULAR_TOLERANCE) goto warn1;
  
  /* Allocate local used memory and set value pointers.*/

  if (dim > 3)
  {
     S_u = newarray(5*dim+2,double);
     if (S_u == SISL_NULL) goto err101;
  }
  else
     S_u  = local;
  
  S_v     = S_u+dim;
  S_uxS_v = S_v+dim;
  s_d     = S_uxS_v+dim;
  N       = s_d+dim;
  d_uv    = N+dim;
  
  s6norm(s_u,dim,S_u,&kstat);
  if (kstat == 0)  goto warn1;
  s6norm(s_v,dim,S_v,&kstat);
  if (kstat == 0)  goto warn1;
  s6crss(S_u,S_v,S_uxS_v);
  a1 = s6scpr(S_u,S_v,dim);
  a2 = s6scpr(f_t,S_u,dim);
  a3 = s6scpr(f_t,S_v,dim);
  if ((a4 = s6scpr(S_uxS_v,S_uxS_v,dim)) < SINGULAR) goto warn1;
  
  d_uv[0] = (a2 - a1*a3)/a4;
  d_uv[1] = (a3 - a1*a2)/a4;
  s6norm(d_uv,2,d_uv,&kstat);
  if (kstat == 0)  goto warn1;
  
  a1 = d_uv[0]*d_uv[0];
  a2 = d_uv[1]*d_uv[1];
  a3 = 2*d_uv[0]*d_uv[1];
  
  for (ki=0; ki<dim; ki++)
     s_d[ki] = a1*s_uu[ki] + a3*s_uv[ki] + a2*s_vv[ki];
  
  for (ki=0; ki<dim; ki++)
     N[ki] = diff[ki]/dist;  
  
  a1 = s6scpr(N,f_tt,dim) - s6scpr(N,s_d,dim);

  return_val = a1 > 1.0e-10;
  goto out;

  /* Error in allocation */

  err101: 
    *jstat = -101;
    s6err("s1772_s6local_pretop",*jstat,0);
    return_val = 0;                 
    goto out;
  
  /* Degenerated system.  */

  warn1: 
    return_val = -1;
    goto out;

  out:
    if (S_u != local && S_u != SISL_NULL) freearray(S_u);
    return return_val;
}

                                    

                                        
