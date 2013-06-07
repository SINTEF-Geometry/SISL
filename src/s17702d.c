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
 *
 *
 */
#define S1770_2D


#include "sislP.h"

#define SINGULAR 1.0e-16
#define copy2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]=(b)[ki]
#define copy3(a,b,c,d) for (ki=0;ki<(d);ki++) (a)[ki]=(b)[ki]=(c)[ki]
#define incr2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]+=(b)[ki]
#define decr2(a,b,c) for (ki=0;ki<(c);ki++) (a)[ki]-=(b)[ki]
#define set_order(a)  {if((a)==1) order=0; else order=1;}

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static
void
   s1770_2D_s9corr(double [],double[],double,double,double,double,int*);
static
void
   s1770_2D_s9dir(double *dist,double diff[],double delta[],
	    double c1[],double c1_t[],double c1_tt[],
	    double c2[],double c2_t[],double c2_tt[],
	    int dim, int second, double* det, int* jstat);
static
void
   s1770_2D_s6sekant1(SISLCurve *pcurve1,SISLCurve *pcurve2,
		 double  par_val[], double delta, double *dist, double aepsge,
		 double astart1,double astart2,double aend1,double aend2,
		 double c1[], double c2[], double norm[],
		 int *jstat);
static
int
   s1770_2D_s6local_pretop(double dist,double diff[],double normal[],
		     double c1[],double c1_t[],double c1_tt[],
		     double c2[],double c2_t[],double c2_tt[],
		     int dim, int*jstat);
#else
static void s1770_2D_s9corr();
static void s1770_2D_s9dir();
static void s1770_2D_s6sekant1();
static int  s1770_2D_s6local_pretop();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
   s1770_2D(SISLCurve *pcurve1,SISLCurve *pcurve2,double aepsge,
	   double astart1,double astart2,double aend1,double aend2,
	   double anext1,double anext2,double *cpos1,double *cpos2,int *jstat)
#else
void s1770_2D(pcurve1,pcurve2,aepsge,astart1,astart2,
	   aend1,aend2,anext1,anext2,cpos1,cpos2,jstat)
     SISLCurve  *pcurve1;
     SISLCurve  *pcurve2;
     double aepsge;
     double astart1;
     double astart2;
     double aend1;
     double aend2;
     double anext1;
     double anext2;
     double *cpos1;
     double *cpos2;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              two curves to find a closest point or an intersection point.
*
*
* INPUT      : pcurve1 - Pointer to the first curve in the intersection.
*              pcurve2 - Pointer to the second curve in the intersection.
*              aepsge  - Geometry resolution.
*              astart1 - Start parameter value of the first curve.
*              astart2 - Start parameter value of the second curve.
*              aend1   - End parameter value of the first curve.
*              aend2   - End parameter value of the second curve.
*              anext1  - Start parameter value of the iteration on
*                        the first curve.
*              anext2  - Start parameter value of the iteration on
*                        the second curve.
*
*
*
* OUTPUT     : cpos1   - Parameter value of of first curve in intersection
*                        point.
*              cpos2   - Parameter value of of second curve in intersection
*                        point.
*              jstat   - status messages
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter directions.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
* REWISED BY : Arne Laksaa, SINTEF-SI, Sep 1993.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int ki;
  int kleft1=0,kleft2=0;    /* Variables used in the evaluator.            */
  int dim;                  /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int p_dir;                /* Changing direction in par-space.            */
  int g_up,ng_up,g_dir;     /* Changing direction in geometric space.      */
  int order;		    /* Order of methode.			   */
  int sing = 0;		    /* Mark that singularity has ocured.	   */
  int keep_order = 0;
  int max_it = 20;          /* Maximum number of iterations.               */
  double delta[2];          /* Parameter interval of the curves.           */
  double dist;              /* Distance between the positions.             */
  double prev_dist;         /* Previous difference between the curves.     */
  double d[2];		    /* Clipped distances between old and new par.
			       value in the two parameter directions.      */
  double det;
  double c_d[2];	    /* Computed distances ....			   */
  double nc_d[2];	    /* New computed distances ....		   */
  double *c1=SISL_NULL;          /* Value  of first curve.			   */
  double *c1_t;		    /* First derivatiev of curve.		   */
  double *c1_tt;	    /* Second derivatiev of curve.		   */
  double *c2;               /* Value of second curve.   		   */
  double *c2_t;		    /* First derivative of second curve. 	   */
  double *c2_tt;       	    /* Second  derivative of second curve.	   */
  double *diff;             /* Difference between the curves.              */
  double *prev_diff;        /* Previous difference.			   */
  double *norm;		    /* Normal to the second curve.		   */
  double *norm_1;	    /* Normal to the first curve.		   */
  double par_val[2];        /* Parameter values                            */
  double local[48];
  int corr = 0, div2 = 0, quick = *jstat;

  if (quick) max_it = 10;  /* Reduce requirement on exactness. */

  /* Test input.  */

  if (pcurve1->idim != pcurve2->idim) goto err106;

  dim = pcurve1 -> idim;

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  delta[0] = pcurve1->et[pcurve1->in] - pcurve1->et[pcurve1->ik - 1];
  delta[1] = pcurve2->et[pcurve2->in] - pcurve2->et[pcurve2->ik - 1];

  /* Allocate local used memory */

    if (dim > 3)
  {
     c1 = newarray(10*dim,double);
     if (c1 == SISL_NULL) goto err101;
  }
  else
     c1 = local;

  c1_t  = c1 + dim;
  c1_tt = c1_t + dim;
  c2    = c1_tt + dim;
  c2_t  = c2 + dim;
  c2_tt = c2_t + dim;
  diff  = c2_tt + dim;
  prev_diff = diff + dim;
  norm      = prev_diff + dim;
  norm_1    = norm + dim;



  /* Initiate variables.  */


  par_val[0] = anext1;
  par_val[1] = anext2;

  for (ki=1; ki<3; ki++)
  {
     set_order(ki);

     /* Evaluate 0-2.st derivatives of curve 1 */

     if (par_val[0] == aend1)
	s1227(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
     else
	s1221(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
     if (kstat < 0) goto error;


     /* Evaluate 0-2.st derivatives of curve 2 */

     if (par_val[1] == aend2)
	s1227(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
     else
	s1221(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
     if (kstat < 0) goto error;


     /* Compute the distanse vector and value and the new step. */

     s1770_2D_s9dir(&dist,diff,c_d, c1,c1_t,c1_tt,
		    		    c2,c2_t,c2_tt,dim,order,&det,&kstat);
     if (kstat < 0) goto error;
     if (kstat == 1) 		/* Singular matrix. */
     {
	if (order == 1 && dist > aepsge) goto singular;
	else if (order == 1) goto not_singular;
     }
     else break;
  }

  /* Correct if we are not inside the parameter intervall. */

  d[0] = c_d[0];
  d[1] = c_d[1];
  norm_1[0] = -c1_t[1]; norm_1[1] = c1_t[0];
  norm[0]   = -c2_t[1]; norm[1]   = c2_t[0];
  g_up = (s6scpr(diff,norm,dim) >= DZERO) ? 1 : -1;
  g_up += ((s6scpr(diff,norm_1,dim) >= DZERO) ? 10 : -10);
  s1770_2D_s9corr(d,par_val,astart1,aend1,astart2,aend2,&corr);

  prev_dist = dist;
  prev_diff[0] = diff[0];
  prev_diff[1] = diff[1];

  /* Iterate to find the intersection point.  */

  for (knbit = 0; knbit < max_it; knbit++)
  {
     incr2(par_val,d,2);

     while (1)
     {
	/* Evaluate 0-2.st derivatives of curve */

	if (par_val[0] == aend1)
	   s1227(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
	else
	   s1221(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
	if (kstat < 0) goto error;

	if (par_val[1] == aend2)
	   s1227(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
	else
	   s1221(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
	if (kstat < 0) goto error;

	/* Compute the distanse vector and value and the new step. */

	s1770_2D_s9dir(&dist,diff,nc_d,c1,c1_t,c1_tt,c2,c2_t,c2_tt,
		    dim,order,&det,&kstat);
	if (kstat < 0) goto error;
	if (kstat == 1)             /* Singular matrix.  */
	{
	   sing++;
	   if (order == 1 && dist > aepsge) goto singular;
	   else if (order == 1) goto not_singular;
	   else set_order(2);               /* Change order to 2. */
	}
	else
	{
	   norm_1[0] = -c1_t[1]; norm_1[1] = c1_t[0];
	   norm[0]   = -c2_t[1]; norm[1]   = c2_t[0];

	   ng_up = (s6scpr(diff,norm,dim) >= DZERO) ? 1 : -1;
	   ng_up += ((s6scpr(diff,norm_1,dim) >= DZERO) ? 10 : -10);
	   g_dir = (ng_up+g_up != 0);			/* 0 if changed. */
	   p_dir = (c_d[0]*nc_d[0] >= DZERO &&
		    c_d[1]*nc_d[1] >= DZERO);		/* 0 if changed. */

	   if (!order && g_dir && (!p_dir || dist > 0.4*prev_dist)
							&& !keep_order)
	   {
	      if (!quick && div2) div2 = 0;
	      set_order(2);
	   }
	   else if (order && !g_dir)
	   {
	      if (sing && dist > aepsge) goto singular;
	      else if (sing) goto not_singular;
	      if (div2) div2 = 0;
	      set_order(1);
	   }
 	   else
	   {
              keep_order = 0;
	      if (sing) sing = 0;
	      break;
	   }
	}
     }

     if (corr)
	if (!(p_dir && g_dir)) corr = 0;

     if (dist < prev_dist || p_dir)
     {

	/* Corrigate if we are not inside the parameter interval. */

	g_up = ng_up;
	copy3(d,c_d,nc_d,2);
	s1770_2D_s9corr(d, par_val, astart1, aend1, astart2, aend2, &corr);
	prev_dist = dist;
	copy2(prev_diff,diff,dim);

	/* if (corr > 3) break; */

	if (corr > 2 ||
	    ((fabs(d[0]/MAX(par_val[0],delta[0])) <= REL_COMP_RES) &&
	     (fabs(d[1]/MAX(par_val[1],delta[1])) <= REL_COMP_RES))) break;
	if (div2) div2 = 0;

	     if (corr > 1 && order)
	     {
		keep_order = 1;
		set_order(1);
	     }
     }

     else if (corr > 2 ||
	      ((fabs(d[0]/MAX(par_val[0],delta[0])) <= REL_COMP_RES) &&
	       (fabs(d[1]/MAX(par_val[1],delta[1])) <= REL_COMP_RES))) break;
     else
     {
	/* Not converging, corrigate and try again.  */

	if (dist > prev_dist && div2 > 5) break;
	if (quick && dist > prev_dist && div2 > 3) break;
	div2++;
	decr2(par_val,d,2);
	d[0] /= (double)2; d[1] /= (double)2;

/*	printf("XXX %d, dist=%f, orden=%d\n",div2,dist,order); */
     }
  }

  /* Iteration stopped, test if point founds found is within resolution */


  if (dim == 2 && fabs(det)<0.1)
  {
    if (order < 1)
    {
      set_order(2);

      if (par_val[0] == aend1)
	s1227(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
      else
	s1221(pcurve1,1+order,par_val[0],&kleft1,c1,&kstat);
      if (kstat < 0) goto error;

      if (par_val[1] == aend2)
	s1227(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
      else
	s1221(pcurve2,1+order,par_val[1],&kleft2,c2,&kstat);
      if (kstat < 0) goto error;
    }
    goto singular;
  }

  goto not_singular;

  singular:

     /*  if (!quick && dist > aepsge) */
     if (!quick && dist > aepsge && dim == 2)
     {
	ki = s1770_2D_s6local_pretop(dist,diff,norm,c1,c1_t,c1_tt,
				  c2,c2_t,c2_tt,dim,&kstat);

	if (kstat < 0) goto error;
	if (ki == 0)
	{
	   s1770_2D_s6sekant1(pcurve1,pcurve2,par_val,c_d[0],&dist,aepsge,
			   astart1,astart2,aend1,aend2,c1,c2,norm,&kstat);
	   if (kstat < 0) goto error;

	}
     }

not_singular:

  if (dist <= aepsge)
  {
    *jstat = 1;
  }
  else
  {

     s6diff(c1,c2,dim,norm);
     if ((PIHALF-s6ang(c1_t,norm,dim)) < ANGULAR_TOLERANCE &&
         (PIHALF-s6ang(c2_t,norm,dim)) < ANGULAR_TOLERANCE)
	*jstat = 3;
     else
	*jstat = 2;
  }

  *cpos1 = par_val[0];
  *cpos2 = par_val[1];

  /* Iteration completed.  */


  goto out;

  /* Error in allocation */

 err101: *jstat = -101;
  s6err("s1770_2D",*jstat,kpos);
  goto out;

  /* Error in input. Conflicting dimensions.  */

 err106: *jstat = -106;
  s6err("s1770_2D",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1770_2D",*jstat,kpos);
  goto out;

 out:
    if (c1 != local && c1 != SISL_NULL) freearray(c1);

    return;
}

#if defined(SISLNEEDPROTOTYPES)
static
void
   s1770_2D_s9corr(double gd[],double acoef[],double astart1,double aend1,
	     double astart2,double aend2,int *corr)
#else
static void s1770_2D_s9corr(gd,acoef,astart1,aend1,astart2,aend2,corr)
     double gd[];
     double acoef[];
     double astart1;
     double aend1;
     double astart2;
     double aend2;
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
*              astart1 - The lower boorder in first direction.
*              aend1   - The higher boorder in first direction.
*              astart2 - The lower boorder in second direction.
*              aend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gdn   - Old and new step value.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, sep. 1993.
*
*********************************************************************
*/
{
  int lcorr = 0;
  if (acoef[0] + gd[0] < astart1)
    {
       gd[0] = astart1 - acoef[0];
       lcorr=1;
    }
  else if (acoef[0] + gd[0] > aend1)
    {
       gd[0] = aend1 - acoef[0];
       lcorr=1;
    }

  if (acoef[1] + gd[1] < astart2)
    {
       gd[1] = astart2 - acoef[1];
       lcorr=1;
    }
  else if (acoef[1] + gd[1] > aend2)
    {
       gd[1] = aend2 - acoef[1];
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
   s1770_2D_s9dir(double *dist,double diff[],double delta[],
	    double c1[],double c1_t[],double c1_tt[],
	    double c2[],double c2_t[],double c2_tt[],
	    int dim, int second, double* det,int* jstat)
#else
static void s1770_2D_s9dir(dist,diff,delta,c1,c1_t,c1_tt,c2,c2_t,c2_tt,
			dim,second,det,jstat)
     double *dist;
     double diff[];
     double delta[];
     double c1[];
     double c1_t[];
     double c1_tt[];
     double c2[];
     double c2_t[];
     double c2_tt[];
     int    dim;
     int    second;
     double* det;
     int*   jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the distance vector and value beetween
*	       a points on the first curve and a point on the second
*	       curve. And to compute a next step on both curves.
*	       This is equivalent to the nearest way to the
*	       parameter plan in the tangent plan from a point in the
*	       distance surface between two curves.
*
*
* INPUT      : c1    - Value in point on curve 1.
*              c1_t  - First derevative of the curve 1.
*              c1_tt - Second derevative of the curve 1.
*              c2    - Value in point on curve 2.
*              c2_t  - First derevative of the curve 2.
*              c2_tt - Second derevative of the curve 2.
*	       dim   - Dimension of space the curves lie in.
*              second - 1 if we have to use second order method
*			0 if only first order method.
*
*
* OUTPUT     : dist    - The lengt of the different vector.
*	       diff[]  - The differens vector.
*	       delta[] - Relative step parameter values towards intersection on
*                        the curve 1 delta[0] and the curve 2 delta[1].
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to each other.
*	       The differens vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the two
*	       points on the curve is "diff" and the two derevativ vectors
*	       are "der1" and "der2", and the two wanted parameter distance
*	       are "dt1" and "dt2", then we get the following system of
*	       equations:
*		 <dt1*der1+dist-dt2*der2,der1+dt1*der1_2> = 0
*		 <dt1*der1+dist-dt2*der2,der2+dt2*der2_2> = 0
*	       This is futher:
*
*	| -<dist,der1_2>-<der1,der1>   <der2,der1> |  | dt1 |   | <dist,der1> |
*	|                           		   |  |     | = |             |
*	| -<der2,der1>	 <der2,der2>-<dist,der2_2> |  | dt2 |   | <dist,der2> |
*
*	       The solution of this matrix equation is the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989.
* REWICED BY : Arne Laksaa, SINTEF-SI sep. 1993.
*
*********************************************************************
*/
{
  int kstat;			/* Local status variable. 		  */
  double a1,a2,a3;		/* The A matrix, diagonal and A12.	  */
  double b1,b2;			/* The B matrix, diagonal.	  	  */
  double A[4],mat[4];		/* Matrix in linear equation to be solved */
  double h[2];			/* Left side in the equation.		  */
  double x[2];			/* Left side in the equation.		  */
  double r[2];			/* Left side in the equation.		  */
  long double ss,aa,xx,bb;	/* For use in iterative improvement.      */
  int    piv[2];		/* Pivotation array                       */
  int k,k3,j;			/* Counters.				  */


  /* Computing the different vector */

  s6diff(c1,c2,dim,diff);

  /* Computing the length of the different vector. */

  *dist = s6length(diff,dim,&kstat);
  if (kstat<0) goto error;

  if (second || dim != 2)
  {
     a1 = s6scpr(c1_t,c1_t,dim);
     a2 = s6scpr(c2_t,c2_t,dim);
     a3 = s6scpr(c1_t,c2_t,dim);
  }

  if (second)
  {
     b1 = s6scpr(diff,c1_tt,dim);
     b2 = s6scpr(diff,c2_tt,dim);
  }
  else b1=b2=0.0;

  if (second || dim != 2)
  {
     mat[0] = -a1-b1;	mat[1] = a3;
     mat[2] = -a3;	mat[3] = a2-b2;

     h[0] =  s6scpr(diff,c1_t,dim);
     h[1] =  s6scpr(diff,c2_t,dim);
  }
  else
  {
     mat[0] = -c1_t[0];	mat[1] = c2_t[0];
     mat[2] = -c1_t[1];	mat[3] = c2_t[1];

     h[0] =  diff[0];
     h[1] =  diff[1];
  }

  for (k=0;k<4;k++) A[k]=mat[k];
  for (k=0;k<2;k++) x[k]=h[k];

  *det = A[0]*A[3]-A[1]*A[2];
  if (fabs(*det) < 1.0e-16)
  {
     *jstat = 1;
     goto out;
  }

  /* solve the linear 2x2 system */

  s6lufacp(mat,piv,2,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }

  s6lusolp(mat,x,piv,2,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }

  for (k=0;k<2;k++) delta[k] = x[k];


  for (k=k3=0; k<2; k++,k3+=2)
  {
     for (ss=0.0,j=0; j<2; j++)
     {
	aa = A[j+k3];
	xx = x[j];
	ss += aa*xx;
     }
     bb = h[k];
     ss = bb-ss;
     r[k] = ss;
  }
  s6lusolp(mat,r,piv,2,&kstat);
  if (kstat<0) goto error;
  if (kstat == 1)
  {
     *jstat = 1;
     goto out;
  }

  for (k=0;k<2;k++) delta[k] = x[k] + r[k];

  *jstat = 0;
  goto out;

  error :
    *jstat = kstat;
    s6err("s1770_2D_s9dir",*jstat,0);
    goto out;

  out:
    return;
}


#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    s1770_2D_s6sekant1(SISLCurve *pcurve1,SISLCurve *pcurve2,
		 double  par_val[], double delta, double *dist, double aepsge,
		 double astart1,double astart2,double aend1,double aend2,
		 double c1[], double c2[], double norm[],
		 int *jstat)
#else
static void s1770_2D_s6sekant1(pcurve1,pcurve2,par_val,delta,dist,aepsge,
			   astart1,astart2,aend1,aend2,c1,c2,norm,jstat)
     SISLCurve  *pcurve1;
     SISLCurve  *pcurve2;
     double  par_val[];
     double  delta;
     double  *dist;
     double aepsge;
     double astart1;
     double astart2;
     double aend1;
     double aend2;
     double c1[];
     double c2[];
     double norm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Sekant methode iteration on the distance function between
*              two curves to find a closest point
*              or an intersection point.
*
*
* INPUT      : pcurve1   - Pointer to the first curve in the intersection.
*              pcurve2   - Pointer to the second curve in the intersection.
*              delta     - Parameter distance on first curve beetveen start values.
*              aepsge    - Geometry resolution.
*              c1        - Array for use in evaluation.
*              c2        - Array for use in evaluation.
*              norm      - Array for use in evaluation.
* INPUT/
* OUTPUT     : par_val[] - Parameter value of the curves in
*                          intersection point.
*              dist      - Distance in space.
* OUTPUT     : jstat     - status messages
*                                = 3   : A minimum distanse found.
*                                = 2   : Nothing found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Sekant mothode in two parameter directions.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, sep 1993.
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
  int cu1_left = 0;	    /* Keep left knot information for evaluator.   */
  int cu2_left = 0;	    /* Keep left knot information for evaluator.   */
  int shift = 0;	    /* Mark that the diriction have been changed.  */

  *jstat = 0;

  /* Test input.  */

  if (pcurve1->idim != pcurve2->idim) goto err106;
  dim = pcurve1->idim;
  diff = c1 + dim;

  if ((pt = newPoint(c1,dim,0)) == SISL_NULL) goto err101;

  if (delta == 0.0) delta =1e-15;

  if (par_val[0] < astart1) par_val[0] = astart1;
  else if (par_val[0] > aend1) par_val[0] = aend1;

  if ((par_val[0] == astart1 && delta < 0.0) ||
      (par_val[0] == aend1   && delta > 0.0))
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


  cu_val[0] = par_val[0];
  s1221(pcurve1,0,cu_val[0],&cu1_left,pt->ecoef,&kstat);
  if (kstat < 0) goto error;
  s1771(pt,pcurve2,aepsge,astart2,aend2,par_val[1],par_val+1,&kstat);
  if (kstat < 0) goto error;
  s1221(pcurve2,1,par_val[1],&cu2_left,c2,&kstat);
  if (kstat < 0) goto error;
  norm[0] = -c2[3]; norm[1] = c2[2];
  for(kj=0; kj<dim; kj++) diff[kj] = c2[kj] - pt->ecoef[kj];
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
  if (cu_val[1] < astart1) cu_val[1] = astart1;
  else if (cu_val[1] > aend1) cu_val[1] = aend1;

  for (ki=0; ki<20; ki++)
  {
    s1221(pcurve1,0,cu_val[1],&cu1_left,pt->ecoef,&kstat);
    if (kstat < 0) goto error;
    s1771(pt,pcurve2,aepsge,astart2,aend2,par_val[1],par_val+1,&kstat);
    if (kstat < 0) goto error;
    s1221(pcurve2,1,par_val[1],&cu2_left,c2,&kstat);
    if (kstat < 0) goto error;
    norm[0] = -c2[3]; norm[1] = c2[2];
    for(kj=0; kj<dim; kj++) diff[kj] = c2[kj] - pt->ecoef[kj];
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
      cu_val[1] = cu_val[0] + delta;
      if (cu_val[1] < astart1) cu_val[1] = astart1;
      else if (cu_val[1] > aend1) cu_val[1] = aend1;
      shift++;
    }
    else if (y[0]*y[1] <= 0.0 || new_y < 0.6) break;
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
     new_cu_val = par_val[0];
     goto out;
  }

  for (knbit=0; knbit < 25; knbit++)
  {
     delta_y = y[0]-y[1];
     if (fabs(delta_y) < REL_COMP_RES) break;

     if (y[0]*y[1] < 0.0 &&
	 (fabs(y[0]) < 6*fabs(y[1]) || fabs(y[1]) < 6*fabs(y[0])))
       new_cu_val = 0.5*(cu_val[1]+cu_val[0]);
     else
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

     s1221(pcurve1,0,new_cu_val,&cu1_left,pt->ecoef,&kstat);
     if (kstat < 0) goto error;
     s1771(pt,pcurve2,aepsge,astart2,aend2,par_val[1],par_val+1,&kstat);
     if (kstat < 0) goto error;
     s1221(pcurve2,1,par_val[1],&cu2_left,c2,&kstat);
     if (kstat < 0) goto error;
     for(kj=0; kj<dim; kj++) diff[kj] = c2[kj] - pt->ecoef[kj];
     norm[0] = -c2[3]; norm[1] = c2[2];
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
    s6err("s1770_2D_s6sekant1",*jstat,kpos);
    goto out;

  /* Error in input. Conflicting dimensions.  */

  err106:
    *jstat = -106;
    s6err("s1770_2D_s6sekant1",*jstat,kpos);
    goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1770_2D_s6sekant1",*jstat,kpos);
    goto out;

  out:
    par_val[0] = new_cu_val;
    if(pt) freePoint(pt);
}


#if defined(SISLNEEDPROTOTYPES)
   static
      int
	    s1770_2D_s6local_pretop(double dist,double diff[],double normal[],
		     double c1[],double c1_t[],double c1_tt[],
		     double c2[],double c2_t[],double c2_tt[],
		     int dim, int*jstat)
#else
static int s1770_2D_s6local_pretop(dist,diff,normal,c1,c1_t,c1_tt,c2,c2_t,c2_tt,
				dim,jstat)
     double dist;
     double diff[];
     double normal[];
     double c1[];
     double c1_t[];
     double c1_tt[];
     double c2[];
     double c2_t[];
     double c2_tt[];
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
* INPUT      : dist     - The lengt of the different vector.
*	       diff[]   - The differens vector.
*	       normal[] - The normal vector of the second curve.
*              c1    	- Value in point on curve.
*              c1_t  	- First derevative of the curve.
*              c1_tt 	- Second derevative of the curve.
*              c2    	- Value in point on curve.
*              c2_t  	- First derivative of the curve.
*              c2_tt 	- Second derivative of the curve.
*	       dim   	- Dimension of space the curves lie in. dim = 2.
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
*   METHOD  : Testing size and direction of the second derivatives
*             of the curve.
*
*   REFERENCES :
*-
*   CALLS      :
*
*   WRITTEN BY : Arne Laksaa and Vibeke Skytt, SINTEF, Oslo, sep. 1993.
*
************************************************************************
*/
{
  int kpos = 0;
  int return_val;	 /* For return value.				*/
  double l_1,l_2;
  double v_1,v_2;
  double k_1,k_2;
  double r_1,r_2;


  *jstat = 0;

  if (dim != 2) goto err101;

  l_1 = s6scpr(c1_tt,diff,dim);
  l_2 = s6scpr(c2_tt,diff,dim);

  if (( l_1 < 0.0 && l_2 > 0.0) || (l_1 > 0.0 && l_2 < 0.0))
  {
    return_val = 1;
    goto out;
  }

  v_1 = s6scpr(c1_t,c1_t,dim);
  v_1 = v_1*sqrt(v_1);
  k_1 = fabs(c1_t[0]*c1_tt[1] - c1_tt[0]*c1_t[1]);
  if (k_1 < REL_COMP_RES)	r_1 = 0.0;
  else				r_1 = v_1/k_1;

  v_2 = s6scpr(c2_t,c2_t,dim);
  v_2 = v_2*sqrt(v_2);
  k_2 = fabs(c2_t[0]*c2_tt[1] - c2_tt[0]*c2_t[1]);
  if (k_2 < REL_COMP_RES)	r_2 = 0.0;
  else				r_2 = v_2/k_2;


  if (( l_1 < 0.0 || l_2 < 0.0) && (r_1 > r_2 + dist))
    return_val = 1;
  else if (( l_1 > 0.0 || l_2 > 0.0) && (r_2 > r_1 + dist))
    return_val = 1;
  else
    return_val = 0;

  goto out;


    /* Error in allocation */

  err101:
    *jstat = -101;
    s6err("s1770_2D_s6local_pretop",*jstat,kpos);
    goto out;

  out:
    return return_val;
}
