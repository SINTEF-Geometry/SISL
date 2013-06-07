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
 * $Id: s1605.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1605

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s1605(SISLCurve *pc,double aepsge,double **gpoint,int *jnbpnt,int *jstat)
#else
void s1605(pc,aepsge,gpoint,jnbpnt,jstat)
     SISLCurve  *pc;
     double aepsge;
     double **gpoint;
     int    *jnbpnt;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate a set of points on a B-spline curve. 
*              The straight lines between the points will not deviate
*              more than aepsge from the B-spline curve at any point.
*             
* INPUT      : pc     - The input B-spline curve.   
*              aepsge - Geometry resolution, maximum distance allowed between
*                       the curve and the straight lines to be calculated. 
*
* OUTPUT     : gpoint - Calculated points
*              jnbpnt - No. of calculated points  
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : The error in approximating a given function f on an interval
*              [a,b] by linear spline approximation is bounded by
*                      (b-a)**2*max(abs(D2(f)))/8
*              where max(abs(D2(f))) denotes the maximum of the absolute value
*              of the 2. derivative of f on the interval [a,b].
*              In order to approximate a derivative curve, we get the bound
*                      (b-a)**2*max(abs(D(2)(f))/8
*              For a spline curve this estimate applies on each of the 
*              component of the geometry space.
*              To determine the sampling points, we start from the left of the
*              given spline and consider one knot interval at the time. On each
*              knot interval the 2. derivative is a spline of order 
*              ik-2 (ik is the order of f), the coefficients of which 
*              is 2. order divided differences of the original coefficients. 
*              The div. diff. of largest absolute value then gives an estimate 
*              of the maximum of the absolute value of  the 2. derivative 
*              so that the error in linear interpolation becomes less than the 
*              tolerance. This can be done in each of the geometry components 
*              and to be on the safe side we choose the smallest of the 
*              resulting step lengths. This step length is then adjusted so 
*              that the sampling points can be distributed uniformly in 
*              the interval, and always with one sampling point at 
*              the beginning of the interval and one at the end, i.e. 
*              all original knots become sampling points.
*              At break points, two sampling points will be placed.
*
* USE        : Array gpoint is allocated inside this routine and must 
*              be made free after use.
*
* REFERENCES :
*
*-                                                 
* CALLS      : s1221,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 10. Nov 1988
*
* REWRITTEN BY : Bjoern Olav Hoset, SI, Oslo, Norway, feb. 1992.
*                The code is identical to s1613 which Arne Laksaa 
*                wrote,.
*********************************************************************
*/
{
  int kstat = 0;           /* Status variable.  */
  int kpos  = 0;           /* Position of error reported. */
  int ki,kj,kk,kh,kr,kl;   /* Counters.         */
  int kikr;                /* Order of current derivative curve. */
  int kih;                 /* Current location in sdd2.          */
  int kihn;                /* Next location in sdd2.             */
  int kjh;                 /* Index in sdd2.                     */
  int kstop;               /* Number of divided differences to compute.   */
  int kdim = pc->idim;     /* Dimension of space in which the curve lies. */
  int kncoef = pc->in;     /* Number of vertices of curve.  */
  int korder = pc->ik;     /* Order of curve.               */
  int kordnew = korder-2;  /* Order of derivative curve of which 
			      to find an upper bound.         */
  int knmbel = 100;        /* Number of elements with which the output
			      array is to be increased.                   */       
  int kmaxpar = knmbel;    /* Length of output parameter array.           */
  int kpar = 0;            /* Number of parameter values computed.        */
  int kant;                /* Number of sampling points at a knot interval. */
  int left = 0;            /* Help index to evaluator.                    */
  double tant;             /* Number of sampling points at a knot interval. */
  double tfac;             /* Factor used in taking divided differences.    */
  double ta;               /* Start value of parameter interval.          */
  double tb;               /* End value of parameter interval.            */
  double th;               /* Distance between output parameter values.   */
  double *st;              /* Pointer to knot vector of curve.            */
  double *par = SISL_NULL;      /* Array used to store parameter values.       */
  double *sh  = SISL_NULL;      /* Work array.    */
  double *sh1 = SISL_NULL;      /* Work array.    */
  double *sdd = SISL_NULL;      /* Work array used to compute divided differences. */
  double *sdd2 = SISL_NULL;     /* Array used to store final divided differences.  */
  double *smaxd = SISL_NULL;    /* Array used to store maximum divided differences.*/
  
  /* Test input.  */
  
  if (korder < 1) goto err110;
  if (kncoef < korder) goto err111;
  if (kdim < 1) goto err102;
  
  if (korder == 1)
    {
      kpar = kncoef;
      if ((*gpoint = newarray(kpar*kdim,DOUBLE)) == SISL_NULL) goto err101;
      memcopy(gpoint,pc->ecoef,kpar*kdim,DOUBLE);
      *jnbpnt = kpar;
      *jstat = 0;
      goto out;
    }
  
  
  /* Set local pointer to knot vector of curve.  */
  
  st = pc->et;
  
  /* Allocate some scratch for output array.  */
  
  if ((par = newarray(knmbel,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Allocate scratch for internal arrays.  */
  
  if ((sdd = new0array((kordnew+6)*kdim,DOUBLE)) == SISL_NULL) goto err101;
  sdd2 = sdd+3*kdim;
  smaxd = sdd2+kordnew*kdim;
  sh = smaxd+kdim;
  sh1 = sh+kdim;
  
  /* Main loop to determine the interpolation points in each knot interval
     by computing 2. order differences and taking the maximum of the 
     appropriate kordnew 2. order differences on each knot interval. The
     number of uniform interplation points in the knot interval can then
     by estimated using the bound described in METHOD.
     In the following the index ki will point at the next coefficient to be 
     included in the computation of the 2. order derivatives.
     kpar points to the location in the array par where the parameter value
     of the next interpolation point is to be stored.
     kih points to the location in sdd2 where the next 2. order difference
     is to be stored.        */
  
  kpar=0;
  tb = st[korder-1];
  
  for (kih=0, ki=0; ki<kncoef; ki=kj)
    {
      /* Compute the index of the first knot greater than tb.  */
      
      for (kj=MAX(korder-1,ki)+1; st[kj]==tb; kj++);
      ta = tb;
      tb = st[kj];
      
      /* Check if the multiplisity of the knot is greater than 
	 korder-2. */
      
      if (kj-ki >= korder-1)
	{
	  if (kpar+1 > kmaxpar)
	    if ((par = increasearray(par,(kmaxpar+=knmbel),DOUBLE)) 
		== SISL_NULL) goto err101;
	  par[kpar] = ta;
	  kpar++;
	}
      
      /* If the 2. order differences may be different from zero, 
	 Compute them.   */
      
      for (kk=ki; korder>2 && kk<kj; kk++)
	{
	  /* Pick the next coefficient.  */
	  
	  memcopy(sh,pc->ecoef+kk*kdim,kdim,DOUBLE);
	  
	  /* Compute the 1, 2 difference at kk.  */
	  
	  kikr = korder - 1;
	  kstop = MIN(2,korder-kj+kk);
	  for (kr=0; kr<kstop; kr++)
	    {
	      tfac = (double)kikr/(st[kk+kikr] - ta);
	      kikr--;
	      for (kh=0; kh<kdim; kh++)
		sh1[kh] = (sh[kh] - sdd[kr*kdim+kh])*tfac;
	      memcopy(sdd+kr*kdim,sh,kdim,DOUBLE);
	      memcopy(sh,sh1,kdim,DOUBLE);
	    }
	  memcopy(sdd+kr*kdim,sh,kdim,DOUBLE);
	  
	  /* Compute the maximum of the 2. order difference. The
	     kordnew previous 2. order differences are stored in
	     sdd2 and the current one (sdd+3*kdim) is to overwrite
	     sdd2+kih*kdim. The fact that only one new element enters
	     sdd2 at a time can be taken advantage of in computing 
	     the maximum difference.   */
	  
	  if (kstop == 2)
	    {
	      for (kihn=(kih+1)%kordnew, kh=0; kh<kdim; kh++)
		{
		  sh[kh] = fabs(sh[kh]);
		  if (sdd2[kih*kdim+kh] < smaxd[kh])
		    smaxd[kh] = MAX(smaxd[kh],sh[kh]);
		  else if (sh[kh] >= smaxd[kh])
		    smaxd[kh] = sh[kh];
		  else
		    {
		      for (kjh=kihn, smaxd[kh]=sh[kh], kl=0;
			   kl<kordnew-1; kl++, kjh=(kjh+1)%kordnew)
			smaxd[kh] = MAX(smaxd[kh],sdd2[kjh*kdim+kh]);
		    }
		  sdd2[kih*kdim+kh] = sh[kh];
		}
	      kih = kihn;
	    }
	}
      
      /* Compute the number of interpolation points.  */
      
      for (kant=0, kh=0; kh<kdim; kh++)
	{
	  tant = (tb-ta)*sqrt(smaxd[kh]/((double)8.0*aepsge));
	  kant = MAX(MAX(kant,(int)tant),1);
	}
      
      /* Make sure that par is great enough.  */
      
      if (kpar+kant >= kmaxpar)
	if ((par = increasearray(par,(kmaxpar+=MAX(kant,knmbel)),DOUBLE)) 
	    == SISL_NULL) goto err101;
      
      /* Compute the parameter values of the interpolation points.  */
      
      for (th=(tb-ta)/(double)(kant+1), kk=0; kk<kant; kpar++,kk++)
	par[kpar] = ta + (double)(kk+1)*th;
      par[kpar] = tb;
      kpar++;
    }
  
  /* Make the points.  */
  
  if ((*gpoint = newarray(kpar*kdim,DOUBLE)) == SISL_NULL) goto err101;
  
  for (kh=kk=0; kk<kpar; kk++,kh+=kdim)
    {
      s1221(pc,0,par[kk],&left,*gpoint+kh,&kstat);
      if (kstat < 0) goto err101;
    }
  
  /* Task performed.  */
  
  *jnbpnt = kpar;
  *jstat = 0;
  goto out;
  
  
  /* Error in scratch allocation.  */
  
  err101 :
    *jstat = -101;
    s6err("s1605",*jstat,kpos);
    goto out;
  
  /* Error in input. Dimension less than one.  */
  
  err102 :
    *jstat = -102;
    s6err("s1605",*jstat,kpos);
    goto out;
  
  /* Error in input. Order less than one.  */
  
  err110 :
    *jstat = -110;
    s6err("s1605",*jstat,kpos);
    goto out;
  
  /* Error in input. Number of coefficients less than order.  */
  
  err111 :
    *jstat = -111;
    s6err("s1605",*jstat,kpos);
    goto out;
  
  out :
    /* Free scratch occupied by local arrays.  */
    
    if (sdd != SISL_NULL) freearray(sdd);
    if (par != SISL_NULL) freearray(par);
}


