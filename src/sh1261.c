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

#define SH1261

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh1261_s9evalbez(double [],int,double,double,double,
			     double *,int *);
#else
static void sh1261_s9evalbez();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
      sh1261(SISLCurve *pcurve1,SISLCurve *pcurve2,double ecoef1[],
	     int ik1,double ecoef2[],int ik2,SISLCurve **rcrtanc,int *jstat)
#else
void sh1261(pcurve1,pcurve2,ecoef1,ik1,ecoef2,ik2,rcrtanc,jstat) 
   SISLCurve *pcurve1,*pcurve2,**rcrtanc; 
     double ecoef1[],ecoef2[]; 
     int ik1,ik2,*jstat; 
#endif     
/*
*********************************************************************
* 
* PURPOSE : Blend two curves to produce a new, higher order curve.
*           Curves to be blended and blending functions are given.  
* 
* 
* 
* INPUT   : pcurve1 - First curve to be blended.  
*           pcurve2 - Second curve to be blended.  
*           ecoef1  - Vertices of first blending function. The function 
*                     is a one-dimensional Bezier curve on the same 
*                     parameter interval as the input curves.  
*           ik1     - Order of first blending function.  
*           ecoef2  - Vertices of second blending function. The function 
*                     is a one-dimensional Bezier curve on the same 
*                      parameter interval as the input curves.  
*           ik2     - Order of second blending function. 
*
* 
* OUTPUT :  rcrtanc - Produced curve.  
*           jstat - status messages 
*                   > 0 : warning 
*                   = 0 : ok 
*                   < 0 : error 
* 
* 
* METHOD :  Compute order and knot vector of the curve to 
*           produce. Then the number of vertices of the output 
*           curve is given. Evaluate the expression 
*           (first blending function)*(first input curve) 
*           + (second blending function)*(second input curve)
*           in a number of points equal to the number of vertices. Interpolate
*           in these pointes to find the vertices of the cross output curve.
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : s1221  - Evaluate B-spline curve.       
*              s1363  - Pick parametrization of curve.  
*              s1891 - Interpolate.                    
*              s6takeunion - Union of two ordered vectors.   
*              newCurve    - Create new curve.           
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int kstat = 0;       /* Status variable.  */
  int ki,kj;           /* Counters.         */
  int kk;              /* Number increased in order.     */
  int kleft = 0;       /* Parameter used in s1221.  */
  int korder;          /* Order of cross tangent curve.  */
  int knunion;         /* Number of knots in union knot vector.           */
  int knvert;          /* Number of vertices of cross derivative curve.   */
  int kn;              /* Number of vertices of cross derivative curve.   */
  int kder = 0;        /* Number of derivatives of curve to compute.      */
  int kind = pcurve1->ikind; /* Kind of curve.                     */
  int kdim = pcurve1->idim;  /* Dimension of geometry space.       */
  int kcopy = 1;             /* Parameter to newCurve.             */
  int knlr = 0;              /* Parameter used in s1891.           */
  int knrc = 0;              /* Parameter used in s1891.           */
  int kopen = SISL_CRV_OPEN; /* Parameter used in s1891.           */
  int *lder = 0;             /* Parameter to the interpolation.    */
  double tpar;         /* Parameter value of interpolation point.  */
  double tmin1,tmax1;  /* Parameter interval of first derivative curve.   */
  double tmin2,tmax2;  /* Parameter interval of second derivative curve.  */
  double tval1,tval2;  /* Values of blending functions.   */
  double *stunion = SISL_NULL;  /* Union of knot vectors of derivative curves. */
  double *stcross = SISL_NULL;  /* Knot vector of cross tangent curve.         */
  double *spar = SISL_NULL;     /* Array containing parameter values of 
			      interpolation points.       */
  double *sbcoef = SISL_NULL;   /* Vertices of cross tangent curve.   */
  double *spoint = SISL_NULL;   /* Interpolation points.       */
  double *sder1 = SISL_NULL;
  double *sder2 = SISL_NULL;    /* Value of second derivative curve.  */

  /* Test input dimensions.  */

  if (pcurve2->idim != kdim) goto err106;
  
  /* Test if the parameter intervals of the derivative curves are equal.
     First pick parameter intervals.  */

  s1363(pcurve1,&tmin1,&tmax1,&kstat);
  if (kstat < 0) goto error;
  
  s1363(pcurve2,&tmin2,&tmax2,&kstat);
  if (kstat < 0) goto error;
  
  /* Test parameter intervals.  */

  if (tmin1 != tmin2 || tmax1 != tmax2) goto err121;
  
  /* Find order of cross tangent curve.  */

  korder = MAX(pcurve1->ik+ik1-1,pcurve2->ik+ik2-1);
  
  /* Find union of knot vectors of derivative vectors.  */

  s6takeunion(pcurve1->et,pcurve1->ik+pcurve1->in,pcurve2->et,
	  pcurve2->ik+pcurve2->in,&stunion,&knunion,&kstat);
  if (kstat < 0) goto error;
  
  /* Compute number of vertices of cross derivative curve.  */

  kk = MAX(ik1,ik2) - 1;
  knvert = knunion - korder + kk + kk;
  
  /* Allocate scratch for local arrays.  */

  if ((stcross = newarray(korder+knvert,DOUBLE)) == SISL_NULL) goto err101;
  if ((spar = newarray(knvert,DOUBLE)) == SISL_NULL) goto err101;
  if ((spoint = newarray(kdim*knvert,DOUBLE)) == SISL_NULL) goto err101;
  if ((sder1 = newarray(kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((sder2 = newarray(kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((lder = new0array(knvert,INT)) == SISL_NULL) goto err101;
  
  /* Produce knot vector of cross tangent curve.  */

  memcopy(stcross+kk,stunion,knunion,DOUBLE);
  for (ki=1; ki<=kk; ki++)
    {
      stcross[kk-ki] = stcross[kk];
      stcross[knunion+kk+ki-1] = stcross[knunion+kk-1];
    }
  
  /* Set up interpolation conditions.  */

  for (ki=0; ki<knvert; ki++)
    {
      /* Compute parameter value of coefficient nr ki.  */

      for (tpar=(double)0.0,kj=1; kj<korder; kj++) tpar += stcross[ki+kj];
      tpar /= (korder - 1);
      
      spar[ki] = tpar;
      
      /* Evaluate first blending function.  */

      sh1261_s9evalbez(ecoef1,ik1,tmin1,tmax1,tpar,&tval1,&kstat);
      if (kstat < 0) goto error;
      
      /* Evaluate second blending function.  */
  
		     sh1261_s9evalbez(ecoef2,ik2,tmin2,tmax2,tpar,&tval2,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate first derivative curve.  */

      s1221(pcurve1,kder,tpar,&kleft,sder1,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate second derivative curve.  */

      s1221(pcurve2,kder,tpar,&kleft,sder2,&kstat);
      if (kstat < 0) goto error;
      
      /* Compute interpolation point.  */

      for (kj=0; kj<kdim; kj++)
	spoint[ki*kdim+kj] = tval1*sder1[kj] + tval2*sder2[kj];
    }
  
  /* Perform interpolation.  */

  s1891(spar,spoint,kdim,knvert,1,lder,kopen,stcross,&sbcoef,&kn,korder,
	knlr,knrc,&kstat);
  if (kstat < 0) goto error;
  
  /* Create cross tangent curve.  */
    if ((*rcrtanc = newCurve(kn,korder,stcross,sbcoef,
			   kind,kdim,kcopy)) == SISL_NULL) goto err101;
    (*rcrtanc)->cuopen = MAX(pcurve1->cuopen,pcurve2->cuopen);
  
  /* Cross tangent curve produced. */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Conflicting dimensions.  */

  err106 :
    *jstat = -106;
  goto out;

  /* Error in input. Derivative curves are not represented on the same
     parameter interval.  */

  err121 :
    *jstat = -121;
  goto out;
  

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :     
    /* Free space occupied by local arrays.  */

    if (stunion != SISL_NULL) freearray(stunion);
  if (stcross != SISL_NULL) freearray(stcross);
  if (sbcoef != SISL_NULL) freearray(sbcoef);
  if (spar != SISL_NULL) freearray(spar);
  if (spoint != SISL_NULL) freearray(spoint);
  if (sder1 != SISL_NULL) freearray(sder1);
  if (sder2 != SISL_NULL) freearray(sder2);
  if (lder != SISL_NULL) freearray(lder);
  
  return;
}


#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1261_s9evalbez(double ecoef[],int iorder,double amin,
			     double amax,double apar,double *cvalue,int *jstat)
#else	       
static void sh1261_s9evalbez(ecoef,iorder,amin,amax,apar,cvalue,jstat)
     int iorder,*jstat;
     double ecoef[],amin,amax,apar,*cvalue;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate one-dimensional Bezier curve. Compute value
*              at given parameter value.
*
*
*
* INPUT      : ecoef      - Vertices of curve.
*              iorder     - Order of curve.
*              amin       - Start parameter value of parameter interval.
*              amax       - End parameter value of parameter interval.
*              apar       - Parameter value of point to evaluate.
*                       
*
* OUTPUT     : cvalue     - Value of curve.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int ki,kj;
  double *scoef = SISL_NULL;
  double tdist = amax - amin;
  
  /* Test input.  */

  if ( DEQUAL(tdist,DZERO)  )
    {
      *cvalue = ecoef[0];
      *jstat = 1;
      goto out;
    }
  
  /* Copy coeficient array to local array. First create local array.  */

  if ((scoef = newarray(iorder,DOUBLE)) == SISL_NULL) goto err101;
  memcopy(scoef,ecoef,iorder,DOUBLE);
  
  /* Evaluate Bezier curve by subdividing in the parameter value apar. */

  for (ki=1; ki<iorder; ki++)
    for (kj=0; kj<iorder-ki; kj++)
      scoef[kj] = (apar-amin)*scoef[kj+1]/tdist + (amax-apar)*scoef[kj]/tdist;
  
  /* Curve evaluated. */

  *cvalue = scoef[0];
  *jstat = 0;
  goto out;
  

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  out :
    /* Free space occupied by local array.  */

    if (scoef != SISL_NULL) freearray(scoef);

  return;
}

