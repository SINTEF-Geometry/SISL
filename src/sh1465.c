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
 * $Id: sh1465.c,v 1.2 2001-03-19 15:59:04 afr Exp $
 *
 */


#define SH1465

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
   static void sh1465_s9der2(double [],double [],double [],double [],
			     double [],int,int,double [],int *);
   typedef void (*fshapeProc)(double [],double [],int,int,int *);
#else
   static void sh1465_s9der2();
   typedef void (*fshapeProc)();
#endif  

#if defined(SISLNEEDPROTOTYPES)
void
  sh1465(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
	 double etwist[],double etang[],double eder[],int *jstat)
#else
void sh1465(fshape,vboundc,icurv,etwist,etang,eder,jstat)
     fshapeProc fshape;
     double etwist[],etang[],eder[];
     SISLCurve *vboundc[];
     int icurv,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a vertex region whith an equal number of sides, 
*              evaluate the first blending surface in the corner lying in 
*              the middle of the vertex region. Compute the tangent vectors 
*              in the middle vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 2*icurv.
*              icurv   - Number of sides.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is icurv*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s6norm   - Normalize vector.  
*              s6scpr   - Scalar product between two vectors.  
*              s6length - Lenght of vector.  
*              s6dist   - Distance between two vectors. 
*              s6crss   - Cross product between two vectors. 
*              s6diff   - Difference between two vectors.  
*              s6curvrad - Estimate curvature radius in point. 
*              s1325    - Calculate tangent length given opening
*                         angle of curve segment and curvature radius.
*              s1221    - Evaluate curve.  
*              sh1465_s9der2 - Compute 2. derivatives of the first 
*			       blending patch.    
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  int kstat = 0;      /* Status variable.  */
  int kder = 0;       /* Number of derivatives of curve to evaluate. */
  int kder1 = 1;      /* Number of derivatives of curve to evaluate. */
  int ki,kj;          /* Counters.  */
  int kdim = 3;       /* Dimension of geometry space.  */
  int kleft = 0;      /* Parameter used in curve evaluation.  */
  int kcurv2;         /* Number of edges diveded in two.  */
  double tpar;        /* Parameter value at which to evaluate curve. */
  double t1;          /* Scalar product between tangent in midpoint
			 and normal of vertex region in midpoint.    */
  double trad1,trad2; /* Estimate of curvature radius of curves in
			 endpoints.   */
  double ta1;         /* Opening angle of curve segment.  */
  double tb1,tb2;     /* Length of original tangents in endpoints
			 of curve segments.               */
  double tang1,tang2; /* Length of tangents based on curvature radius. */
  double tscal1;      /* Scaler product between unit tangent vectors.  */
  double tdist;       /* Distance between endpoints of curve segment.  */
  double smid[3];     /* Midpoint of vertex region.       */
  double snorm1[3];   /* Cross product of two tangents in the midpoint. */
  double snorm[3];    /* Normal of vertex region in the midpoint.      */
  double svec[6];     /* Tangent along curve in the midpoint of the two
			 first position curves.           */
  double *sder = SISL_NULL;   /* Value of boundary curves at the midpoints of
			    the curves.                   */
  double *stang = SISL_NULL;  /* Tangent vectors in the midpoint of the region. */
  SISLCurve *qc;      /* Local pointer to edge curve.  */
  
  kcurv2 = icurv/2;

  /* Allocate scratch for values on edge curves. */

  if ((sder = newarray(2*icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((stang = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;

  for (ki=0; ki<icurv; ki++)
    {  
      qc = vboundc[2*ki];

      /* Test dimension of curve. */

      if (qc->idim != kdim) goto err102;

      /* Evaluate boundary curves of current edge in the midpoint.
	 First find midpoint.  */

      tpar = (*(qc->et + qc->ik - 1) + *(qc->et + qc->in))/(double)2.0;
      
      /* Evaluate position curve. */

      s1221(qc,kder1,tpar,&kleft,sder+2*ki*kdim,&kstat);
      if (kstat < 0) goto error;
      
      if (ki < 2) memcopy(svec+ki*kdim,sder+(2*ki+1)*kdim,kdim,DOUBLE);
      if (ki == 1) kder1 = 0;
      
      /* Evaluate cross derivative curve. */

      s1221(vboundc[2*ki+1],kder,tpar,&kleft,sder+(2*ki+1)*kdim,&kstat);
      if (kstat < 0) goto error;
    }
  
  for (ki=0; ki<kcurv2; ki++)
    {
       /* Set new tangent length based on curvature radius. First estimate
         curvature radius.  */

      s6curvrad(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,sder+(2*ki+1)*kdim,
		kdim,&trad1,&kstat);
      s6curvrad(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,sder+(2*(ki+kcurv2)+1)*kdim,
		kdim,&trad2,&kstat);
  
      /* Normalize tangents. */

      tb1 = s6norm(sder+(2*ki+1)*kdim,kdim,sder+(2*ki+1)*kdim,&kstat);
      tb2 = s6norm(sder+(2*(ki+kcurv2)+1)*kdim,kdim,sder+(2*(ki+kcurv2)+1)*kdim,&kstat);
      
      /* Compute distance between endpoints of curve.  */

      tdist = s6dist(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,kdim);
      
      /* Find opening angle of curve segment.  */

      tscal1  = s6scpr(sder+(2*ki+1)*kdim,sder+(2*(ki+kcurv2)+1)*kdim,kdim);

      if (tscal1 >= DZERO)
	tscal1  = MIN((double)1.0,tscal1);
      else
	tscal1  = MAX((double)-1.0,tscal1);

      ta1 = acos(tscal1);

      if (fabs(ta1) < ANGULAR_TOLERANCE) ta1 = DZERO;

      if (DNEQUAL(ta1,DZERO))
	{
	  /*  Make tangents based on radius of curvature */

	  tang1 = s1325(trad1,ta1);
	  tang2 = s1325(trad2,ta1);
	}

      /* Test if the found tangent length can be used. Otherwise
	 adjust the length.   */

      if (DEQUAL(ta1,DZERO) || trad1 < 0) tang1 = tdist/(double)3.0;
      if (DEQUAL(ta1,DZERO) || trad2 < 0) tang2 = tdist/(double)3.0;
      if (tang1 > (double)0.5*tdist || tang2 > (double)0.5*tdist) 
	{
	  tang1 = tb1;
	  tang2 = tb2;
	}
      
      /* Set tangent length of tangents in endpoints of the curves to find. */

      for (kj=0; kj<kdim; kj++)
	{
	  sder[(2*ki+1)*kdim+kj] *= tang1;
	  sder[(2*(ki+kcurv2)+1)*kdim+kj] *= tang2;
	} 
    }  
  

  /* Estimate midpoint of region and tangent in midpoint.  */

  for (kj=0; kj<kdim; kj++) 
    {
      smid[kj] = (double)0.0;
      snorm[kj] = (double)0.0;
      
      for (ki=0; ki<kcurv2; ki++)
	{
	  /* Estimate midpoint.  */

	  smid[kj] += (sder[2*ki*kdim+kj] + sder[2*(ki+kcurv2)*kdim+kj])/(double)2.0
	    + (sder[(2*ki+1)*kdim+kj] 
	       + sder[2*(ki+kcurv2)*kdim+kdim+kj])/(double)8.0;
	  
	  /* Compute tangent.  */

	  stang[(ki+kcurv2)*kdim+kj] = (double)1.5*(sder[2*(ki+kcurv2)*kdim+kj]
					- sder[2*ki*kdim+kj])
	    + (double)0.25*(sder[2*(ki+kcurv2)*kdim+kdim+kj]
			   - sder[(2*ki+1)*kdim+kj]);
	  stang[ki*kdim+kj] = - stang[(ki+kcurv2)*kdim+kj];
	}
      smid[kj] /= (double)kcurv2;
    }

  /* Find medium plane given by the tangents.  */

  for (ki=0; ki<kcurv2; ki++)
    {
      s6crss(stang+ki*kdim,stang+(ki+1)*kdim,snorm1);
  
      for (kj=0; kj<kdim; kj++) 
	snorm[kj] += snorm1[kj]/(double)kcurv2;
    }
  
  (void)s6norm(snorm,kdim,snorm,&kstat);
  
  /* Project the tangents into this plane.  */

  for (ki=0; ki<icurv; ki++)
    {
      t1 = -s6scpr(stang+ki*kdim,snorm,kdim);
      for (kj=0; kj<kdim; kj++) stang[ki*kdim+kj] += t1*snorm[kj];
    }
  
  /* Application driven routine to alter the midpoint and tangents in the
     midpoint.  */

  fshape(smid,stang,kdim,icurv,&kstat);
  if (kstat < 0) goto error;
  
  /* Compute second order derivatives of first surface in the midpoint.  */

  sh1465_s9der2(sder,smid,stang,snorm,svec,icurv,kdim,eder+3*kdim,&kstat);
  if (kstat < 0) goto error;
  
  /* Set lengths of tangents and second derivatives according to patches
     whith side length equal to half the curvelengths considered in this
     routine.  */

  for (ki=0; ki<icurv*kdim; ki++) stang[ki] *= (double)0.5;
  for (ki=3*kdim; ki<6*kdim; ki++) eder[ki] *= (double)0.25;
  
  /* Copy position and tangent information into the output array giving
     derivatives of the first blending surface in the midpoint. */

  memcopy(eder,smid,kdim,DOUBLE);
  memcopy(eder+kdim,stang,2*kdim,DOUBLE);
  
  /* Copy tangents into output array containing tangents.  */

  memcopy(etang,stang+kdim,(icurv-1)*kdim,DOUBLE);
  memcopy(etang+(icurv-1)*kdim,stang,kdim,DOUBLE);
					
  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */

  err102 :
    *jstat = -102;
  goto out;
  
  /* Error in a lower level function.  */

 error:
  *jstat = kstat;
  goto out;
  
  out :

    /* Free space occupied by local arrays.  */

    if (sder) freearray(sder);
  if (stang) freearray(stang);
  
  return;
}


#if defined(SISLNEEDPROTOTYPES)
static void
  sh1465_s9der2(double ebound[],double epoint[],double etang[],
		double enorm[],double evec[],int icurv,
		int idim,double eder2[],int *jstat)
#else	       
static void sh1465_s9der2(ebound,epoint,etang,enorm,evec,icurv,
			  idim,eder2,jstat)
     int       icurv,idim,*jstat;
     double    ebound[],epoint[],etang[],enorm[],evec[],eder2[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a vertex region with an equal number of sides, 
*              estimate the 2. derivatives of the first blending surface 
*              in the midpoint of the region.
*
*
*
* INPUT      : ebound  - Position of boundary curve and cross tangent at the 
*                        midpoint of each edge. Dimension is 2*icurv*idim.
*              epoint  - The midpoint of the vertex region. Dimension is idim.
*              etang   - Tangents in epoint, pointing towards the midpoints of
*                        the edges. Dimension is icurv*idim.
*              enorm   - Normal to surface in midpoint of region. Dimension is idim.
*              evec    - The tangent at the midpoint of the boundary curves at the
*                        two first edges. Dimension is 2*idim.
*              icurv   - Number of sides. icurv is an equal number.
*              idim    - Dimension of the geometry space.
*                       
*
* OUTPUT     : eder2   - Second derivative of the first blending surface in 
*                        the corner at the midpoint. The sequence is the
*                        following : 2. derivative in the 1. parameter direction, 
*                        mixed derivative and 2. derivative in the 2. parameter 
*                        direction. Dimension is 3*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221     - Curve evaluator.  
*              s1334     - Curve interpolation. 
*              s6lufacp  - LU-factorizing of matrix. 
*              s6lusolp  - Solve to LU-factorized equation system. 
*              s6curvature - Compute curvature vector of curve. 
*              s6scpr    - Scalar product between two vectors. 
*              s6length  - Length of vector.  
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Status variable.  */
  int ki,kj;         /* Counters.         */
  int kleft = 0;     /* Parameter to curve evaluator.  */
  int knbpnt = 6;    /* Number of interpolation conditions. */
  int kcnsta = 0;    /* No extra condition on startpoint in interpolation. */
  int kcnend = 0;    /* No extra condition on endpoint in interpolation. */
  int kopen = 1;     /* Produce open curve.   */
  int kcurv2;        /* Number of edges divided into 2. */
  int kord = 6;      /* Order of interpolated curve.    */
  int knpar;         /* Number of different parameter values of 
			interpolation conditions.       */
  int kder = 2;      /* Number of derivatives to evaluate. */
  int ll[3];         /* Pivoting array in solving equation system. */
  double tstpar = (double)0.0;   /* Start parameter of interpolated curve. */
  double tendpar;    /* End parameter of interpolated curve. */
  double tpar;       /* Parameter value at which to evaluate. */
  double tncurv;     /* Normal curvature at midpoint of curve. */
  double te,tf,tg;   /* Coefficients of first fundamental form of surface. */
  double tl,tm,tn;   /* Coefficients of second fundamental form of surface. */
  double talfa;      /* Angle between the two parameter directions of a
			rectangular patch in the parameter area of the 
			n-sided vertex region.       */
  double tcos;       /* Cosinus of the angle talfa.  */
  double tdudt;      /* Factor in parameter change.  */
  double tdvdt;      /* Factor in parameter change.  */
  double tform1;     /* First fundamental form.      */
  double spoint[18]; /* Interpolation conditions.    */
  double stype[6];   /* Type of interpolation conditions. */
  double *spar = SISL_NULL;  /* Parameter value of interpolation conditions. */
  double sder[18];   /* Value and derivatives of curve in the midpoint. */
  double scurv[3];   /* Curvature vector in midpoint of curve. */
  double smat[9];    /* Matrix in equation system to compute mixed derivative. */
  SISLCurve *qc = SISL_NULL;  /* Curve across the vertex region between midpoint
			       of edge position curves through the midpoint
			       of the region.      */
  
  kcurv2 = icurv/2;
  
  /* Test input.  */

  if (idim != 3) goto err104;
  
  /* Make curves that limits first blending surface. */

  for (ki=0; ki<2; ki++)
    {
      /* Set up interpolation conditions of curve. */

      for (kj=0; kj<idim; kj++)
	{
	  spoint[kj] = ebound[2*(kcurv2+ki)*idim+kj];
	  spoint[idim+kj] = ebound[(2*(kcurv2+ki)+1)*idim+kj];
	  spoint[2*idim+kj] = epoint[kj];
	  spoint[3*idim+kj] = etang[ki*idim+kj];
	  spoint[4*idim+kj] = ebound[2*ki*idim+kj];
	  spoint[5*idim+kj] = -ebound[(2*ki+1)*idim+kj];
	}
      
      /* Type of interpolation conditions.  */

      stype[0] = (double)1.0;
      stype[1] = (double)4.0;
      stype[2] = (double)1.0;
      stype[3] = (double)4.0;
      stype[4] = (double)1.0;
      stype[5] = (double)4.0;

      /* Interpolate curve.  */

      s1334(spoint,knbpnt,idim,stype,kcnsta,kcnend,kopen,kord,tstpar,&tendpar,
	    &qc,&spar,&knpar,&kstat);
      if (kstat < 0) goto error;
      
      /* Evaluate curve in midpoint.  */
				
      tpar = spar[1];
      s1221(qc,kder,tpar,&kleft,sder+3*ki*idim,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Copy 2. derivatives of the two curves to the 2. derivatives of the
     surface in the midpoint in the 1. and 2. parameter direction. */

  memcopy(eder2,sder+2*idim,idim,DOUBLE);
  memcopy(eder2+2*idim,sder+5*idim,idim,DOUBLE);
  
  /* Compute curvature vector in the midpoint of the first curve.  */

  s6curvature(sder,idim,scurv,&kstat);

  /* Compute normal curvature of the first curve at the midpoint of the region. */
  
  tncurv = s6scpr(scurv,enorm,idim);

  /* Compute parameter direction of the curve compared to that of the
     first blending surface. */

  talfa = TWOPI/(double)icurv;
  tcos = cos(talfa);
  tdudt = (DEQUAL(tcos+(double)1.0,(double)1.0)) ? (double)0.0 : (double)1.0/tcos;
  tdvdt = (double)1.0;
  
  /* Compute coefficients of the first fundamental form of the surface. */

  te = s6scpr(sder+idim,sder+idim,idim);
  tf = s6scpr(sder+idim,sder+4*idim,idim);
  tg = s6scpr(sder+4*idim,sder+4*idim,idim);
  
  /* Compute the first fundamental form.  */

  tform1 = te*tdudt*tdudt + (double)2.0*tf*tdudt*tdvdt + tg*tdvdt*tdvdt;
  
  /* Compute 1. and 3. coefficient of the second fundamental form of the surface. */

  tl = s6scpr(sder+2*idim,enorm,idim);
  tn = s6scpr(sder+5*idim,enorm,idim);
  
  /* Compute 2. coefficient of the second fundamental form which is equal to
     the length of the component of the twist along the surface normal.  */

  eder2[idim] = tm = (tncurv*tform1 - tl*tdudt*tdudt - tn*tdvdt*tdvdt)/(double)2.0;
  
  /* Set the length of the component of the twist along the derivative in
     the first parameter direction equal to zero.   */

  eder2[idim+1] = (double)0.0;
  
  /* Set the length of the component of the twist along the derivative 
     in the second parameter direction equal to zero.   */

  eder2[idim+2] = (double)0.0;
  
  /* Compute twist vector at the midpoint.  */

  memcopy(smat,enorm,idim,DOUBLE);
  memcopy(smat+idim,sder+idim,idim,DOUBLE);
  memcopy(smat+2*idim,sder+4*idim,idim,DOUBLE);
  
  s6lufacp(smat,ll,3,&kstat);
  if (kstat < 0) goto error;
  
  s6lusolp(smat,eder2+idim,ll,3,&kstat);
  if (kstat < 0) goto error;
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */

  err104 :
    *jstat = -104;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :
    return;
}     


