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
 * $Id: sh1461.c,v 1.2 2001-03-19 15:59:03 afr Exp $
 *
 */


#define SH1461

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
typedef void (*fshapeProc)(double [],double [],int,int,int *);
typedef void (*initProc)(fshapeProc,SISLCurve *[],int,double [],
			 double [],double [],int *);
#else
typedef void (*fshapeProc)();
typedef void (*initProc)();
#endif

#if defined(SISLNEEDPROTOTYPES)
static void sh1461_s9coef(double [],double [],double [],double [],int,
			  double *,double *,double *,double *,int *);
static void sh1461_s9hermit(double [],int,int,int *);
static void sh1461_s9chcoor(double [],int,double [],int,double [],
			    int,double [],double [],double [],int,
			    double [],int *,int *);
static void sh1461_s9mult(double [],double [],int,int,double [],int *);
static void sh1461_s9comder(int,int,double [],int,double,double,double,
			    double,double [],int *);
static double sh1461_s9ang(double [],double [],int);
#else
static void sh1461_s9coef();
static void sh1461_s9hermit(); 
static void sh1461_s9chcoor();
static void sh1461_s9mult(); 
static void sh1461_s9comder(); 
static double sh1461_s9ang();   
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh1461(fshapeProc fshape,initProc f_initmid,
	 SISLCurve *vboundc[],int icurv,SISLSurf *vsurf[],int *jstat)
#else	 
void sh1461(fshape,f_initmid,vboundc,icurv,vsurf,jstat)
     fshapeProc  fshape;
     initProc    f_initmid;
     SISLCurve   *vboundc[];
     SISLSurf    *vsurf[];
     int         icurv;
     int         *jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Make a blend consisting of n tensor-product surfaces over
*              a vertex region with n edges, n>=3. Along each edge position
*              and cross-tangent conditions are given. The blend is supposed
*              to be G1.
*
*
*
* INPUT      : fshape    - Application driven routine that gives the user an
*                          ability to change the middle point of the region
*                          (the vertex at which the blending surfaces meet),
*                          and the tangent vectors in the middle point along
*                          the curves which divedes the region. 
*              f_initmid - Function used to compute the position and 
*                          derivatives of the first blending surface in
*                          the middle of the vertex region. This function
*                          is dependant on the number of edges of the region.
*              vboundc   - Pointers to curves describing edge-conditions.
*                          For each edge, two curves are given, one describing
*                          position and one cross-derivatives. The edges are
*                          sorted counter-clockwise around the region, and the
*                          curves are oriented counter-clockwise. The array has
*                          dimension 2*icurv, i.e. 6.
*              icurv     - Number of boundary curves. icurv >= 3.
*                       
*
* OUTPUT     : vsurf      - Pointers to blending surfaces. The array has dimension
*                           icurv.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Hahn's method. Split the region in n 4-sided regions. Define
*              position and cross tangent conditions along the inner edges
*              such that G1-continuity is satisfied. Define the blending
*              patches over the 4-sided regions as Coons patches.
*
*
* REFERENCES : Joerg Hahn : "Filling Polygonal Holes with Rectangular Patches"
*              Theory and Practice of Geometric Modeling, Blackburn Oct 1988.
*              
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221      - Evaluate curve.
*              s1706      - Turn orientation of curve.
*              s1401 - Rectangular blending.
*              newCurve   - Create new curve.
*              freeCurve  - Free space occupied by a curve.
*              sh1461_s9coef - Compute coefficients used for computing 
*                              derivatives in the midpoint.  
*              sh1461_s9hermit - Compute vertices of Bezier curve.  
*              sh1461_s9chcoor - Perform coordinate change on derivative curve.  
*              sh1461_s9mult   - Multiply two 4-order Bezier curves.  
*              sh1461_s9comder - Compute given 2. derivative of patch in midpoint
*			         of region.
*              sh1461_s9ang    - The angle between two vectors.  
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 03.90.
*
*********************************************************************
*/
{
  int kstat = 0;                /* Status variable.                */
  int kdim = 3;                 /* Dimension of geometry space.    */
  int ki,kj,kk,kh;              /* Counters.                       */
  int kder = 2;                 /* Number of derivatives to
                                   compute while evaluating curve. */
  int knmb = 0;                 /* Number of derivatives computed
                                   in midpoint of patch.           */
  int kleft = 0;                /* Parameter of curve evaluation.  */
  int kkcrt2;                   /* Order of cross tangent curve.   */
  int lder[4];
  double tpar;                  /* Parameter value at which to 
                                   evaluate curve.                 */
  double tro00,tro01,tro10,tro11;  /* Coefficients used to compute
                                      derivatives at midpoint of region. */
  double *sder = SISL_NULL;          /* Array of derivatives at midpoint    */
  double *stwist = SISL_NULL;        /* Twist vectors of corners of region. */
  double *stang = SISL_NULL;         /* Tangent vectors at midpoint.        */
  double spos[15];              /* Conditions and vertices of position
                                   curve along inner edge of region.   */
  double scrt1[21];             /* Conditions and vertices of cross tangent
                                   curve along inner edge of region.   */
  double scrt2[12];             /* Conditions and vertices of cross tangent
                                   curve along inner edge of region.   */
  double stpos[10];             /* Knot vector corresponding to spos.  */
  double stcrt1[14];            /* Knot vector corresponding to scrt1. */
  double stcrt2[8];             /* Knot vector corresponding to scrt2. */
  double sblend[4];             /* Vertices of blending function.      */
  double sdum[6];               /* Value and 1. derivative of position
                                   boundary curve.                     */
  double sdum2[6];              /* Value and 1. derivative of cross
                                   tangent boundary curve.             */
  double *stwist2 = SISL_NULL;       /* Twist vectors of 4-sided region.    */
  double *st,*st2;              /* Pointers into knot vectors. Used to
                                   change parametrization of curve.    */
  double tstart;                /* Start parameter of knot vector.     */
  SISLCurve *qc;                /* Pointer to curve.                   */
  SISLCurve **qbound = SISL_NULL;    /* Boundary conditions of 4-sided regions. */
  
  for (ki=0; ki<4; ki++) lder[ki] = 2;
  
  for (ki=0; ki<=kder; ki++) knmb += ki + 1;
  
  /* Initiate knot vectors of position curve and cross tangent
     curves of inner edges.  */

  for (ki=0; ki<5; ki++)
    {
      stpos[ki] = (double)0.0;
      stpos[5+ki] = (double)1.0;
    }
  for (ki=0; ki<4; ki++)
    {
      stcrt2[ki] = (double)0.0;
      stcrt2[4+ki] = (double)1.0;
    }
  for (ki=0; ki<7; ki++)
    {
      stcrt1[ki] = (double)0.0;
      stcrt1[7+ki] = (double)1.0;
    }
      
      
  /* Allocate scratch for local arrays.  */

  if ((sder = new0array(icurv*kdim*knmb,DOUBLE)) == SISL_NULL) goto err101;
  if ((stwist = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((stwist2 = newarray(4*icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((stang = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((qbound = newarray(8*icurv,SISLCurve*)) == SISL_NULL) goto err101;
  
  for (ki=0; ki<icurv; ki++)
    {
      /* Test dimension of curve.   */

      if (vboundc[ki]->idim != kdim) goto err102;
      
      /* Initiate twist vector of n-sided region.  */

      tpar = (double)0.0;
      s1221(vboundc[2*ki+1],1,tpar,&kleft,sdum,&kstat);
      if (kstat < 0) goto error;
      
      memcopy(stwist+ki*kdim,sdum+kdim,kdim,double);
    }
 
  /* Initiate those twist vectors of the rectangular blending surfaces
     that coincide with the twist vectors of the n-sided region.  */

  for (ki=0; ki<kdim; ki++)
    for (kj=0; kj<icurv; kj++)
      {
	kk = (kj + 1) % icurv;
	stwist2[(kj*4+2)*kdim+ki] = stwist[kk*kdim+ki]/(double)4.0;
      }
  
      
  /* Set up blending function.  */

  sblend[0] = (double)1.0;
  sblend[1] = sblend[2] = sblend[3] = (double)0.0;
  sh1461_s9hermit(sblend,4,1,&kstat);
  
  /* Set up value and derivative of the 1. blending surface in the
     midpoint of the vertex region.  */

  f_initmid(fshape,vboundc,icurv,stwist,stang,sder,&kstat);
  if (kstat < 0) goto error;
  
  for (ki=1; ki<icurv; ki++)
    {
      /* Compute value and derivatives of the other blending surfaces
	 in the midpoint of the vertex region.  */

      /* Copy the value in the midpoint.  */

      memcopy(sder+ki*knmb*kdim,sder,kdim,DOUBLE);
      
      /* Copy the first derivatives of the current patch.  */

      memcopy(sder+(ki*knmb+1)*kdim,stang+(ki-1)*kdim,kdim,DOUBLE);
      memcopy(sder+(ki*knmb+2)*kdim,stang+ki*kdim,kdim,DOUBLE);
      
      /* The second derivative in the first parameter direction 
	 is equal to the first derivative in the second parameter
	 direction of the previous surface.    */

      memcopy(sder+(ki*knmb+3)*kdim,sder+((ki-1)*knmb+5)*kdim,kdim,DOUBLE);
      
      /* Compute help variables used to define the rest of the derivatives. */

      sh1461_s9coef(stang+(icurv-1)*kdim,stang,stang+(ki-1)*kdim,
	     stang+ki*kdim,kdim,&tro00,&tro01,&tro10,&tro11,&kstat);
      if (kstat < 0) goto error;
      
      /* Compute mixed derivative of the patch.  */

		     sh1461_s9comder(1,1,sder+3*kdim,kdim,tro00,tro01,tro10,tro11,
	       sder+(ki*knmb+4)*kdim,&kstat);
	  
      /* Compute the second derivative in the second parameter direction. */

		     sh1461_s9comder(0,2,sder+3*kdim,kdim,tro00,tro01,tro10,tro11,
	       sder+(ki*knmb+5)*kdim,&kstat);
    }  

  /* Set up conditions for 4-edged blending.  */

  for (ki=0; ki<icurv; ki++)
    {
      /* Copy value and derivatives of the current patch of the midpoint of
	 the vertex region into arrays containing interpolation conditions 
	 and twist vector of current blending surface. */

      memcopy(spos,sder+ki*knmb*kdim,kdim,DOUBLE);
      memcopy(spos+kdim,sder+(ki*knmb+2)*kdim,kdim,DOUBLE);
      memcopy(spos+2*kdim,sder+(ki*knmb+5)*kdim,kdim,DOUBLE);
      memcopy(scrt2,sder+(ki*knmb+1)*kdim,kdim,DOUBLE);
      memcopy(scrt2+kdim,sder+(ki*knmb+4)*kdim,kdim,DOUBLE);
      memcopy(stwist2+4*ki*kdim,sder+(ki*knmb+4)*kdim,kdim,DOUBLE);
      
      /* Compute value and derivatives of the current inner edge curve at 
	 the boundary of the vertex region.  */

      kj = (ki < icurv-1) ? ki+1 : 0;
      
      qc = vboundc[2*kj];
      tpar = (*(qc->et + qc->ik - 1) + *(qc->et + qc->in))/(double)2.0;
      
      s1221(qc,1,tpar,&kleft,sdum,&kstat);
      if (kstat < 0) goto error;
      
      s1221(vboundc[2*kj+1],1,tpar,&kleft,sdum2,&kstat);
      if (kstat < 0) goto error;
      
      /* Copy value and derivatives at the current edge into arrays containing 
	 interpolation conditions and twist vector of current blending surface. */

      memcopy(spos+4*kdim,sdum,kdim,DOUBLE);
      for (kh=0; kh<kdim; kh++)
	{
	  spos[3*kdim+kh] = sdum2[kh]/(double)2.0;
	  scrt2[3*kdim+kh] = -sdum[kdim+kh]/(double)2.0;
	  stwist2[(ki*4+3)*kdim+kh] = scrt2[2*kdim+kh] = -sdum2[kdim+kh]/(double)4.0;
	  stwist2[(kj*4+1)*kdim+kh] = -stwist2[(ki*4+3)*kdim+kh];
	}
      
      /* Perform Hermit interpolation of position and first cross tangent curve. */

      sh1461_s9hermit(spos,5,kdim,&kstat);
      if (kstat < 0) goto error;
      
		     sh1461_s9hermit(scrt2,4,kdim,&kstat);
      if (kstat < 0) goto error;
      
      /* Construct second cross tangent curve. */

      sh1461_s9chcoor(sblend,4,spos,5,scrt2,4,sder+(knmb*ki+1)*kdim,sder+(knmb*ki+2)*kdim,
	       sder+(kj*knmb+2)*kdim,kdim,scrt1,&kkcrt2,&kstat);
      if (kstat < 0) goto error;
      
      /* Represent inner boundary conditions as curves. */

      if ((qbound[8*kj] = newCurve(5,5,stpos,spos,1,kdim,1)) == SISL_NULL) 
	goto err101;
      if ((qbound[8*kj+1] = newCurve(kkcrt2,kkcrt2,stcrt1,scrt1,1,kdim,1)) == SISL_NULL)
	goto err101;
      if ((qbound[8*ki+6] = newCurve(5,5,stpos,spos,1,kdim,1)) == SISL_NULL) 
	goto err101;
      if ((qbound[8*ki+7] = newCurve(4,4,stcrt2,scrt2,1,kdim,1)) == SISL_NULL)
	goto err101;
      
      /* Split current boundary of the vertex region.  */

      s1710(qc,tpar,qbound+8*ki+4,qbound+8*kj+2,&kstat);
      if (kstat < 0) goto error;
      
      s1710(vboundc[2*kj+1],tpar,qbound+8*ki+5,qbound+8*kj+3,&kstat);
      if (kstat < 0) goto error;

      /* Turn direction of curves corresponding to standard edge 3.  */

      s1706(qbound[8*ki+4]);
      s1706(qbound[8*ki+5]);

      /* Adjust length of derivative curves. */

      for (kk=0; kk<kdim*(qbound[8*ki+5]->in); kk++)
	*(qbound[8*ki+5]->ecoef+kk) *= (double)0.5;
      
      for (kk=0; kk<kdim*(qbound[8*kj+3]->in); kk++)
	*(qbound[8*kj+3]->ecoef+kk) *= (double)0.5;

      /* Reparametrize boundary curves in order to get correct tangent length.  */

      for (kh=2; kh<6; kh++)
	{
	  kk = (kh > 3) ? 8*ki : 8*kj;
	  kk += kh;
	  
	  for (st=qbound[kk]->et,tstart=st[0],
	       st2=st+qbound[kk]->in+qbound[kk]->ik;
	       st<st2; st++)
	    *st = (double)2.0*(*st - tstart);
	}
      
    }
  
  for (ki=0; ki<icurv; ki++)
    {
      /* Perform rectangular blending.  */

       s1401(qbound+8*ki,stwist2+4*ki*kdim,vsurf+ki,&kstat);  
       /* s1390(qbound+8*ki,vsurf+ki,lder,&kstat); */
      if (kstat < 0) goto error;
    }
  
  /* Blending performed.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
/* Dimension not equal to 3.  */

  err102 :
    *jstat = -102;
  goto out;
  
  /* Error in lower level function.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :

    /* Free space occupied by local arrays.  */

    if (sder) freearray(sder);
  if (stwist) freearray(stwist);
  if (stwist2) freearray(stwist2);
  if (stang) freearray(stang);
  if (qbound)
    {
      for (ki=0; ki<8*icurv; ki++)
	if (qbound[ki]) freeCurve(qbound[ki]);
      freearray(qbound);      
    }
  
    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1461_s9hermit(double econd[],int icond,int idim,int *jstat)
#else	       
static void sh1461_s9hermit(econd,icond,idim,jstat)
     int icond,idim,*jstat;
     double econd[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Hermite interpolation of position and icond-3 derivatives
*              in one endpoint and position and derivative in the other
*              endpoint, represented as a Bezier curve on the interval [0,1].
*
*
*
* INPUT      : icond      - Number of interpolation conditions. 
*                           icond = 4 or icond = 5.
*              idim       - Dimension of geometry space.
*
*
* INPUT/OUTPUT : econd    - Interpolation conditions as input, Bezier coefficients
*                           as output. The dimension is icond*idim.
*                       
*
* OUTPUT     : jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
  int ki;    /* Index.  */

  /* Test input. The number of conditions has to be 4 or 5.  */

  if (icond != 4 && icond != 5) goto err001;
  
  if (icond == 4)
    {
      /* Hermit interpolation with Bezier curve of order 4.  */

      for (ki=0; ki<idim; ki++)
	{
	  econd[idim+ki] = ONE_THIRD*econd[idim+ki] + econd[ki];
	  econd[2*idim+ki] = ONE_THIRD*econd[2*idim+ki] + econd[3*idim+ki];
	}
    }
  
  if (icond == 5)
    {
      /* Hermit interpolation with Bezier curve of order 5.  */

      for (ki=0; ki<idim; ki++)
	{
	  econd[idim+ki] = ONE_FOURTH*econd[idim+ki] + econd[ki];
	  econd[2*idim+ki] = econd[2*idim+ki]/(double)12.0 
	    + (double)2.0*econd[idim+ki] - econd[ki];
	  econd[3*idim+ki] = ONE_FOURTH*econd[3*idim+ki] + econd[4*idim+ki];	  
	}
    }
   
  *jstat = 0;
  goto out;
  
  /* Error in input. Number of coefficients not 4 or 5.  */

  err001 :
    *jstat = -1;
  goto out;
  
  out :
    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1461_s9coef(double evec1[],double evec2[],double evec3[],
		double evec4[],int idim,double *cro00,double *cro01,
		double *cro10,double *cro11,int *jstat)
#else
static void sh1461_s9coef(evec1,evec2,evec3,evec4,idim,cro00,cro01,
			  cro10,cro11,jstat)
     int idim,*jstat;
     double evec1[],evec2[],evec3[],evec4[];
     double *cro00,*cro01,*cro10,*cro11;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute factors of expression to find derivatives of 
*              patches in the midpoint of the vertex region.
*
*
*
* INPUT      : evec1   - Derivative in first parameter direction of first patch.
*              evec2   - Derivative in second parameter direction of first patch.
*              evec3   - Derivative in first parameter direction of current patch.
*              evec4   - Derivative in second parameter direction of current patch.
*                        NB! All these derivative vectors are expected to lie in
*                            the same plane.
*              idim    - Dimension of geometry space.
*
*
* OUTPUT     : cro00   - First factor.
*              cro01   - Second factor.
*              cro10   - Third factor.
*              cro11   - Fourth factor.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* CALLS      : s6length  - Length of vector.
*              s6scpr    - Scalar product between two vectors.
*              s6crss    - Cross product between two vectors.
*
*********************************************************************
*/
{
  int kstat = 0;                       /* Status variable.              */
  int ksin,ksin1,ksin2,ksin3,ksin4;    /* Sign of sinus of angles.      */
  double tang,tang1,tang2,tang3,tang4; /* Angles between input vectors. */
  double tl1,tl2,tl3,tl4;              /* Lengths of input vectors.     */
  double tsin,tsin1,tsin2,tsin3,tsin4; /* Sinus of angles.              */
  double snorm[3];                     /* Normal of plane spanned by
                                          tangent vectors.              */
  double svec[3];                      /* Vector in tangent plane normal
                                          to a specific tangent vector. */
  double tlvec;                        /* Length of the vector svec.    */
  
  /* Compute the normal to the tangent plane spanned by the input vectors.  */

  s6crss(evec1,evec2,snorm);
  
  /* Compute the lengths of the vectors evec1 to evec4. */

  tl1 = s6length(evec1,idim,&kstat);
  tl2 = s6length(evec2,idim,&kstat);
  tl3 = s6length(evec3,idim,&kstat);
  tl4 = s6length(evec4,idim,&kstat);
  
  /* Compute the vector lying in the tangent plane normal to evec1. */

  s6crss(snorm,evec1,svec);
  tlvec = s6length(svec,idim,&kstat);
  
  /* Compute the sinus of angles between evec1 and the other input vectors. */

  tsin = s6scpr(svec,evec2,idim)/(tlvec*tl2);
  tsin1 = s6scpr(svec,evec3,idim)/(tlvec*tl3);
  tsin2 = s6scpr(svec,evec4,idim)/(tlvec*tl4);
  
  /* Compute the vector lying in the tangent plane normal to evec2. */

  s6crss(snorm,evec2,svec);
  tlvec = s6length(svec,idim,&kstat);

  /* Compute the sinus of angles between evec2 and the later input vectors. */

  tsin3 = s6scpr(svec,evec3,idim)/(tlvec*tl3);
  tsin4 = s6scpr(svec,evec4,idim)/(tlvec*tl4); 
  
  /* Fetch the sign of the sinuses of the angles.  */

  ksin = (tsin < DZERO) ? -1 : 1;
  ksin1 = (tsin1 < DZERO) ? -1 : 1;
  ksin2 = (tsin2 < DZERO) ? -1 : 1;
  ksin3 = (tsin3 < DZERO) ? -1 : 1;
  ksin4 = (tsin4 < DZERO) ? -1 : 1;
  
  /* Compute the angles in a more stable way. */

  tang = sh1461_s9ang(evec1,evec2,idim);
  tang1 = sh1461_s9ang(evec1,evec3,idim);
  tang2 = sh1461_s9ang(evec1,evec4,idim);
  tang3 = sh1461_s9ang(evec2,evec3,idim);
  tang4 = sh1461_s9ang(evec2,evec4,idim);
  
  /* Compute the sinuses of the angles.  */

  tsin = ksin*sin(tang);
  tsin1 = ksin1*sin(tang1);
  tsin2 = ksin2*sin(tang2);
  tsin3 = ksin3*sin(tang3);
  tsin4 = ksin4*sin(tang4);
  
  /* Compute the output factors.  */

  *cro00 = (tl3*tsin1)/(tl2*tsin);
  *cro01 = (tl4*tsin2)/(tl2*tsin);
  *cro10 = -(tl3*tsin3)/(tl1*tsin);
  *cro11 = -(tl4*tsin4)/(tl1*tsin);
  
  *jstat = 0;
  
  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1461_s9chcoor(double eblend[],int iordblend,double epos[],
		  int iordpos,double ecrt1[],int iordcrt1,
		  double evec1[],double evec2[],double evec3[],
		  int idim,double ecrt2[],int *jordcrt2,int *jstat)
#else	       
static void sh1461_s9chcoor(eblend,iordblend,epos,iordpos,ecrt1,iordcrt1,
			    evec1,evec2,evec3,idim,ecrt2,jordcrt2,jstat)
     int iordblend,iordpos,iordcrt1,idim,*jordcrt2,*jstat;
     double evec1[],evec2[],evec3[];
     double eblend[],epos[],ecrt1[],ecrt2[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute that cross tangent curve belonging to one of the 
*              blending patches, that is a mapping of the derivative
*              curves of the previous patch.
*
*
*
* INPUT      : eblend    - Vertices of blending function. The dimension
*                          of the array is iordblend.
*              iordblend - Number of vertices of blending function. 
*                          iordblend = 4.
*              epos      - Vertices of position curve. This curve is
*                          between the corners (0,0) and (0,1) of the
*                          previous patch, and the corners (0,0) and 
*                          (1,0) of the current patch.
*              iordpos   - Number of vertices of the position curve.
*                          iordpos = 5.
*              ecrt1     - Vertices of input cross tangent curve. This
*                          curve is the derivative curve in the 1. parameter
*                          direction along the edge from (0,0) to (0,1)
*                          of the previous patch.
*              iordcrt1  - Number of vertices of input cross tangent curve.
*                          iordcrt1 = 4.
*              evec1     - Tangent vector at the corner (0,0) in the 1. 
*                          parameter direction of the previous patch.
*              evec2     - Tangent vector at the corner (0,0) in the 2. 
*                          parameter direction of the previous patch, and
*                          in the 1. parameter direction of the current patch.
*              evec3     - Tangent vector at the corner (0,0) in the 2. 
*                          parameter direction of the current patch.
*              idim      - Dimension of geometry space.
*
*
* OUTPUT     : ecrt2     - Vertices of the produced cross tangent curve. 
*                          This curve is the derivative curve in the 2.
*                          parameter direction along the edge from (0,0) to
*                          (1,0) of the current patch.
*              jordcrt2  - Number of vertices of produced cross tangent curve.
*                          *jordcrt2 = 7.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Let pc1 be the input cross tangent curve, pc2 be the
*              derivative curve along the position curve and alpha the 
*              blending function. Then the output cross tangent curve is 
*              given by :
*              pc3(s) = ((my+1)*alpha(s)-1)*pc1(s) + lambda*alpha(s)*pc2(s)
*              my and lambda are constants depending on the lengths of 
*              input vectors, evec1, evec2 and evec3, and the angles 
*              between them. See the paper of Hahn.
*
* CALLS     : s6length - Length of vector.
*
*********************************************************************
*/
{
  int kstat;                       /* Status variable.              */
  int ki;                          /* Counters.                     */
  double tl1,tl2,tl3;              /* Lengths of the input vectors. */
  double tang1,tang2;              /* Angles between input vectors. */
  double tsin1,tsin2,tsin3;        /* Sinus of angles between input
                                      vectors.                      */
  double tmy,tmy1;                 /* Coefficients dependant on the
                                      input vectors.                */
  double tlambda;                  /* Coefficient dependant on the
                                      input vectors.                */
  double scoef[12];                /* Vertices of the derivative
                                      curve along the position curve. */
  double sc1[21],sc2[21];          /* Coefficients of products of curves. */
  double sblend2[4];               /* Coefficients of curve equal to
                                      a factor times the blending 
                                      curve minus 1.                 */
  
  if (iordblend != 4 || iordcrt1 != 4) goto err002;   

  *jordcrt2 = iordblend + iordcrt1 - 1;
  
  /* Compute constants.  */

  tl1 = s6length(evec1,idim,&kstat);
  tl2 = s6length(evec2,idim,&kstat);
  tl3 = s6length(evec3,idim,&kstat);
  tang1 = sh1461_s9ang(evec1,evec2,idim);
  tang2 = sh1461_s9ang(evec2,evec3,idim);
  tsin1 = sin(tang1);
  tsin2 = sin(tang2);
  tsin3 = sin(tang1+tang2);
  tmy = - (tl3*tsin2)/(tl1*tsin1);
  tlambda = (tl3*tsin3)/(tl2*tsin1);
  tmy1 = tmy + (double)1.0;
  
  /* Compute coefficients of the curve (my+1)alpha(s)-1.  */

  for (ki=0; ki<4; ki++) sblend2[ki] = tmy1*eblend[ki] - (double)1.0;
  
  /* Compute coefficients of derivative curve along position curve. */

  for (ki=0; ki<4*idim; ki++)
    scoef[ki] = (double)4.0*(epos[idim+ki] - epos[ki]);

   /* Compute first part of the second cross tangent curve. */

  sh1461_s9mult(eblend,scoef,4,idim,sc1,&kstat);

   /* Compute second part of the second cross tangent curve. */

  sh1461_s9mult(sblend2,ecrt1,4,idim,sc2,&kstat);
  
  /* Compute cross tangent curve.  */

  for (ki=0; ki<7*idim; ki++)
    ecrt2[ki] = (tlambda*sc1[ki] + sc2[ki]);
  
  *jstat = 0;
  goto out;
  
  /* Error in input.  */

  err002 :
    *jstat = -2;
  goto out;
  
  out :
    return;
}


#if defined(SISLNEEDPROTOTYPES)
static void
  sh1461_s9mult(double eblend[],double ecoef[],int iord,
		int idim,double ecoefnew[],int *jstat)
#else	       
static void sh1461_s9mult(eblend,ecoef,iord,idim,ecoefnew,jstat)
     int iord,idim,*jstat;
     double eblend[],ecoef[],ecoefnew[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the product of two Bezier curves of order 4,
*              one of which is a blending function of dimension 1.
*
*
*
* INPUT      : eblend    - Vertices of blending function. The dimension
*                          of the array is iord, i.e. 4.
*              ecoef     - Vertices of the other Bezier curve. The 
*                          dimension of the array is iord*idim, i.e. 4*idim.
*              iord      - Order of the Bezier curves. iord = 4.
*              idim      - Dimension of geometry space.
*
*
* OUTPUT     : ecoefnew  - Vertices of the product curve. The array is
*                          allocated outside this routine, and has
*                          dimension (2*iord-1)*idim, i.e. 7*idim
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* USE        : Used in Hahn's method when computing the cross tangent
*              curve which is a mapping of the derivative curves of the
*              previous patch. Called from s9chcoor.
*
*********************************************************************
*/
{
  int ki;   /* Index.   */
  
  /* Test if the order of the curves is equal to 4.  */

  if (iord != 4) goto err001;
  
  for (ki=0; ki<idim; ki++)
    {
      /* Compute the vertices of the product curve. */

      ecoefnew[ki] = eblend[0]*ecoef[ki];
      ecoefnew[idim+ki] = (eblend[0]*ecoef[idim+ki] 
			   + eblend[1]*ecoef[ki])/(double)2.0;
      ecoefnew[2*idim+ki] = (eblend[0]*ecoef[2*idim+ki]
			     + (double)3.0*eblend[1]*ecoef[idim+ki]
			     + eblend[2]*ecoef[ki])/(double)5.0;
      ecoefnew[3*idim+ki] = (eblend[0]*ecoef[3*idim+ki] + eblend[3]*ecoef[ki]
			     + (double)9.0*(eblend[1]*ecoef[2*idim+ki] + 
					    eblend[2]*ecoef[idim+ki]))/(double)20.0;
      ecoefnew[4*idim+ki] = (eblend[1]*ecoef[3*idim+ki]
			     + (double)3.0*eblend[2]*ecoef[2*idim+ki]
			     + eblend[3]*ecoef[idim+ki])/(double)5.0;      
      ecoefnew[5*idim+ki] = (eblend[2]*ecoef[3*idim+ki] 
			     + eblend[3]*ecoef[2*idim+ki])/(double)2.0;
      ecoefnew[6*idim+ki] = eblend[3]*ecoef[3*idim+ki];
    }

  *jstat = 0;
  goto out;
  
  /* Error in input. The order of the curves is not equal to 4.  */

  err001 :
    *jstat = -1;
  goto out;
  
  out :
    return;
}

#if defined(SISLNEEDPROTOTYPES)
static double
  sh1461_s9ang(double evec1[],double evec2[],int idim)
#else	       
static double sh1461_s9ang(evec1,evec2,idim)
     double evec1[];
     double evec2[];
     int    idim;
#endif     
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : sh1461_s9ang   - Angle in radians between vectors
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   - Scalar product between two vectors.
*              s6length - Length of vector.
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*              Vibeke Skytt SI, 90-04.
*
*********************************************************************
*/                                     
{
  double tscpr,tang,tlength1,tlength2,tcos;
  int    kstat1,kstat2;
  
  tscpr = s6scpr(evec1,evec2,idim);
  
  tlength1 = s6length(evec1,idim,&kstat1);
  tlength2 = s6length(evec2,idim,&kstat2);
  
  if (!kstat1 || !kstat2)
    tang = DZERO;
  else
    {
      tcos = tscpr/(tlength1*tlength2);
      tcos = MIN((double)1.0,tcos);
      tcos = MAX(-(double)1.0,tcos);
      tang = acos(tcos);
    }
      
  return(tang);
}

#if defined(SISLNEEDPROTOTYPES)
static  void
  sh1461_s9comder(int ider1,int ider2,double ederprev[],int idim,
		  double aro00,double aro01,double aro10,
		  double aro11,double eder[],int *jstat)
#else	    
static void sh1461_s9comder(ider1,ider2,ederprev,idim,aro00,aro01,aro10,
			    aro11,eder,jstat)
     int ider1,ider2,idim,*jstat;
     double ederprev[],aro00,aro01,aro10,aro11,eder[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute given 2. derivative of a blending surface at the 
*              midpoint of the vertex region when the 2. derivatives of the 
*              first blending patch are known.
*
*
*
* INPUT      : ider1    - Order of differentiation in first parameter direction.
*              ider2    - Order of differentiation in second parameter direction.
*              ederprev - 2. derivatives of 1. patch stored in the following 
*                         order, (2,0)-, (1,1)- and (0,2)-derivative.
*              idim     - Dimension of geometry space.
*              aro00    - First factor of coordinate transformation.
*              aro01    - Second factor of coordinate transformation.
*              aro10    - Third factor of coordinate transformation.
*              aro11    - Fourth factor of coordinate transformation.
*
*
* OUTPUT     : eder    - Actual 2. derivative.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
  int kj,kk,kh;            /* Counters.  */
  double t00,t01,t10,t11;  /* Coefficients used to compute
			      derivatives at midpoint of region. */
  double tfac1,tfac2;      /* Factor used to compute derivatives 
			      at midpoint of region.  */
  
  /* Test input.  */

  if (ider1 + ider2 != 2) goto err001;
  
  /* Compute requested 2. derivative.  */

  t00 = t10 = (double)1.0;
  for (kj=0; kj<ider1; kj++) t00 *= aro00;

  for (kj=0; kj<=ider1; kj++)
    {
      tfac1 = (ider1 > 0) ? (double)((kj % ider1) + 1) : (double)1.0;

      t01 = t11 = (double)1.0;
      for (kk=0; kk<ider2; kk++) t01*=aro01;

      for (kk=0; kk<=ider2; kk++)
	{
	  tfac2 = (ider2 > 0) ? (double)((kk % ider2) + 1) : (double)1.0;

	  for (kh=0; kh<idim; kh++)
	    eder[kh] += tfac1*tfac2*t00*t01*t10*t11*ederprev[(2-kj-kk)*idim+kh];
	      
	  if (aro01 != DZERO) t01 /= aro01;
	  else if (kk == ider2-1) t01 = (double)1.0;
	  t11 *= aro11;
	}
      if (aro00 != DZERO) t00 /= aro00;
      else if (kj == ider1-1) t00 = (double)1.0;
      t10 *= aro10;	  
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Wrong order of differentiation.  */

  err001 :
    *jstat = -1;
  goto out;
  
  out :
    return;
}
