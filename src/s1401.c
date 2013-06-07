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
 * $Id: s1401.c,v 1.5 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1401

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void s1401_s9blend(double [],int,int,int *);
static void s1401_s9basis1(int,int,double [],int,int,int,double [],
				 int,int,double [],double [],double **,int *);
static void s1401_s9basis2(int,int,double [],int,int,int,double [],
				 int,int,double [],double [],double **,int *);
#else
  static void   s1401_s9blend();
  static void   s1401_s9basis1();
  static void   s1401_s9basis2();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
   s1401(SISLCurve *vcurve[],double etwist[],
	      SISLSurf **rsurf,int *jstat)
#else
void s1401(vcurve,etwist,rsurf,jstat)
     double etwist[];
     SISLCurve *vcurve[];
     int *jstat;
     SISLSurf **rsurf;
#endif
/*
*********************************************************************
*
* PURPOSE    : Compute a Gordon patch given position and cross tangent
*              conditions at the boundary of a squared region and the
*              twist vector in the corners/
*
*
*
* INPUT      : vcurve - Position and cross-tangent curves around the square
*                       region. For each edge of the region position and cross-
*                       tangent curves are given. The dimension of the array is 8.
*                       The orientation is as follows :
*
*                                     3 ---->
*                               -------------------
*                               |                 |
*                           /|\ |                 | /|\
*                            |  |                 |  |
*                            |  |                 |  |
*                            4  |                 |  2
*                               |                 |
*                               -------------------
*                                     1 ---->
*
*              etwist - Twist-vectors of the corners of the vertex region. The
*                       first element of the array is the twist in the corner
*                       before the first edge, etc. The dimension of the array
*                       is 4*kdim.
*
*
* OUTPUT     : rsurf  - Gordons patch represented as a B-spline patch.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
* REFERENCES : I. D. Faux and M. J. Pratt :
*                        Computational Geometry for Design and Manufacture
*              J. Hahn : Filling Polygonal Holes with Rectangular Patches
*
* USE        : 3D non-rational geometry only.
*
*-
* CALLS      : s1221 - Evaluate curve at a given parameter value.
*              s1931unit
*              s1750
*              freeCurve
*              newSurf
*
* WRITTEN BY : Vibeke Skytt, SI, April 90.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Updated
*              comments from Coons to Gordons patches.
*
*********************************************************************
*/
{
  int kstat = 0;
  int kdim;
  int kder;
  int kcurve;
  int ki,kj,kk,kh,kl;
  int kleft = 0;
  int kn1;
  int kk1;
  int kn2;
  int kk2;
  int korder = 4;

  double tpar;
  double salpha[16];
  double sknot0[8];
  double *sder = SISL_NULL;
  double *sknot1 = SISL_NULL;
  double *scoef1 = SISL_NULL;
  double *sknot2 = SISL_NULL;
  double *scoef2 = SISL_NULL;
  double *ssurf1 = SISL_NULL;
  double *ssurf2 = SISL_NULL;
  double *ssurf3 = SISL_NULL;
  double *ssurf4 = SISL_NULL;
  double *smat = SISL_NULL;
  double *ssurf1new = SISL_NULL;
  double *ssurf2new = SISL_NULL;
  double *ssurf3new = SISL_NULL;

  SISLCurve *qc1[4];
  SISLCurve *qc2[4];

  for (ki=0; ki<4; ki++)
  {
     qc1[ki] = SISL_NULL;
     qc2[ki] = SISL_NULL;
  }

  /* Test dimension.  */

  kdim = vcurve[0] -> idim;
  for (ki=1; ki<8; ki++)
    if (vcurve[ki]->idim != kdim) goto err102;

  /* Allocate scratch for local array.  */

  if ((sder = newarray(3*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((smat = newarray(16*kdim,DOUBLE)) == SISL_NULL) goto err101;

  /* TEST !!! */

  /*
  kder = 1;
  for (ki=0; ki<4; ki++)
    {
      kj = (ki + 1) % 4;

      tpar = (ki >= 2) ? *(vcurve[2*ki]->et + vcurve[2*ki]->ik - 1) :
	*(vcurve[2*ki]->et + vcurve[2*ki]->in);

      s1221(vcurve[2*ki+1],kder,tpar,&kleft,sder,&kstat);
      if (kstat < 0) goto error;

      memcopy(sder,sder+kdim,kdim,DOUBLE);

      tpar = (kj < 2) ? *(vcurve[2*kj]->et + vcurve[2*kj]->ik - 1) :
	*(vcurve[2*kj]->et + vcurve[2*kj]->in);

      s1221(vcurve[2*kj+1],kder,tpar,&kleft,sder+kdim,&kstat);
      if (kstat < 0) goto error;
    }
  */


  /* Set up blending functions.  */

  salpha[0] = (double)1.0;
  salpha[1] = salpha[2] = salpha[3] = (double)0.0;
  s1401_s9blend(salpha,4,1,&kstat);

  salpha[4] = salpha[5] = salpha[6] = (double)0.0;
  salpha[7] = (double)1.0;
  s1401_s9blend(salpha+4,4,1,&kstat);

  salpha[9] = (double)1.0;
  salpha[8] = salpha[10] = salpha[11] = (double)0.0;
  s1401_s9blend(salpha+8,4,1,&kstat);

  salpha[14] = (double)1.0;
  salpha[12] = salpha[13] = salpha[15] = (double)0.0;
  s1401_s9blend(salpha+12,4,1,&kstat);

  /* Initiate knot vector of blending functions.  */

  for (ki=0; ki<4; ki++)
    {
      sknot0[ki] = (double)0.0;
      sknot0[4+ki] = (double)1.0;
    }

  /* Evaluate the boundary curves in the corners of the region.  */

  for (ki=0; ki<4; ki++)
    {
      tpar = (ki < 2) ? *(vcurve[2*ki]->et + vcurve[2*ki]->ik - 1) :
	*(vcurve[2*ki]->et + vcurve[2*ki]->in);

      kder = 1;

      s1221(vcurve[2*ki],kder,tpar,&kleft,sder,&kstat);
      if (kstat < 0) goto error;

      kder = 0;
      s1221(vcurve[2*ki+1],kder,tpar,&kleft,sder+2*kdim,&kstat);
      if (kstat < 0) goto error;

      /* Copy into matrix of correction term.  */

      kj = (ki > 1);
      kk = MIN(ki % 3,1);

      memcopy(smat+(kk*4+kj)*kdim,sder,kdim,DOUBLE);
      memcopy(smat+((kk+2)*4+kj)*kdim,sder+(((kk+kj)%2)+1)*kdim,kdim,DOUBLE);
      memcopy(smat+(kk*4+kj+2)*kdim,sder+(2-((kk+kj)%2))*kdim,kdim,DOUBLE);

    }

  for (kh=0; kh<kdim; kh++)
    {
      smat[3*kdim+kh] *= -(double)1.0;
      smat[13*kdim+kh] *= -(double)1.0;
    }

  /* Copy twist vector into matrix of correction term.  */

  memcopy(smat+10*kdim,etwist,kdim,DOUBLE);
  memcopy(smat+11*kdim,etwist+3*kdim,kdim,DOUBLE);
  memcopy(smat+14*kdim,etwist+kdim,kdim,DOUBLE);
  memcopy(smat+15*kdim,etwist+2*kdim,kdim,DOUBLE);

  /* Put all curves along the standard edges 1 and 3 on the same knotvector.
     First fetch the actual curves.  */

  kcurve = 4;
  qc1[0] = vcurve[0];
  qc1[1] = vcurve[4];
  qc1[2] = vcurve[1];
  qc1[3] = vcurve[5];

  for (ki=0; ki<4; ki++)
    {
      if (qc1[ki]->ik < 4)
	{
	  /* Raise the order of the curve. */

	  s1750(qc1[ki],korder,qc1+ki,&kstat);
	  if (kstat < 0) goto error;
	}
    }

  /* Express the curves in a common basis.  */

  s1931unit(kcurve,qc1,&sknot1,&scoef1,&kn1,&kk1,&kstat);
  if (kstat < 0) goto error;

  /* Put all curves along the standard edges 2 and 4 on the same knotvector.
     First fetch the actual curves.  */

  qc2[0] = vcurve[6];
  qc2[1] = vcurve[2];
  qc2[2] = vcurve[7];
  qc2[3] = vcurve[3];

  for (ki=0; ki<4; ki++)
    {
      if (qc2[ki]->ik < 4)
	{
	  /* Raise the order of the curve. */

	  s1750(qc2[ki],korder,qc2+ki,&kstat);
	  if (kstat < 0) goto error;
	}
    }

  /* Express the curves in a common basis.  */

  s1931unit(kcurve,qc2,&sknot2,&scoef2,&kn2,&kk2,&kstat);
  if (kstat < 0) goto error;

  /* Allocate scratch for coefficients of blending surfaces.  */

  if ((ssurf1 = new0array(4*kn1*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((ssurf2 = new0array(kn2*4*kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((ssurf3 = new0array(16*kdim,DOUBLE)) == SISL_NULL) goto err101;

  /* Compute coefficients of 1. blending surface.  */

  for (ki=0; ki<4; ki++)
    for (kj=0; kj<4; kj++)
      for (kk=0; kk<kn1; kk++)
	for (kh=0; kh<kdim; kh++)
	  ssurf1[(kk*4+kj)*kdim+kh] +=
	    salpha[ki*4+kj]*scoef1[(ki*kn1+kk)*kdim+kh];

  /* Compute coefficients of 2. blending surface.  */

  for (ki=0; ki<4; ki++)
    for (kk=0; kk<kn2; kk++)
      for (kj=0; kj<4; kj++)
	for (kh=0; kh<kdim; kh++)
	  ssurf2[(kj*kn2+kk)*kdim+kh] +=
	    scoef2[(ki*kn2+kk)*kdim+kh]*salpha[ki*4+kj];

  /* Compute the correction surface.  */

  for (ki=0; ki<4; ki++)
    for (kj=0; kj<4; kj++)
      for (kk=0; kk<4; kk++)
	for (kl=0; kl<4; kl++)
	  for (kh=0; kh<kdim; kh++)
	    ssurf3[(kk*4+kl)*kdim+kh] +=
	      salpha[ki*4+kk]*salpha[kj*4+kl]*smat[(ki*4+kj)*kdim+kh];

  /* Express all 3 surfaces in the same basis.  */

  s1401_s9basis1(kk2,kn2,sknot2,kdim,4,4,sknot0,kk1,kn1,sknot1,ssurf1,
	  &ssurf1new,&kstat);
  if (kstat < 0) goto error;

  s1401_s9basis2(kk1,kn1,sknot1,kdim,kk2,kn2,sknot2,4,4,sknot0,ssurf2,
	  &ssurf2new,&kstat);
  if (kstat < 0) goto error;

  s1401_s9basis1(kk2,kn2,sknot2,kdim,4,4,sknot0,4,4,sknot0,ssurf3,
	  &ssurf4,&kstat);
  if (kstat < 0) goto error;

  s1401_s9basis2(kk1,kn1,sknot1,kdim,kk2,kn2,sknot2,4,4,sknot0,ssurf4,
	  &ssurf3new,&kstat);
  if (kstat < 0) goto error;

  /* Compute the coefficients of Gordons patch.  */

  for (ki=0; ki<kn1*kn2*kdim; ki++)
    ssurf1new[ki] += ssurf2new[ki] - ssurf3new[ki];

  /* Store Gordons patch.  */

  if ((*rsurf = newSurf(kn2,kn1,kk2,kk1,sknot2,sknot1,ssurf1new,1,kdim,1)) == SISL_NULL)
    goto err101;

  /* Gordons patch defined.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;

  /* Error in input. Dimension of curves not equal.  */

  err102 :
    *jstat = -102;
  goto out;

  /* Error in lower level function.  */

  error :
    *jstat = kstat;
  goto out;

  out :

    /* Free space occupied by local arrays.  */

    if (sder != SISL_NULL) freearray(sder);
  if (smat != SISL_NULL) freearray(smat);
  if (qc1[0] && qc1[0] != vcurve[0]) freeCurve(qc1[0]);
  if (qc1[1] && qc1[1] != vcurve[4]) freeCurve(qc1[1]);
  if (qc1[2] && qc1[2] != vcurve[1]) freeCurve(qc1[2]);
  if (qc1[3] && qc1[3] != vcurve[5]) freeCurve(qc1[3]);
  if (qc2[0] && qc2[0] != vcurve[6]) freeCurve(qc2[0]);
  if (qc2[1] && qc2[1] != vcurve[2]) freeCurve(qc2[1]);
  if (qc2[2] && qc2[2] != vcurve[7]) freeCurve(qc2[2]);
  if (qc2[3] && qc2[3] != vcurve[3]) freeCurve(qc2[3]);
  if (sknot1 != SISL_NULL) freearray(sknot1);
  if (sknot2 != SISL_NULL) freearray(sknot2);
  if (scoef1 != SISL_NULL) freearray(scoef1);
  if (scoef2 != SISL_NULL) freearray(scoef2);
  if (ssurf1 != SISL_NULL) freearray(ssurf1);
  if (ssurf2 != SISL_NULL) freearray(ssurf2);
  if (ssurf3 != SISL_NULL) freearray(ssurf3);
  if (ssurf4 != SISL_NULL) freearray(ssurf4);
  if (ssurf1new != SISL_NULL) freearray(ssurf1new);
  if (ssurf2new != SISL_NULL) freearray(ssurf2new);
  if (ssurf3new != SISL_NULL) freearray(ssurf3new);

    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
   s1401_s9blend(double econd[],int icond,int idim,int *jstat)
#else
static void s1401_s9blend(econd,icond,idim,jstat)
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
*              idim       - Dimension of geometry space.
*
*
* INPUT/OUTPUT : econd    - Interpolation conditions as input, Bezier coefficients
*                           as output.
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
  int ki;

  if (icond != 4) goto err001;

  for (ki=0; ki<idim; ki++)
    {
      econd[idim+ki] = ONE_THIRD*econd[idim+ki] + econd[ki];
      econd[2*idim+ki] = ONE_THIRD*econd[2*idim+ki] + econd[3*idim+ki];
    }

  *jstat = 0;
  goto out;

  err001 :
    *jstat = -1;
  goto out;

  out :
    return;
}


#if defined(SISLNEEDPROTOTYPES)
static void
   s1401_s9basis1(int ik1new,int in1new,double et1new[],
		       int idim,int ik1,int in1,double et1[],
		       int ik2,int in2,double et2[],double ecoef[],
		       double **gcoefnew,int *jstat)
#else
static void s1401_s9basis1(ik1new,in1new,et1new,idim,ik1,in1,et1,
				 ik2,in2,et2,ecoef,gcoefnew,jstat)
     int ik1new,in1new,idim,ik1,in1,ik2,in2,*jstat;
     double et1new[],et1[],et2[],ecoef[],**gcoefnew;
#endif
/*
*********************************************************************
*
* PURPOSE    : Express a surface on an extended knot vector in the 1.
*              parameter direction.
*
*
*
* INPUT      : ik1new     - New order in the 1. parameter direction.
*              in1new     - New number of vertices in 1. parameter direction.
*              et1new     - New knot vector in the 1. parameter direction.
*              ik1        - Current order in 1. parameter direction.
*              in1        - Current number of vertices in 1. parameter direction.
*              et1        - Current knot vector in 1. parameter direction.
*              ik2        - Order in 2. parameter direction.
*              in2        - Number of vertices in 2. parameter direction.
*              et2        - Knot vector in 2. parameter direction.
*              ecoef      - Vertices of current surface.
*              idim       - Dimension of geometry space.
*
*
*
* OUTPUT     : gcoefnew   - Vertices of the new surface.
*              jstat      - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{

  int kstat = 0;
  int kdim = idim*in2;
  int knbcrv = 1;
  int kkind = 1;
  int kcopy = 1;
  double tstart,tstop;
  double *scoef = SISL_NULL;
  double *scoef2 = SISL_NULL;
  SISLCurve *qc = SISL_NULL;

  tstart = et1new[ik1new-1];
  tstop = et1new[in1new];

  /* Allocate scratch for coefficient arrays.  */

  if ((scoef = newarray(in1*in2*idim,DOUBLE)) == SISL_NULL) goto err101;
  if ((*gcoefnew = newarray(in1new*in2*idim,DOUBLE)) == SISL_NULL)
     goto err101;

  /* Change parameter directions of surface.  */

  s6chpar(ecoef,in1,in2,idim,scoef);

  /* Express the surface with interchanged parameter directions
     as a curve.  */

  if ((qc = newCurve(in1,ik1,et1,scoef,kkind,kdim,kcopy)) == SISL_NULL) goto err101;

  /* Express the surface interpreted as a curve on the extended
     knot vector.  */

  s1932(knbcrv,&qc,tstart,tstop,et1new,in1new,ik1new,&scoef2,&kstat);
  if (kstat < 0) goto error;

  /* Change parameter directions of output surface. */

  s6chpar(scoef2,in2,in1new,idim,*gcoefnew);

  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;

  /* Error in lower level function.  */

  error :
    *jstat = kstat;
  goto out;

  out :
     /* Free scratch.  */

    if (scoef != SISL_NULL) freearray(scoef);
    if (scoef2 != SISL_NULL) freearray(scoef2);
    if (qc != SISL_NULL) freeCurve(qc);

    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
   s1401_s9basis2(int ik2new,int in2new,double et2new[],
		       int idim,int ik1,int in1,double et1[],
		       int ik2,int in2,double et2[],double ecoef[],
		       double **gcoefnew,int *jstat)
#else
static void s1401_s9basis2(ik2new,in2new,et2new,idim,ik1,in1,et1,
				 ik2,in2,et2,ecoef,gcoefnew,jstat)
     int ik2new,in2new,idim,ik1,in1,ik2,in2,*jstat;
     double et2new[],et1[],et2[],ecoef[],**gcoefnew;
#endif
/*
*********************************************************************
*
* PURPOSE    : Express a surface on an extended knot vector in the 2.
*              parameter direction.
*
*
*
* INPUT      : ik2new     - New order in the 2. parameter direction.
*              in2new     - New number of vertices in 2. parameter direction.
*              et2new     - New knot vector in the 2. parameter direction.
*              ik1        - Order in 1. parameter direction.
*              in1        - Number of vertices in 1. parameter direction.
*              et1        - Knot vector in 1. parameter direction.
*              ik2        - Current order in 2. parameter direction.
*              in2        - Current number of vertices in 2. parameter direction.
*              et2        - Current knot vector in 2. parameter direction.
*              ecoef      - Vertices of current surface.
*              idim       - Dimension of geometry space.
*
*
*
* OUTPUT     : gcoefnew   - Vertices of the new surface.
*              jstat      - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* CALLS     : s1932, newCurve
*
*********************************************************************
*/
{
  int kstat = 0;
  int kdim = idim*in1;
  int knbcrv = 1;
  int kkind = 1;
  int kcopy = 1;
  double tstart,tstop;
  SISLCurve *qc = SISL_NULL;

  tstart = et2new[ik2new-1];
  tstop = et2new[in2new];

  /* Express the surface as a curve.  */

  if ((qc = newCurve(in2,ik2,et2,ecoef,kkind,kdim,kcopy)) == SISL_NULL) goto err101;

  /* Express the surface interpreted as a curve on the extended
     knot vector.  */

  s1932(knbcrv,&qc,tstart,tstop,et2new,in2new,ik2new,gcoefnew,&kstat);
  if (kstat < 0) goto error;

  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;

  /* Error in lower level function.  */

  error :
    *jstat = kstat;
  goto out;

  out :

     /* Free scratch.  */

     if (qc != SISL_NULL) freeCurve(qc);

    return;
}
