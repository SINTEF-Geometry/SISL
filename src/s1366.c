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
 * $Id: s1366.c,v 1.7 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1366

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1366(SISLSurf *ps, double aoffset, double aepsge, double amax, int idim,
	   double *eknot13, int in13, int ik13,
	   double *eknot24, int in24, int ik24,
	   SISLSurf **rs, int *jstat)
#else
void s1366(ps,aoffset,aepsge,amax,idim,eknot13,in13,ik13,
	   eknot24,in24,ik24,rs,jstat)
     SISLSurf   *ps;
     double aoffset;
     double aepsge;
     double amax;
     int    idim;
     double *eknot13;
     int    in13;
     int    ik13;
     double *eknot24;
     int    in24;
     int    ik24;
     SISLSurf   **rs;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To create a B-spline approximating the offset surface of
*              a B-spline surface.
*
*
*
* INPUT      : ps     - The input B-spline surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              aepsge - Maximal deviation allowed between true offset surface
*                       and the approximated offset surface.
*              amax   - Maximal stepping length. Is negleceted if amax<=aepsge
*              idim   - The dimension of the space (2 or 3).
*              eknot13- Pointer to common knot-vector along first parameter
*                       direction.
*              in13   - Number of vertices of knot-vector along first
*                       parameter direction.
*              ik13   - Order of knot-vector along first parameter direction.
*              eknot24- Pointer to common knot-vector along second parameter
*                       direction.
*              in24   - Number of vertices of knot-vector along second
*                       parameter direction.
*              ik24   - Order of knot-vector along second parameter direction.
*
* OUTPUT     :
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer the approximated offset surface
*
* METHOD     : Points, derivatives, cross derivatives and normals are calcu-
*              lated on selected points on the surface. The selected points
*              corresponds to parameter-pairs of all specific values of the
*              knot-vector. The offset values of these points are calculated.
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*-
* CALLS      : s1367     - Evaluate the surface at a given parameter pair
*                          value.
*              s1347     - Interpolate a bicubic hermite spline-approximation
*                          surface.
*
* WRITTEN BY : Per Evensen,  SI, 89-4.
* REWISED BY : Per Evensen,  SI, 90-9 Corrected pointers.
* REWISED BY : Paal Fugelli, SINTEF, 1994-07 Added free'ing of allocated
*              memory (at end) to fix memory leakage problem.
*              1994-09 Added error checks after calls to increasearray() macros
*                      and changed usage style of 'kmaxik'.
*********************************************************************
*/
{
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int ki;            /* Loop controller. */
  int kk1;           /* Loop controller. */
  int kk2;           /* Loop controller. */

  int kder = 1;      /* Number of derivatives necessary to achieve.      */
  int klfs=0;       /* Pointer to knot interval.                        */
  int klft=0;       /* Pointer to knot interval.                        */
  int knbpnt=0;      /* Number of points stored. */
  int km1=1;         /* Counter of points in 1. parameter direction. */
  int km2=1;         /* Counter of points in 2. parameter direction. */
  int kain=0;        /* Array index. */
  int kmaxis=512;    /* Number of vertices space is allocated for       */
  int kmaxik1=512;   /* Number of vertices along first parameter direction which
			space is allocated for. */
  int kmaxik2=512;   /* Number of vertices along second parameter direction which
			space is allocated for. */
  int kpar=3;        /* Flag determining the parametrization of the data-points.
			= 1: Mean accumulated cord length parametrization.
			= 2: Uniform parametrization.
			= 3: Parametrization given by spar1 and spar2. */
  int lend[4];       /* Array containing the no. of derivatives to be kept fixed
			along each of the edges of teh surface. */
  int kopt=3;        /* Flag indicating the order in which the data-reduction
			is to be performed. */
  int ktmax=10;      /* Max. number of itarations. */
  double *spar1;     /* Array (length km1) containing a parametrization in the
			1. parameter direction. */
  double *spar2;     /* Array (length km2) containing a parametrization in the
			2. parameter direction. */
  double seps[3];    /* Array containing the maximum deviation which is
			acceptable in each of the idim components of the
			surface (except possibly along the edges). */
  double sedgeps[12];/* Array containing the max. deviation which is acceptable
			along the edges of teh surface. */
  double tepsco;     /* Computer resolution. */
  double smxerr[3];  /* Array containing an upper bound for the error comitted
			in each component during the data reduction. */
  double *spnt=SISL_NULL; /* Pointer to storage of point info. */
  double *stng1=SISL_NULL;/* Pointer to derivatives in 1. parameter direction. */
  double *stng2=SISL_NULL;/* Pointer to derivatives in 2. parameter direction. */
  double *scrss=SISL_NULL; /* Pointer to cross derivatives. */
  double spar[2];    /* Parameter value at where to evaluate surface.    */
  double sder[27];   /* Pointer to array containing the derivatives.     */

  /* Initialization of variables */
  tepsco = (double)0.000001;

  for (ki=0; ki<3; ki++) seps[ki] = aepsge;

  for (ki=0; ki<4; ki++) lend[ki] = 1;

  for (ki=0; ki<12; ki++) sedgeps[ki] = aepsge;

  /* Allocate space for storage of points, derivatives and cross-derivatives. */

  spnt = newarray(idim*kmaxis,DOUBLE);
  if (spnt == SISL_NULL) goto err101;

  stng1 = newarray(idim*kmaxis,DOUBLE);
  if (stng1 == SISL_NULL) goto err101;

  stng2 = newarray(idim*kmaxis,DOUBLE);
  if (stng2 == SISL_NULL) goto err101;

  scrss = newarray(idim*kmaxis,DOUBLE);
  if (scrss == SISL_NULL) goto err101;

  /* Allocate space for storage of parametrization in 1. direction. */

  spar1 = newarray(kmaxik1,DOUBLE);
  if (spar1 == SISL_NULL) goto err101;

  /* Allocate space for storage of parametrization in 2. direction. */

  spar2 = newarray(kmaxik2,DOUBLE);
  if (spar2 == SISL_NULL) goto err101;

  /* Initiate parameter value in first and second parameter direction. */
  spar[0] = eknot13[ik13-1];
  spar[1] = eknot24[ik24-1];

  /* Evaluate the surface (ps) at the parameter value spar.  */
  s1367(ps,aoffset,aepsge,idim,spar,kder,&klfs,&klft,sder,&kstat);
  if (kstat<0) goto error;

  /* Store information about points, derivatives and cross-derivative. */
  for (ki=0; ki<idim; ki++)
  {
    spnt[kain] = sder[ki];
    stng1[kain] = sder[ki+idim];
    stng2[kain] = sder[ki+2*idim];
    scrss[kain] = sder[ki+3*idim];
    kain += 1;
  }

  /* Store parametrization along both parameter directions. */
  spar1[km1-1] = spar[0];
  spar2[km2-1] = spar[1];

  knbpnt += 1;

  for (kk1=ik13; kk1<=in13; kk1++)
  {

    /* check if current knot-value in first parameter direction is
       equal the previous or not. */
    if (DNEQUAL(spar[0],eknot13[kk1]))
    {

      /* Evaluate the surface (ps) at the parameter value spar.  */
      spar[0] = eknot13[kk1];
      s1367(ps,aoffset,aepsge,idim,spar,kder,&klfs,&klft,sder,&kstat);
      if (kstat<0) goto error;

      /* Store information about points, derivatives and cross-derivative. */
      knbpnt += 1;
      km1 += 1;

      /* Possibly increase size of arrays. */
      if (knbpnt>=kmaxis)
      {
	kmaxis += 512;
	spnt  = increasearray(spnt,idim*kmaxis,DOUBLE);
	if ( spnt == SISL_NULL )  goto err101;
	stng1 = increasearray(stng1,idim*kmaxis,DOUBLE);
	if ( stng1 == SISL_NULL )  goto err101;
	stng2 = increasearray(stng2,idim*kmaxis,DOUBLE);
	if ( stng2 == SISL_NULL )  goto err101;
	scrss = increasearray(scrss,idim*kmaxis,DOUBLE);
	if ( scrss == SISL_NULL )  goto err101;
      }

      /* Possibly increase size of arrays. */
      if (km1>=kmaxik1)
      {
	kmaxik1 = km1 + 512;	/* kmaxik += 10; (PFU 19/09-94) */
	spar1 = increasearray(spar1,kmaxik1,DOUBLE);
	if ( spar1 == SISL_NULL )  goto err101;
      }

      for (ki=0; ki<idim; ki++)
      {
	spnt[kain] = sder[ki];
	stng1[kain] = sder[ki+idim];
	stng2[kain] = sder[ki+2*idim];
	scrss[kain] = sder[ki+3*idim];
	kain += 1;
      }

      /* Store parametrization along 1. parameter direction. */
      spar1[km1-1] = spar[0];
    }
  }


  for (kk2=ik24; kk2<=in24; kk2++)
  {

    /* check if current knot-value in second parameter direction is
       equal the previous or not. */
    if (DNEQUAL(spar[1],eknot24[kk2]))
    {

      /* Initiate parameter value in first parameter direction. */
      spar[0] = eknot13[ik13-1];
      klfs = 0;

      /* Evaluate the surface (ps) at the parameter value spar.  */
      spar[1] = eknot24[kk2];
      s1367(ps,aoffset,aepsge,idim,spar,kder,&klfs,&klft,sder,&kstat);
      if (kstat<0) goto error;

      /* Store information about points, derivatives and cross-derivative. */
      knbpnt += 1;
      km2 += 1;

      /* Possibly increase size of arrays. */
      if (knbpnt>=kmaxis)
      {
	kmaxis += 512;
	spnt  = increasearray(spnt,idim*kmaxis,DOUBLE);
	if ( spnt == SISL_NULL )  goto err101;
	stng1 = increasearray(stng1,idim*kmaxis,DOUBLE);
	if ( stng1 == SISL_NULL )  goto err101;
	stng2 = increasearray(stng2,idim*kmaxis,DOUBLE);
	if ( stng2 == SISL_NULL )  goto err101;
	scrss = increasearray(scrss,idim*kmaxis,DOUBLE);
	if ( scrss == SISL_NULL )  goto err101;
      }

      /* Possibly increase size of arrays. */
      if (km2>=kmaxik2)
      {
	kmaxik2 = km2 + 512;  /* kmaxik += 10; (PFU 19/09-94) */
	spar2 = increasearray(spar2,kmaxik2,DOUBLE);
	if ( spar1 == SISL_NULL )  goto err101;
      }

      for (ki=0; ki<idim; ki++)
      {
	spnt[kain] = sder[ki];
	stng1[kain] = sder[ki+idim];
	stng2[kain] = sder[ki+2*idim];
	scrss[kain] = sder[ki+3*idim];
	kain += 1;
      }

      /* Store parametrization along 1. parameter direction. */
      spar2[km2-1] = spar[1];

      for (kk1=ik13; kk1<=in13; kk1++)
      {

	/* check if current knot-value in first parameter direction is
	   equal the previous or not. */
	if (DNEQUAL(spar[0],eknot13[kk1]))
	{
	  spar[0] = eknot13[kk1];

	  /* Evaluate the surface (ps) at the parameter value spar.  */
	  s1367(ps,aoffset,aepsge,idim,spar,kder,
		&klfs,&klft,sder,&kstat);
	  if (kstat<0) goto error;

	  /* Store information about points, derivatives and
	     cross-derivative. */
	  knbpnt += 1;

	  /* Possibly increase size of arrays. */
	  if (knbpnt>=kmaxis)
	  {
	    kmaxis += 512;
	    spnt  = increasearray(spnt,idim*kmaxis,DOUBLE);
	    if ( spnt == SISL_NULL )  goto err101;
	    stng1 = increasearray(stng1,idim*kmaxis,DOUBLE);
	    if ( stng1 == SISL_NULL )  goto err101;
	    stng2 = increasearray(stng2,idim*kmaxis,DOUBLE);
	    if ( stng2 == SISL_NULL )  goto err101;
	    scrss = increasearray(scrss,idim*kmaxis,DOUBLE);
	    if ( scrss == SISL_NULL )  goto err101;
	  }

	  for (ki=0; ki<idim; ki++)
	  {
	    spnt[kain] = sder[ki];
	    stng1[kain] = sder[ki+idim];
	    stng2[kain] = sder[ki+2*idim];
	    scrss[kain] = sder[ki+3*idim];
	    kain += 1;
	  }
	}
      }
    }
  }

  /* Compute a bicubic hermite spline-approximation to the position and
     derivative data given by spnt, stng1, stng2 and scrss. */

  s1347(spnt,stng1,stng2,scrss,km1,km2,idim,kpar,spar1,spar2,seps,lend,
	sedgeps,tepsco,kopt,ktmax,rs,smxerr,&kstat);
  if (kstat<0) goto error;

  /* Surface approximated. */

  *jstat = 0;
  goto out;

  /* Error in memory allocation */

err101:
  *jstat = -101;
  s6err("s1366",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

error:
  *jstat = kstat;
  s6err("s1366",*jstat,kpos);
  goto out;

out:
  if (spnt != SISL_NULL) freearray(spnt);
  if (stng1 != SISL_NULL) freearray(stng1);
  if (stng2 != SISL_NULL) freearray(stng2);
  if (scrss != SISL_NULL) freearray(scrss);
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);

  return;
}
