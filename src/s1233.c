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
 * $Id: s1233.c,v 1.4 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1233

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1233(SISLCurve *pc,double afak1,double afak2,SISLCurve **rc,int *jstat)
#else
void s1233(pc,afak1,afak2,rc,jstat)
     SISLCurve  *pc;
     double     afak1;
     double     afak2;
     SISLCurve  **rc;
     int    	*jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To extend a B-spline curve (i.e. NOT rationals) at the start
*              and/or the end of the curve by continuing the polynomial
*              behaviour of the curve.
*
* INPUT      : pc     - Pointer to the curve that will be extended.
*              afak1  - How much the curve is to be stretched at the
*                       start of the curve. The length of the stretched
*                       curve will be equal to (1+afak1) times the
*                       input curve. afak1 >= 0 and will be set to 0 if
*                       negative.
*              afak2  - How much the curve is to be stretched at the
*                       end of the curve. The length of the stretched
*                       curve will be equal to (1+afak2) times the
*                       input curve. afak2 >= 0 and will be set to 0 if
*                       negative.
*
* OUTPUT     : rc     - Pointer to the extended curve.
*              jstat  - status messages
*                       > 0      : warning
*                       = 1      : Stretching factors less than 0 - readjusted
*                                  factor(s) have been used.
*                       = 0      : ok
*                       < 0      : error
*
* METHOD     : 1. The new knot vector is made
*              2. The transformation matrix for the ik first vertices between
*                 the new and the old knot vector is made.
*              3. The transformation matrix is inverted and used to update the
*                 ik first vertices.
*              4. The transformation matrix for the ik last vertices between
*                 the new and the old knot vector is made.
*              5. The transformation matrix is inverted and used to update the
*                 ik last vertices.
*              6. The knot vector is updated.
*-
* CALLS      : make_cv_kreg,s1219,s1701,s6lufacp,s6lusolp,s6err.
*
* Written by : Paal Fugelli, SINTEF, Oslo, Norway, Sept-1992. Adapted from
*              a FORTRAN based s1233() written by Vibeke Skytt and s1333_cyclic()
*              written by Tor Dokken.
*
*********************************************************************
*/
{
  double *ext = SISL_NULL;                     /* Extended version of knot vector */
  double *smatrix = SISL_NULL;                 /* Matrix converting between basises */
  double *salloc = SISL_NULL;                  /* Matrix for memory allocation */
  double *salfa = SISL_NULL;                   /* The values of a discrete B-spline
                                             calculation */
  double *spek = SISL_NULL;                    /* Pointer used in traversing arrays */
  double *sb = SISL_NULL;                      /* Right hand side of equation */
  double *sfrom, *sto;
  double *sp;                             /* Help array for s1701 */
  double *stx = SISL_NULL;                     /* Knot vector after insertion of knots
                                             at start */
  int    *mpiv=SISL_NULL;                      /* Pointer to pivotation array */
  SISLCurve *kreg;                        /* k-regular curve */

  double *st = SISL_NULL;                      /* Internal version of et */
  double *scoef = SISL_NULL;                   /* Copy of the vertices of the surface */
  int    kdim = pc->idim;
  int    kk = pc->ik;
  int    kn = pc->in;

  double tlen;                            /* Length of curve parameterization */
  double tstart;                          /* New start knot value */
  double tend;                            /* New end knot value */
  int    ki, kj;
  int    knst;
  int    knstx;
  int    kleft = 0;
  int    kpl, kfi, kla;

  int    kstat;                           /* Error status from lower level */
  int    kpos = 0;



  /* Ensure reasonable return value */
  *rc = SISL_NULL;


  /* Test input */

  if ( kk < 1 ) goto err110;

  if ( kn < kk ) goto err111;

  if ( afak1 < DZERO || afak2 < DZERO )
  {
    /* Warning - so correct the factor(s) */
    *jstat = 1;

    if ( afak1 < DZERO ) afak1 = DZERO;
    if ( afak2 < DZERO ) afak2 = DZERO;
  }


  /* Ensure k-regular curve */

  make_cv_kreg(pc, &kreg, &kstat);
  if ( kstat < 0 ) goto error;


  /* Alloocate array for pivotation vector */

  mpiv = new0array(2*kk, INT);
  if ( mpiv == SISL_NULL ) goto err101;

  /* Allocate space (en bloc) for the local vectors */

  salloc = new0array(3*kn + 9*kk + 4*kk*kk + kdim*kn, DOUBLE);
  if ( salloc == SISL_NULL ) goto err101;

  ext = salloc;                    /* Size kn+kk */
  smatrix = ext + kn + kk;         /* Max size 4*kk*kk */
  salfa = smatrix + 4*kk*kk;       /* Size kk */
  scoef = salfa + kk;              /* Size kdim*kn */
  sb    = scoef + kdim*kn;         /* Size 2*kk */
  sp    = sb + 2*kk;               /* Size kk */
  st    = sp + kk;                 /* Size kn + 2*kk */
  stx   = st + kn + 2*kk;          /* Size kn + kk */



  /* Copy knots and vertices, to avoid destruction of curve */

  memcopy(ext, kreg->et, kn + kk, DOUBLE);
  memcopy(scoef, kreg->ecoef, kdim*kn, DOUBLE);


  /* Make extended knot vector */

  tlen = ext[kn] - ext[kk-1];

  if ( afak1 > DZERO )
  {
    /* Extend the basis at the start of the curve */

    tstart = ext[kk-1] - tlen*afak1;
    for ( ki = 0; ki < kk; ki++ )  ext[ki] = tstart;
  }

  if ( afak2 > DZERO )
  {
    /* Extend the basis at the end of the curve */

    tend = ext[kn] + tlen*afak2;
    for ( ki = kn; ki < kn+kk; ki++ )  ext[ki] = tend;
  }


  /* s1701 expects et to be a refinement of ext, thus we have to make a new
     version of e1 with the extra kk new knots before the start and after the
     end and one intermediate version with only kk at the start */

  memcopy(st, ext, kk - 1, DOUBLE);
  memcopy(st + kk - 1, kreg->et, kn + kk, DOUBLE);
  memcopy(st + 2*kk - 1 + kn, ext + kn + 1, kk - 1, DOUBLE);
  knst = kn + 2*(kk - 1);

  memcopy(stx, ext ,kn, DOUBLE);
  memcopy(stx + kn, st + kk - 1 + kn, 2*kk - 1, DOUBLE);
  knstx = kn + kk - 1;

  /* STEP 2 Make matrix going between bases, only the kk first and last knots
     are to be changed.  */


  /* Now we have two cases. We know that only the kk+1 first and kk+1
     last vertices are to be changed. However 2*(kk+1) might be equal to kn.
     Thus we have to change all vertices if kn<=2*(kk+1) */


  /* Make two steps one for the start and one for the end of the curve */


  /* Make matrix for the kk first vertices */

  for ( ki=kk-1, spek=smatrix; ki < 2*kk-1; ki++, spek+=kk )
  {
    /* we use kn instead of knstx since s1219 expects et[in-1] != et[in], we only
       address vertices at the start so this does not matter */

    s1219(stx, kk, kn, &kleft, st[ki], &kstat);
    if ( kstat < 0 ) goto error;

    s1701(ki, kleft, kk, knstx, &kpl, &kfi, &kla, st, stx, sp, salfa, &kstat);
    if( kstat < 0 ) goto error;

    /* Copy the discrete B-splines into the right position */

    memcopy(spek + kfi, salfa + kpl + kfi, kla - kfi + 1, DOUBLE);
  }



  /* Do the factorisation of the matrix */

  s6lufacp(smatrix, mpiv, kk, &kstat);
  if ( kstat < 0 ) goto error;

  /* The vertices in the curve are ordered in the sequence
     (x1,y1,z1),..,(xi,yi,zi), i=1,..,in.
     The only vertices affected by this backsubstitution is the kk first.
     We want to treat the back substitution as idim(=3) backsubstitutions.
     Thus we have to copy the proper parts of the vertices into a temporary
     array. Do backsubstitution and copy back into the curve object */


  for ( ki=0; ki < kdim; ki++ )
  {
    for ( kj=0,sfrom=(kreg->ecoef)+ki,sto=sb; kj < kk; kj++,sfrom+=kdim,sto++ )
      *sto = *sfrom;

    /* sb now contains the parts of vertices to be backsubsituted */

    s6lusolp(smatrix, sb, mpiv, kk, &kstat);
    if ( kstat < 0 ) goto error;

    /* Copy the backsubsituted vertices back into scoef */

    for ( kj=0,sto=scoef+ki,sfrom=sb;  kj < kk; kj++,sfrom++,sto+=kdim )
      *sto = *sfrom;
  }


  /* Make matrix for the kk last vertices */

  for ( ki=0,spek=smatrix; ki < kk*kk; ki++,spek++ ) *spek = DZERO;


  for ( ki=kn-kk,spek=smatrix; ki < kn; ki++,spek+=kk )
  {
    s1219(ext, kk, kn, &kleft, stx[ki], &kstat);
    if ( kstat < 0 ) goto error;

    s1701(ki, kleft, kk, kn, &kpl, &kfi, &kla, stx, ext, sp, salfa, &kstat);
    if ( kstat < 0 ) goto error;

    /* Copy the discrete B-splines into the right position */

    memcopy(spek+kfi-(kn-kk), salfa+kpl+kfi, kla-kfi+1, DOUBLE);
  }


  /* Do the factorisation of the matrix */

  s6lufacp(smatrix, mpiv, kk, &kstat);
  if ( kstat < 0 ) goto error;

  /* The vertices in the curve are ordered in the sequence
     (x1,y1,z1),..,(xi,yi,zi), i=1,..,in1. The only vertices
     affected by this backsubstitution is the kk-1 last rows.
     We want to treat the back substitution as idim(=3) backsubstitutions.
     Thus we have to copy the proper parts of the vertices into a temporary
     array. Do backsubstitution and copy back into the surface object */

  for ( ki=0; ki < kdim; ki++ )
  {
    for ( kj=0,sfrom=scoef+kdim*(kn-kk)+ki,sto=sb; kj < kk; kj++,sfrom+=kdim,sto++ )
      *sto = *sfrom;

    /* sb now contains the vertices to be backsubsituted */

    s6lusolp(smatrix, sb, mpiv, kk, &kstat);
    if ( kstat < 0 ) goto error;

    /* Copy the backsubsituted vertices back into scoef */

    for ( kj=0,sto=scoef+kdim*(kn-kk)+ki,sfrom=sb; kj < kk; kj++,sto+=kdim,sfrom++ )
      *sto = *sfrom;
  }


  /* Copy knots and vertices into the surface object */

  memcopy(kreg->ecoef, scoef, kdim*kn, DOUBLE);
  memcopy(kreg->et, ext, kn+kk, DOUBLE);

  /* Set periodicity flag */
  kreg->cuopen = SISL_CRV_OPEN;


  /* Task done */

  *rc = kreg;

  *jstat = 0;
  goto out;


  /* Error in allocation. */

 err101:
  *jstat = -101;
  s6err("s1233",*jstat,kpos);
  goto out;


  /* Error in curve description - order less than 1 */

 err110:
  *jstat = -110;
  s6err("s1233",*jstat,kpos);
  goto out;


  /* Error in curve desctiption - number of vertices is less than the order */

 err111:
  *jstat = -111;
  s6err("s1233",*jstat,kpos);
  goto out;


  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1233",*jstat,kpos);
  goto out;


 out:

  /* Free allocated scratch  */

  if (salloc != SISL_NULL) freearray(salloc);
  if (mpiv != SISL_NULL) freearray(mpiv);

  return;

}
