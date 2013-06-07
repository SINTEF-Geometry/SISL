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



#define S1941

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1941(SISLCurve *pcurve,int icont,int *jstat)
#else
void s1941(pcurve,icont,jstat)
     SISLCurve  *pcurve;
     int        icont;
     int    	*jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe the B-spline curve (i.e. not NURBS) pcurve with a
*              cyclic basis of continuity icont.
*
* INPUT      : pcurve - Pointer to the curve
*              icont  - The required continuity
*
* OUTPUT     : jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              pcurve - The curve with modified description
*
* METHOD     : 1. The cyclic knot vector with the right continuity is made
*              2. The transformation matrix for the ik2 first vertices between
*                 the cyclic and the old knot vector is made.
*              3. The transformation matrix is inverted and used to update the
*                 ik2 first vertices.
*              4. The transformation matrix for the ik2 last vertices between
*                 the cyclic and the old knot vector is made.
*              5. The transformation matrix is inverted and used to update the
*                 ik2 last vertices.
*              6. The knot vector is updated.
*-
* CALLS      : s1701,s6err.
*
* Written  by : Tor Dokken, SI, Oslo, Norway,feb-1992
* Renamed and adapted for curves by : Vibeke Skytt, SINTEF Oslo, 01.95.
*
*********************************************************************
*/
{
  double *scycl=SISL_NULL;                    /* Cyclic version of knot vector */
  double *smatrix=SISL_NULL;                   /* Matrix converting between baes */
  double *salloc=SISL_NULL;                    /* Matrix for memory allocation */
  double *salfa=SISL_NULL;                     /* The values of a discrete B-spline
                                             calculation */
  double *spek=SISL_NULL;                      /* Pointer used in traversing arrays */
  double *scoef=SISL_NULL;                     /* Copy of the vertices of the surface */
  double *sb=SISL_NULL;                        /* Right hand side of equation */
  double *sfrom,*sto;
  double *sp;                             /* Hlep array for s1701 */
  double *st1=SISL_NULL;                       /* Internal version of et */
  double *stx=SISL_NULL;                       /* Knot vector after insertion of knots
                                             at start */

  int    kdim = pcurve->idim;
  int    kk = pcurve->ik;
  int    kn = pcurve->in;

  int    kcont;                           /* Checked continuity */
  int    kmult;                           /* Multiplicity of knot kk-1 and kn*/
  int    ki,kj,kl;
  int    kperlength;
  int    kant;
  int    kleft=0;                         /* Pointer into knot vector */
  int    kpl,kfi,kla;                     /* Pointers into conversion matrix */
  int    kstat;
  int    *mpiv=SISL_NULL;                      /* Pointer to pivotation array */
  int    kpos = 0;
  int    knst1;                           /* NUmber of basis functions in st1 */
  int    knstx;                           /* Number of basis functions in stx */




  /* Test continuity */

  if (icont < 0) goto finished;
  kcont = icont;
  if (icont >= kk-2) icont = kk-2;

  /* Make multiplicty to be used at value et[ik-1] and et[in] */
  kmult = kk - kcont - 1;

  /* Make the number of knots to be changed at the start and the end, this
     is also equal to extra knot to be inserted in internal version of array et */

  kant = kk-kmult;


  /* Alloocate array for pivotation vector */

  mpiv = new0array(2*kk,INT);
  if (mpiv == SISL_NULL) goto err101;

  salloc = new0array(3*kn+9*kk+4*kk*kk+kdim*kn,DOUBLE);
  if (salloc == SISL_NULL) goto err101;
  scycl = salloc;                  /* Size kn+kk */
  smatrix = scycl + kn + kk;  /* Max size 4*kk*kk */
  salfa = smatrix + 4*kk*kk;     /* Size kk */
  scoef = salfa + kk;           /* Size kdim*kn */
  sb    = scoef + kdim*kn;    /* Size 2*kk */
  sp    = sb + 2*kk;              /* Size kk */
  st1   = sp + kk;                /* Size kn + 2*kk */
  stx   = st1 + kn + 2*kk;       /* Size kn + 2*kk */



  /* Copy vertices, to avoid destruction of the curve */

  memcopy(scoef,pcurve->ecoef,kdim*kn,DOUBLE);



  /* Make cyclic knot vector */


  /* First copy all knots */

  memcopy(scycl,pcurve->et,kn+kk,DOUBLE);

  /* The change the ik first and the ik last knots to make a cyclic basis */

  kperlength = kn - kk + kmult;

  /* Make knots 0 to ik - kmult - 1 */

  for (ki=kk-kmult-1 ; 0<=ki ; ki--)
    {
      scycl[ki] = scycl[kk-1] - (scycl[kn] - scycl[ki+kperlength]);
    }


  /* Make knots kn1 + kmult to kn1 + kk1 -1 */

  for (ki=kmult ; ki < kk ; ki++)
    {
      scycl[kn+ki] = scycl[kn] + (scycl[kk+ki-kmult] - scycl[kk-1]);

    }
      /* s1701 expects et1 to be a refinement of scyclic, thus we have to make a new
	 version of et1 with the extra kk1-kmult new knots before the start and
	 after the end and one intermediate version with only kk1-kmult at the start */

  memcopy(st1,scycl,kant,DOUBLE);
  memcopy(st1+kant,pcurve->et,kn+kk,DOUBLE);
  memcopy(st1+kant+kk+kn,scycl+kn+kk-kant,kant,DOUBLE);
  knst1 = kn + 2*kant;

  memcopy(stx,scycl,kn,DOUBLE);
  memcopy(stx+kn,st1+kant+kn,kk+kant,DOUBLE);
  knstx = kn + kant;

  /* STEP 2 Make matrix going between bases, only the kk2-kmult first and last knots
     are to be changed.  */


  /* Now we have two cases. We know that only the kk1-kmult first and kk1-kmult
     last vertices are to be changed. However 2*(kk1-kmult) might be a bigger
     number than kn1. Thus we have to change all vertices if kn1<=2(kk1-kmult) */


  /* Make two steps one for the start and one for the end of the surface */


  /* Make matrix for the kk1 first vertices */

  for (ki=kant,spek=smatrix ; ki <kk+kant ; ki++, spek+=kk)
    {
      /* we use kn instead of knstx since s1219 expects et[in-1] != et[in], we only
         address vertices at the start so this does not matter*/

      s1219(stx,kk,kn,&kleft,st1[ki],&kstat);
      if (kstat<0) goto error;

      s1701(ki,kleft,kk,knstx,&kpl,&kfi,&kla,st1,stx,sp,salfa,&kstat);
      if(kstat<0) goto error;

      /* Copy the discrete B-splines into the right position */

      memcopy(spek+kfi,salfa+kpl+kfi,kla-kfi+1,DOUBLE);
    }



  /* Do the factorisation of the matrix */

  s6lufacp(smatrix,mpiv,kk,&kstat);
  if (kstat<0) goto error;

  /* The vertices in the curve are ordered in the sequence
     (x1,y1,z1),..,(xi,yi,zi), i=1,..,in. The only vertices
     affected by this backsubstitution is the kant first ones.
     We want to treat the back substitution
     as idim(=3) backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the curve object */


  for (kl=0 ; kl<kdim ; kl++)
  {
     for (kj=0, sfrom=(pcurve->ecoef)+kl,sto=sb ;
      kj<kk ; kj++,sfrom+=kdim,sto++)
	*sto = *sfrom;
     
     /* sb now contains the vertices to be backsubsituted */
     
     s6lusolp(smatrix,sb,mpiv,kk,&kstat);
     if (kstat<0) goto error;
     
     /* Copy the backsubsituted vertices back into scoef */
     
     for (kj=0, sto=scoef+kl,sfrom=sb ;
      kj<kk ; kj++,sfrom++,sto+=kdim)
	*sto = *sfrom;
  }


  /* Make matrix for the kk last vertices */


  for (ki=0,spek=smatrix ; ki<kk*kk ; ki++,spek++) *spek = DZERO;


  for (ki=kn-kk ,spek=smatrix ; ki <kn ; ki++, spek+=kk)
    {
      s1219(scycl,kk,kn,&kleft,stx[ki],&kstat);
      if (kstat<0) goto error;

      s1701(ki,kleft,kk,kn,&kpl,&kfi,&kla,stx,scycl,sp,salfa,&kstat);
      if(kstat<0) goto error;

      /* Copy the discrete B-splines into the right position */

      memcopy(spek+kfi-(kn-kk),salfa+kpl+kfi,kla-kfi+1,DOUBLE);
    }



  /* Do the factorisation of the matrix */

  s6lufacp(smatrix,mpiv,kk,&kstat);
  if (kstat<0) goto error;

  /* The vertices in the curve are ordered in the sequence
     (x1,y1,z1),..,(xi,yi,zi), i=1,..,in. The only vertices
     affected by this backsubstitution is the kant last ones.
     We want to treat the back substitution
     as idim(=3) backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the curve object */

  for (kl=0 ; kl<kdim ; kl++)
  {
     for (kj=0, sfrom=scoef+kdim*(kn-kk)+kl,sto=sb ;
      kj<kk ; kj++,sfrom+=kdim,sto++)
	*sto = *sfrom;
     
     /* sb now contains the vertices to be backsubsituted */
     
     s6lusolp(smatrix,sb,mpiv,kk,&kstat);
     if (kstat<0) goto error;
     
     /* Copy the backsubsituted vertices back into scoef */
     
     for (kj=0, sto=scoef+kdim*(kn-kk)+kl,sfrom=sb ;
      kj<kk ; kj++,sto+=kdim,sfrom++)
	*sto = *sfrom;
  }


  /* Copy knots and vertices into the curve object */

  memcopy(pcurve->ecoef,scoef,kdim*kn,DOUBLE);
  memcopy(pcurve->et,scycl,kn+kk,DOUBLE);

  /* Set periodicity flag */
  pcurve->cuopen = SISL_CRV_PERIODIC;


  /* Task done */

 finished:

  *jstat = 0;
  goto out;

  /* Error in allocation. */

 err101:
  *jstat = -101;
  s6err("s1941",*jstat,kpos);
  goto out;



  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  s6err("s1941",*jstat,kpos);
  goto out;
 out:

  /* Free allocated scratch  */
  if (salloc != SISL_NULL) freearray(salloc);
  if (mpiv != SISL_NULL) freearray(mpiv);

  return;

}
