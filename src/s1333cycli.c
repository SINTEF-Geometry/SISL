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
 * $Id: s1333cycli.c,v 1.4 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1333_CYCLIC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1333_cyclic(SISLSurf *vsurf,int icont,int *jstat)
#else
void s1333_cyclic(vsurf,icont,jstat)
     SISLSurf   *vsurf;
     int        icont;
     int    	*jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe the B-spline surface (i.e. not NURBS) vsurf with a
*              cyclic basis of continuity icont in the first parameter direction.
*
* INPUT      : vsurf  - Pointer to the surface
*              icont  - The required continuity
*
* OUTPUT     : jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              vsurf  - The surface with modified description
*
* METHOD     : 1. The cyclic knot vector with the right continuity is made
*              2. The transformation matrix for the ik2 first vertices between
*                 the cyclic and the old knot vector in second paramter direction
*                 is made.
*              3. The transformation matrix is inverted and used to update the
*                 ik2 first vertices.
*              4. The transformation matrix for the ik2 last vertices between
*                 the cyclic and the old knot vector in second paramter direction
*                 is made.
*              5. The transformation matrix is inverted and used to update the
*                 ik2 last vertices.
*              6. The knot vector is updated.
*-
* CALLS      : s1701,s6err.
*
* Written  by : Tor Dokken, SI, Oslo, Norway,feb-1992
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
  double *st1=SISL_NULL;                       /* Internal version of et1 */
  double *stx=SISL_NULL;                       /* Knot vector after insertion of knots
                                             at start */

  int    kdim = vsurf->idim;
  int    rat = (vsurf->ikind == 1 || vsurf->ikind == 3) ? 0 : 1;
  double *sourcecoef = (rat) ? vsurf->rcoef : vsurf->ecoef;
  int    kdim2 = kdim + rat;
  int    kk1 = vsurf->ik1;
  int    kn1 = vsurf->in1;
  int    kn2 = vsurf->in2;
  int    knumb = kdim2*kn1;

  int    kcont;                           /* Checked continuity */
  int    kmult;                           /* Multiplicity of knot kk2-1 and kn2*/
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
  if (icont >= kk1-2) icont = kk1-2;

  /* Make multiplicty to be used at value et1[ik1-1] and et1[in1] */
  kmult = kk1 - kcont - 1;

  /* Make the number of knots to be changed at the start and the end, this
     is also equal to extra knot to be inserted in internal version of array et */

  kant = kk1-kmult;


  /* Alloocate array for pivotation vector */

  mpiv = new0array(2*kk1,INT);
  if (mpiv == SISL_NULL) goto err101;

  salloc = new0array(3*kn1+9*kk1+4*kk1*kk1+kdim2*kn1*kn2,DOUBLE);
  if (salloc == SISL_NULL) goto err101;
  scycl = salloc;                  /* Size kn1+kk1 */
  smatrix = scycl + kn1 + kk1;  /* Max size 4*kk1*kk1 */
  salfa = smatrix + 4*kk1*kk1;     /* Size kk1 */
  scoef = salfa + kk1;           /* Size kdim2*kn1*kn2 */
  sb    = scoef + kdim2*kn1*kn2;    /* Size 2*kk1 */
  sp    = sb + 2*kk1;              /* Size kk1 */
  st1   = sp + kk1;                /* Size kn1 + 2*kk1 */
  stx   = st1 + kn1 + 2*kk1;       /* Size kn1 + 2*kk1 */



  /* Copy vertices, to avoid destruction of surface */

  memcopy(scoef,sourcecoef,kdim2*kn1*kn2,DOUBLE);



  /* Make cyclic knot vector */


  /* First copy all knots */

  memcopy(scycl,vsurf->et1,kn1+kk1,DOUBLE);

  /* The change the ik1 first and the ik1 last knots to make a cyclic basis */

  kperlength = kn1 - kk1 + kmult;

  /* Make knots 0 to ik - kmult - 1 */

  for (ki=kk1-kmult-1 ; 0<=ki ; ki--)
    {
      scycl[ki] = scycl[kk1-1] - (scycl[kn1] - scycl[ki+kperlength]);
    }


  /* Make knots kn1 + kmult to kn1 + kk1 -1 */

  for (ki=kmult ; ki < kk1 ; ki++)
    {
      scycl[kn1+ki] = scycl[kn1] + (scycl[kk1+ki-kmult] - scycl[kk1-1]);

    }
      /* s1701 expects et1 to be a refinement of scyclic, thus we have to make a new
	 version of et1 with the extra kk1-kmult new knots before the start and
	 after the end and one intermediate version with only kk1-kmult at the start */

  memcopy(st1,scycl,kant,DOUBLE);
  memcopy(st1+kant,vsurf->et1,kn1+kk1,DOUBLE);
  memcopy(st1+kant+kk1+kn1,scycl+kn1+kk1-kant,kant,DOUBLE);
  knst1 = kn1 + 2*kant;

  memcopy(stx,scycl,kn1,DOUBLE);
  memcopy(stx+kn1,st1+kant+kn1,kk1+kant,DOUBLE);
  knstx = kn1 + kant;

  /* STEP 2 Make matrix going between bases, only the kk2-kmult first and last knots
     are to be changed.  */


  /* Now we have two cases. We know that only the kk1-kmult first and kk1-kmult
     last vertices are to be changed. However 2*(kk1-kmult) might be a bigger
     number than kn1. Thus we have to change all vertices if kn1<=2(kk1-kmult) */


  /* Make two steps one for the start and one for the end of the surface */


  /* Make matrix for the kk1 first vertices */

  for (ki=kant,spek=smatrix ; ki <kk1+kant ; ki++, spek+=kk1)
    {
      /* we use kn1 instead of knstx since s1219 expects et[in-1] != et[in], we only
         address vertices at the start so this does not matter*/

      s1219(stx,kk1,kn1,&kleft,st1[ki],&kstat);
      if (kstat<0) goto error;

      s1701(ki,kleft,kk1,knstx,&kpl,&kfi,&kla,st1,stx,sp,salfa,&kstat);
      if(kstat<0) goto error;

      /* Copy the discrete B-splines into the right position */

      memcopy(spek+kfi,salfa+kpl+kfi,kla-kfi+1,DOUBLE);
    }



  /* Do the factorisation of the matrix */

  s6lufacp(smatrix,mpiv,kk1,&kstat);
  if (kstat<0) goto error;

  /* The vertices in the surface are ordered in the sequence
     (x11,y11,z11),..,(xij,yij,zij), i=1,..,in1 and j=1,..,in2.
     i is running faster than j. The only vertices
     affected by this backsubstitution is the kant first rows.
     We want to treat the back substitution
     as idim(=3)*in2 backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the surface object */


  for (ki=0 ; ki<kn2 ; ki++)
    for (kl=0 ; kl<kdim2 ; kl++)
      {
	for (kj=0, sfrom=sourcecoef+ki*knumb+kl,sto=sb ;
	     kj<kk1 ; kj++,sfrom+=kdim2,sto++)
	  *sto = *sfrom;

	/* sb now contains the vertices to be backsubsituted */

	s6lusolp(smatrix,sb,mpiv,kk1,&kstat);
	if (kstat<0) goto error;

	/* Copy the backsubsituted vertices back into scoef */

	for (kj=0, sto=scoef+ki*knumb+kl,sfrom=sb ;
	     kj<kk1 ; kj++,sfrom++,sto+=kdim2)
	  *sto = *sfrom;
      }


  /* Make matrix for the kk1 last vertices */


  for (ki=0,spek=smatrix ; ki<kk1*kk1 ; ki++,spek++) *spek = DZERO;


  for (ki=kn1-kk1 ,spek=smatrix ; ki <kn1 ; ki++, spek+=kk1)
    {
      s1219(scycl,kk1,kn1,&kleft,stx[ki],&kstat);
      if (kstat<0) goto error;

      s1701(ki,kleft,kk1,kn1,&kpl,&kfi,&kla,stx,scycl,sp,salfa,&kstat);
      if(kstat<0) goto error;

      /* Copy the discrete B-splines into the right position */

      memcopy(spek+kfi-(kn1-kk1),salfa+kpl+kfi,kla-kfi+1,DOUBLE);
    }



  /* Do the factorisation of the matrix */

  s6lufacp(smatrix,mpiv,kk1,&kstat);
  if (kstat<0) goto error;

  /* The vertices in the surface are ordered in the sequence
     (x11,y11,z11),..,(xij,yij,zij), i=1,..,in1 and j=1,..,in2.
     i is running faster than j The only vertices
     affected by this backsubstitution is the kant last rows.
     We want to treat the back substitution
     as idim(=3)*in2 backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the surface object */

  for (ki=0 ; ki<kn2 ; ki++)
    for (kl=0 ; kl<kdim2 ; kl++)
      {
	for (kj=0, sfrom=scoef+ki*knumb+kdim2*(kn1-kk1)+kl,sto=sb ;
	     kj<kk1 ; kj++,sfrom+=kdim2,sto++)
	  *sto = *sfrom;

	/* sb now contains the vertices to be backsubsituted */

	s6lusolp(smatrix,sb,mpiv,kk1,&kstat);
	if (kstat<0) goto error;

	/* Copy the backsubsituted vertices back into scoef */

	for (kj=0, sto=scoef+ki*knumb+kdim2*(kn1-kk1)+kl,sfrom=sb ;
	 kj<kk1 ; kj++,sto+=kdim2,sfrom++)
	  *sto = *sfrom;
      }


  /* Copy knots and vertices into the surface object */

  memcopy(sourcecoef,scoef,kdim2*kn1*kn2,DOUBLE);
  memcopy(vsurf->et1,scycl,kn1+kk1,DOUBLE);

  /* Set periodicity flag */
  vsurf->cuopen_1 = SISL_SURF_PERIODIC;

    /* Update divided coefficients */
    if (rat)
      {
	for (ki=0; ki<kn1*kn2; ++ki)
	  {
	    for (kj=0; kj<kdim; ++kj)
	      vsurf->ecoef[ki*kdim+kj] = 
		vsurf->rcoef[ki*kdim2+kj]/vsurf->rcoef[ki*kdim2+kdim];
	  }
      }

  /* Task done */

 finished:

  *jstat = 0;
  goto out;

  /* Error in allocation. */

 err101:
  *jstat = -101;
  s6err("s1333_cyclic",*jstat,kpos);
  goto out;



  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  s6err("s1333_cyclic",*jstat,kpos);
  goto out;
 out:

  /* Free allocated scratch  */
  if (salloc != SISL_NULL) freearray(salloc);
  if (mpiv != SISL_NULL) freearray(mpiv);

  return;

}
