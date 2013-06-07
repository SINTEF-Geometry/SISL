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
 * $Id: mk_cv_cycl.c,v 1.3 2005-02-28 09:04:47 afr Exp $
 *
 */


#define MAKE_CV_CYCLIC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
make_cv_cyclic(SISLCurve *pcurve,int icont,int *jstat)
#else
void make_cv_cyclic(pcurve,icont,jstat)
     SISLCurve  *pcurve;
     int        icont;
     int    	*jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe the curve pcurve with a cyclic basis of continuity
*              icont.
*
* INPUT      : pcurve - Pointer to the curve
*              icont  - The required continuity
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              pcurve  - The curve with modified description
*
* METHOD     : 1. The cyclic knot vector with the right continuity is made
*              2. The transformation matrix for the ik first vertices between
*                 the cyclic and the old knot vector is made.
*              3. The transformation matrix is inverted and used to update the
*                 ik first vertices.
*              4. The transformation matrix for the ik last vertices between
*                 the cyclic and the old knot vector is made.
*              5. The transformation matrix is inverted and used to update the
*                 ik last vertices.
*              6. The knot vector is updated.
*-
* CALLS      : s1701,s6err.
*
* Written  by : Vibeke Skytt, SI, 05.92 based on a routine by
*               Tor Dokken, SI, Oslo, Norway,feb-1992
*
*********************************************************************
*/
{
  double *scycl=SISL_NULL;                    /* Cyclic version of knot vector */
  double *smatrix=SISL_NULL;                   /* Matrix converting between baes */
  /* Pointers to two conversion matrices  @afr: Commented out */
  /*  double *smatr1=SISL_NULL; */
  /*  double *smatr2=SISL_NULL; */ 
  double *salloc=SISL_NULL;                    /* Matrix for memory allocation */
  double *salfa=SISL_NULL;                     /* The values of a discrete B-spline
                                             calculation */ 
  double *spek=SISL_NULL;                      /* Pointer used in traversing arrays */
  double *scoef=SISL_NULL;                     /* Copy of the vertices of the surface */
  double *sb=SISL_NULL;                        /* Right hand side of equation */
  double *sfrom,*sto;
  double *sp;                             /* Help array for s1701 */  
  double *st1=SISL_NULL;                       /* Internal version of et */  
  double *stx=SISL_NULL;                       /* Knot vector after insertion of knots
                                             at start */ 
 
  int    kdim = pcurve->idim; 
  int    kk = pcurve->ik;
  int    kn = pcurve->in;
  int    rat = (pcurve->ikind == 1 || pcurve->ikind == 3) ? 0 : 1;
  double *sourcecoef = (rat) ? pcurve->rcoef : pcurve->ecoef;
  int    kdim2 = kdim + rat;
  
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
  if (icont >= kk-2) icont = kk-2;
  
  /* Make multiplicty to be used at value et[ik-1] and et[in] */
  kmult = kk - kcont - 1;
  
  /* Make the number of knots to be changed at the start and the end, this
     is also equal to extra knot to be inserted in internal version of array et */
  
  kant = kk-kmult;
  
  
  /* Alloocate array for pivotation vector */
  
  mpiv = new0array(2*kk,INT);
  if (mpiv == SISL_NULL) goto err101;
  
  salloc = new0array(3*kn+9*kk+4*kk*kk+kdim2*kn,DOUBLE);
  if (salloc == SISL_NULL) goto err101;
  scycl = salloc;                  /* Size kn+kk */
  smatrix = scycl + kn + kk;  /* Max size 4*kk*kk */
  salfa = smatrix + 4*kk*kk;     /* Size kk */
  scoef = salfa + kk;           /* Size kdim2*kn */
  sb    = scoef + kdim2*kn;    /* Size 2*kk */  
  sp    = sb + 2*kk;              /* Size kk */
  st1   = sp + kk;                /* Size kn + 2*kk */
  stx   = st1 + kn + 2*kk;       /* Size kn + 2*kk */
  
  
  
  /* Copy vertices, to avoid destruction of curve */
  
  memcopy(scoef, (rat) ? pcurve->rcoef : pcurve->ecoef, kdim2*kn, DOUBLE);
  
  
  
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
  
  
  /* Make knots kn + kmult to kn + kk -1 */
  
  for (ki=kmult ; ki < kk ; ki++)
    {
      scycl[kn+ki] = scycl[kn] + (scycl[kk+ki-kmult] - scycl[kk-1]);
      
    }
      /* s1701 expects et to be a refinement of scyclic, thus we have to make a new
	 version of et with the extra kk-kmult new knots before the start and
	 after the end and one intermediate version with only kk-kmult at the start */
 
  memcopy(st1,scycl,kant,DOUBLE);
  memcopy(st1+kant,pcurve->et,kn+kk,DOUBLE);
  memcopy(st1+kant+kk+kn,scycl+kn+kk-kant,kant,DOUBLE);
  knst1 = kn + 2*kant;

  memcopy(stx,scycl,kn,DOUBLE);
  memcopy(stx+kn,st1+kn+kant,kk+kant,DOUBLE);
  knstx = kn + kant;
  
  /* STEP 2 Make matrix going between bases, only the kk-kmult first and last 
     knots are to be changed.  */
  
  
  /* Now we have two cases. We know that only the kk-kmult first and kk-kmult
     last vertices are to be changed. However 2*(kk-kmult) might be a bigger
     number than kn. Thus we have to change all vertices if kn<=2(kk-kmult) */
  
  
  /* Make two steps one for the start and one for the end of the surface */
  
  
  /* Make matrix for the kk first vertices */
  
  for (ki=kant,spek=smatrix ; ki <kk+kant ; ki++, spek+=kk)
    {
      
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
  
  /* TThe only vertices of the curve
     affected by this backsubstitution is the kant first.
     We want to treat the back substitution
     as idim(=3) backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the curve object */
  
  
    for (kl=0 ; kl<kdim2 ; kl++)
      {
	for (kj=0, sfrom=(sourcecoef)+kl,sto=sb ;
	     kj<kk ; kj++,sfrom+=kdim2,sto++)
	  *sto = *sfrom;
	
	/* sb now contains the vertices to be backsubsituted */
	
	s6lusolp(smatrix,sb,mpiv,kk,&kstat);
	if (kstat<0) goto error;
	
	/* Copy the backsubsituted vertices back into scoef */
	
	for (kj=0, sto=scoef+kl,sfrom=sb ;
	     kj<kk ; kj++,sfrom++,sto+=kdim2)
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
  
  /* The only vertices of the curve
     affected by this backsubstitution is the kant last.
     We want to treat the back substitution
     as idim(=3) backsubstitutions. Thus we have to copy the proper
     parts of the vertices into a temporary array. Do backsubstitution and
     copy back into the curve object */
  
    for (kl=0 ; kl<kdim2 ; kl++)  
      {
	for (kj=0, sfrom=scoef+kdim2*(kn-kk)+kl,sto=sb ;
	     kj<kk ; kj++,sfrom+=kdim2,sto++)
	  *sto = *sfrom;
	
	/* sb now contains the vertices to be backsubsituted */
	
	s6lusolp(smatrix,sb,mpiv,kk,&kstat);
	if (kstat<0) goto error;
	
	/* Copy the backsubsituted vertices back into scoef */
	
	for (kj=0, sto=scoef+kdim2*(kn-kk)+kl,sfrom=sb ;
	     kj<kk ; kj++,sto+=kdim2,sfrom++)
	  *sto = *sfrom;
      }
  
  
    /* Copy knots and vertices into the curve object */

    memcopy((rat) ? pcurve->rcoef : pcurve->ecoef, scoef, kdim2*kn, DOUBLE);
    memcopy(pcurve->et,scycl,kn+kk,DOUBLE); 
    pcurve->cuopen = SISL_CRV_PERIODIC;

    /* Update divided coefficients */
    if (rat)
      {
	for (ki=0; ki<kn; ++ki)
	  {
	    for (kj=0; kj<kdim; ++kj)
	      pcurve->ecoef[ki*kdim+kj] = 
		pcurve->rcoef[ki*kdim2+kj]/pcurve->rcoef[ki*kdim2+kdim];
	  }
      }

  
  /* Task done */
  
 finished:
  
  *jstat = 0;
  goto out;
  
  /* Error in allocation. */
  
 err101: 
  *jstat = -101;
  s6err("make_cv_cyclic",*jstat,kpos);
  goto out;
  
  
  
  /* Error in lower level routine.  */
  
  error : 
    *jstat = kstat;     
  s6err("make_cv_cyclic",*jstat,kpos);
  goto out;
 out:
  
  /* Free allocated scratch  */
  if (salloc != SISL_NULL) freearray(salloc);  
  if (mpiv != SISL_NULL) freearray(mpiv);
  
  return;
  
}
