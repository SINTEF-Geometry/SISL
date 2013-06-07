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
 * $Id: s1731.c,v 1.3 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1731

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1731(SISLSurf *ps,SISLSurf **rsnew,int *jstat)
#else
void s1731(ps,rsnew,jstat)
     SISLSurf *ps;
     SISLSurf **rsnew;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a B-spline surface to Bezier surfaces.
*
*
* INPUT      : ps     - Surface to convert.
*
*
*
* OUTPUT     : rsnew     - The new Bezier represented surface.
*              jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Inserting knots until all knots
*              have multiplisity pc->ik.
*
*
* REFERENCES :
*
*-
* CALLS      : newarray  - Allocate space for array of given type.
*              new0array - Allocate space whith zero values.
*              freearray - Free space occupied by given array.
*              newSurf   - Allocate space for a new surf-object.
*              freeSurf  - Free space occupied by given surf-object.
*              S1701.C   - Making the knot-inserten-transformation matrix.
*              make_sf_kreg - Ensure that the input surface is k-regular.
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, May 1992 (Introduced NURBS).
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994. Moved free'ing
*              of 'qkreg' from 'outfree:' to 'out:' to remove memory leak.
*
**********************************************************************/
{
  register int ki,ki1,ki2;       /* Control variable in loop.                   */
  register int kj,kj1,kj2,kj3;   /* Control variable in loop.                   */
  int kstat;                     /* Local status variable.                      */
  int kpos=0;                    /* Position of error.                          */
  int kmy;                       /* An index to the knot-vector.                */
  int kpl,kfi,kla;               /* To posisjon elements in trans.-matrix.      */
  int kk1=ps->ik1;               /* Order one of the input surface.             */
  int kk2=ps->ik2;               /* Order two of the input surface.             */
  int kn=ps->in1;                /* Number of the vertices in input curves.     */
  int kdim=ps->idim;             /* Dimensjon of the space in whice surf lies.  */
  int kn1,kn2;                   /* Number of vertices in the new surface.      */
  double *s1,*s2,*s3;            /* Pointers used in loop.                      */
  double *st1=SISL_NULL;              /* The new knot-vector.                        */
  double *st2=SISL_NULL;              /* The new knot-vector.                        */
  double *sp=SISL_NULL;               /* To use in s1701.c                           */
  double *salfa=SISL_NULL;            /* A line of the trans.-matrix.                */
  double *scoef=SISL_NULL;            /* The new vertice.                            */
  double *scoefh=SISL_NULL;           /* A new vertice for help.                     */
  SISLSurf *q1=SISL_NULL;             /* Pointer to new surf-object.                 */

  double *rcoef;                 /* Potential rational vertices.                */
  int rdim;                      /* Potential rational dimension.               */
  SISLSurf *qkreg=SISL_NULL;          /* Input surface made k-regular.               */

  /* Check that we have a surface to treat. */

  if (!ps) goto err150;

  /* Make sure that the surface is k-regular.  */

  if (ps->cuopen_1 == SISL_SURF_PERIODIC ||
      ps->cuopen_2 == SISL_SURF_PERIODIC)
  {
     make_sf_kreg(ps,&qkreg,&kstat);
     if (kstat < 0) goto err153;
  }
  else qkreg = ps;


  /* Check if the surface is rational. */

  if (qkreg->ikind == 2 || qkreg->ikind == 4)
  {
     rcoef = qkreg->rcoef;
     rdim = kdim + 1;
  }
  else
  {
     rcoef = qkreg->ecoef;
     rdim = kdim;
  }

  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix, and space for new knots
     to use in s1701.c */

  if ((salfa=newarray(kk1+kk2,double))==SISL_NULL) goto err101;
  if ((sp=newarray(kk1+kk2,double))==SISL_NULL) goto err101;

  /* Find the number of vertices in the first direction
     in the new surface. */

  for(ki=0,kn1=0;ki<kn+kk1;ki+=kj,kn1+=kk1)
    for(kj=1;ki+kj<kn+kk1 && (qkreg->et1[ki] == qkreg->et1[ki+kj]);kj++);
  kn1 -= kk1;

  /* Find the number of vertices in the second direction
     in the new surface. */

  for(kn=qkreg->in2,ki=0,kn2=0;ki<kn+kk2;ki+=kj,kn2+=kk2)
    for(kj=1;ki+kj<kn+kk2 && (qkreg->et2[ki] == qkreg->et2[ki+kj]);kj++);
  kn2 -= kk2;

  /* Allocating the new arrays to the new surface. */

  if ((st1=newarray(kn1+kk1,double))==SISL_NULL) goto err101;
  if ((st2=newarray(kn2+kk2,double))==SISL_NULL) goto err101;
  if ((scoefh=new0array(kn1*kn*rdim,double))==SISL_NULL) goto err101;
  if ((scoef=new0array(kn1*kn2*rdim,double))==SISL_NULL) goto err101;

  /* Making the new knotvectors in the first direction */

  for(kn=qkreg->in1,ki=0,ki1=0;ki<kn+kk1;ki+=kj)
    {
      for(kj=1;ki+kj<kn+kk1 && (qkreg->et1[ki] == qkreg->et1[ki+kj]);kj++);
      for(kj1=0;kj1<kk1;kj1++,ki1++) st1[ki1] = qkreg->et1[ki];
    }

  /* Making the new knotvectors in the second direction. */

  for(kn=qkreg->in2,ki=0,ki1=0;ki<kn+kk2;ki+=kj)
    {
      for(kj=1;ki+kj<kn+kk2 && (qkreg->et2[ki] == qkreg->et2[ki+kj]);kj++);
      for(kj1=0;kj1<kk2;kj1++,ki1++) st2[ki1] = qkreg->et2[ki];
    }

  /* Updating the coefisientvector to the new surface in the first
     direction.*/

  for(s1=scoefh,ki2=kn1*kn*rdim,ki=0,kmy=0;ki<kn1;ki++)
    {
      /* Here we compute a new line with line number ki of
	 the knot inserten matrix. */

      while(qkreg->et1[kmy+1] <= st1[ki]) kmy++;
      s1701(ki,kmy,kk1,qkreg->in1,&kpl,&kfi,&kla,st1,qkreg->et1,sp,salfa,&kstat);
      if (kstat) goto err153;

      /* Compute the kn2*rdim vertices with the same "index". */

      for (kj=0; kj<rdim; kj++,s1++)
	for (s2=s1,s3=s2+ki2,kj3=kj;s2<s3;s2+=rdim*kn1,kj3+=rdim*qkreg->in1)
	  for (*s2=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla; kj1++,kj2++)
	    *s2 += salfa[kj2] * rcoef[kj1*rdim+kj3];
    }

  /* Updating the coefisientvector to the new surface in the second
     direction.*/

  for(s1=scoef,ki2=kn1*rdim,ki=0,kmy=0;ki<kn2;ki++,s1+=kn1*rdim)
    {
      /* Here we compute a new line with line number ki of
	 the knot inserten matrix. */

      while(qkreg->et2[kmy+1] <= st2[ki]) kmy++;
      s1701(ki,kmy,kk2,kn,&kpl,&kfi,&kla,st2,qkreg->et2,sp,salfa,&kstat);
      if (kstat) goto err153;

      /* Compute the kn1*rdim vertices with the same "index". */

      for (kj=0; kj<rdim; kj++)
	for (s2=s1+kj,s3=s2+ki2,kj3=kj;s2<s3;s2+=rdim,kj3+=rdim)
	  for (*s2=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla; kj1++,kj2++)
	    *s2 += salfa[kj2] * scoefh[kj1*kn1*rdim+kj3];
    }


  /* Allocating new surface-objects.*/

  if ((q1=newSurf(kn1,kn2,kk1,kk2,st1,st2,scoef,qkreg->ikind,kdim,2)) == SISL_NULL)
    goto err101;

  /* Updating output. */

  *rsnew = q1;
  *jstat = 0;
  goto out;

  /* Error. Subrutine error. */

 err153: *jstat = kstat;
  goto outfree;

  /* Error. No surface to treat.  */

 err150: *jstat = -150;
  s6err("s1731",*jstat,kpos);
  goto out;

  /* Error. Allocation error, not enough memory.  */

 err101: *jstat = -101;
  s6err("s1731",*jstat,kpos);
  goto outfree;

 outfree:
  if(q1) freeSurf(q1);
  else
    {
      if (st1)   freearray(st1);
      if (st2)   freearray(st2);
      if (scoef) freearray(scoef);
    }

  /* Free local used memory. */

 out:
  if (qkreg != SISL_NULL && qkreg != ps) freeSurf(qkreg);
  if (salfa)  freearray(salfa);
  if (sp)     freearray(sp);
  if (scoefh) freearray(scoefh);
}
