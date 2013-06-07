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
 * $Id: s1730.c,v 1.3 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1730

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1730(SISLCurve *pc,SISLCurve **rcnew,int *jstat)
#else
void s1730(pc,rcnew,jstat)
     SISLCurve *pc;
     SISLCurve **rcnew;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a B-spline curve to a sequence of Bezier
*              curves.
*
*
* INPUT      : pc        - SISLCurve to convert.
*
*
*
* OUTPUT     : rcnew     - The new Bezier represented curve .
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
*              newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              S1701.C   - Making the knot-inserten-transformation matrix.
*              make_cv_kreg - Ensure that the input curve is k-regular.
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, May 1992 (Inroduced NURBS)
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994. Moved free'ing
*              of 'qkreg' from label 'outfree' to label 'out'.
*
**********************************************************************/
{
  register int ki,ki1,ki2; /* Control variable in loop.                */
  register int kj,kj1,kj2; /* Control variable in loop.                */
  register double *s1;   /* Pointers used in loop.                     */
  int kstat;             /* Local status variable.                     */
  int kpos=0;            /* Position of error.                         */
  int kmy;               /* An index to the knot-vector.               */
  int kpl,kfi,kla;       /* To posisjon elements in trans.-matrix.     */
  int kk=pc->ik;         /* Order of the input curve.                  */
  int kn=pc->in;         /* Number of the vertices in input curves.    */
  int kdim=pc->idim;     /* Dimensjon of the space in whice curve lies.*/
  int kn1;               /* Number of vertices in the new curve .      */
  double *st=SISL_NULL;       /* The new knot-vector.                       */
  double *sp=SISL_NULL;       /* To use in s1701.c                          */
  double *salfa=SISL_NULL;    /* A line of the trans.-matrix.               */
  double *scoef=SISL_NULL;    /* The new vertice.                           */
  SISLCurve *qkreg = SISL_NULL;   /* Input curve made k-regular.            */
  SISLCurve *q1=SISL_NULL;    /* Pointer to new curve-object.               */

  double *rcoef;         /* Potential rational coefficients.           */

  /* Check that we have a curve to treat. */

  if (!pc) goto err150;

  if (pc->cuopen == SISL_CRV_PERIODIC)
  {
     make_cv_kreg(pc,&qkreg,&kstat);
     if (kstat < 0) goto err153;
  }
  else qkreg = pc;

  /* Check if this is a rational curve. */

  if (qkreg->ikind == 2 || qkreg->ikind == 4)
    {
       kdim++;
       rcoef = qkreg->rcoef;
    }
  else
    rcoef = qkreg->ecoef;

  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix, and space for new knots
     to use in s1701.c */

  if ((salfa=newarray(kk,double))==SISL_NULL) goto err101;
  if ((sp=newarray(kk,double))==SISL_NULL) goto err101;

  /* Find the number of vertices in the new curve. */

  for(ki=0,kn1=0;ki<kn+kk;ki+=kj,kn1+=kk)
    for(kj=1;ki+kj<kn+kk && (qkreg->et[ki] == qkreg->et[ki+kj]);kj++);
  kn1 -= kk;

  /* Allocating the new arrays to the new curve. */

  if ((st=newarray(kn1+kk,double))==SISL_NULL) goto err101;
  if ((scoef=new0array(kn1*kdim,double))==SISL_NULL) goto err101;

  /* Making the new knotvectors. */

  for(ki=0,ki1=0;ki<kn+kk;ki+=kj)
    {
      for(kj=1;ki+kj<kn+kk && (qkreg->et[ki] == qkreg->et[ki+kj]);kj++);
      for(kj1=0;kj1<kk;kj1++,ki1++) st[ki1] = qkreg->et[ki];
    }

  /* Updating the coefisientvector to the new curve.*/

  for(s1=scoef,ki=0,kmy=0;ki<kn1;ki++)
    {
      /* Here we compute a new line with line number ki of
	 the knot inserten matrix. */

      while(qkreg->et[kmy+1] <= st[ki]) kmy++;
      s1701(ki,kmy,kk,kn,&kpl,&kfi,&kla,st,qkreg->et,sp,salfa,&kstat);
      if (kstat) goto err153;

      /* Compute the kdim vertices with the same "index". */

      for (kj=0; kj<kdim; kj++,s1++)
	for (*s1=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla; kj1++,kj2++)
	  {
	    ki2=kj1*kdim+kj;
	    *s1 += salfa[kj2] * rcoef[ki2];
	  }
    }

  /* Allocating new curve-objects.*/

  if ((q1=newCurve(kn1,kk,st,scoef,qkreg->ikind,qkreg->idim,2)) == SISL_NULL) goto err101;

  /* Updating output. */

  *rcnew = q1;
  *jstat = 0;
  goto out;


  /* Error. Subrutine error. */

 err153: *jstat = kstat;
  goto outfree;


  /* Error. No curve to treat.  */

 err150: *jstat = -150;
  s6err("s1730",*jstat,kpos);
  goto out;



  /* Error. Allocation error, not enough memory.  */

 err101: *jstat = -101;
  s6err("s1730",*jstat,kpos);
  goto outfree;


 outfree:
  if(q1) freeCurve(q1);
  else
    {
      if (st) freearray(st);
      if (scoef) freearray(scoef);
    }


  /* Free local used memory. */

 out:
  if (qkreg != SISL_NULL && qkreg != pc) freeCurve(qkreg);
  if (salfa) freearray(salfa);
  if (sp) freearray(sp);
}
