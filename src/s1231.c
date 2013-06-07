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
 * $Id: s1231.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1231
#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1231(SISLCurve *pc1,double apar,SISLCurve **rcnew1,SISLCurve **rcnew2,int *jstat)
#else
void s1231(pc1,apar,rcnew1,rcnew2,jstat)
     SISLCurve  *pc1;
     double apar;
     SISLCurve  **rcnew1;
     SISLCurve  **rcnew2;
     int    *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Subdivide a B-spline curve at a given parameter-value.
*              NOTE: When the curve is periodic (ie. cuopen is
*                    set to SISL_CRV_PERIODIC and open basis
*                    with order-mult repeated knots and coeffic.)
*                    this function will return only ONE curve
*                    through rcnew1. This curve is the same
*                    geometric curve as pc1, but is represented on
*                    a closed basis with start and end at pc(apar).
*                    Cuopen is set to SISL_CRV_CLOSED.
*                    jstat equals 2 when this occurs.
*
*
*
* INPUT      : pc1      - SISLCurve to subdivide.
*              apar     - Parameter-value at which to subdivide.
*
*
*
* OUTPUT     : rcnew1    - First part of the subdivided curve.
*              rcnew2    - Second part of the subdivided curve.
*                          If the parameter value is at the end of a
*                          curve SISL_NULL pointers might be returned
*              jstat     - status messages
*                                         = 2      : pc1 periodic, rcnew2=SISL_NULL
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Inserting knots at the subdividing-point until
*              we have a ktuple-knot. Then we may separate the
*              curve into two parts.
*
*
* REFERENCES :
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              S1700.C   - Making the knot-inserten-transformation matrix.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* MODIFIED BY : Mike Floater, SI, 90-12. Subdivide rational curves.
* MODIFIED BY : Ulf J. Krystad, SI, 92-01. Subdivide periodic crvs.
* MODIFIED BY : Arne Laksaa, SI, 92-09. Move apar to closest knot
*		if it is close to a knot. Using D(N)EQUAL().
*
*
**********************************************************************/
{
  int kstat;              /* Local status variable.                     */
  int kpos=0;             /* Position of error.                         */
  int kmy;                /* An index to the knot-vector.               */
  int kv,kv1;             /* Number of knots we have to insert.         */
  int kpl,kfi,kla;        /* To posisjon elements in trans.-matrix.     */
  int kk=pc1->ik;         /* Order of the input curve.                  */
  int kn=pc1->in;         /* Number of the vertices in input curves.    */
  int kdim=pc1->idim;     /* Dimensjon of the space in whice curve lies.*/
  int kind=pc1->ikind;    /* Type of curve pc1 is.                      */
  int kn1,kn2;            /* Number of vertices in the new curves.      */
  int knum;               /* Number of knots less and equal than
			     the intersection point.                    */
  int ki,ki1;             /* Control variable in loop.                  */
  int kj,kj1,kj2;         /* Control variable in loop.                  */
  int newkind=1;          /* Type of curve the subcurves are            */
  double *s1,*s2,*s3,*s4; /* Pointers used in loop.                     */
  double *st1=SISL_NULL;       /* The first new knot-vector.                 */
  double *st2=SISL_NULL;       /* The second new knot-vector.                */
  double *salfa=SISL_NULL;     /* A line of the trans.-matrix.               */
  double *scoef;          /* Pointer to vertices.                       */
  double *scoef1=SISL_NULL;    /* The first new vertice.                     */
  double *scoef2=SISL_NULL;    /* The second new vertice.                    */
  SISLCurve *q1=SISL_NULL;     /* Pointer to new curve-object.               */
  SISLCurve *q2=SISL_NULL;     /* Pointer to new curve-object.               */
  int incr;		  /* No of extra knots copied during periodicity*/
  int mu;		  /* Multiplisity at the k'th knot              */
  int kleft = kk-1;	  /* Knot navigator                             */
  double delta;           /* Period size in knot array.                 */
  double salfa_local[5];  /* Local help array.			        */
  
  *rcnew1 = SISL_NULL;
  *rcnew2 = SISL_NULL;
  
  
 /* if pc1 is rational, do subdivision in homogeneous coordinates */
 /* just need to set up correct dim and kind for the new curves at end of routine */
  if(kind == 2 || kind == 4)
  {
     scoef = pc1->rcoef;
     kdim++;
     newkind++;
  }
  else
  {
     scoef = pc1->ecoef;
  }
  
  /* Check that we have a curve to subdivide. */
  
  if (!pc1) goto err150;
  
  /* Periodic curve treatment, UJK jan 92--------------------------------- */
  if (pc1->cuopen == SISL_CRV_PERIODIC)
  {
     delta = (pc1->et[kn] - pc1->et[kk - 1]);
     
     /* Check that the intersection point is an interior point. */
     if (apar < *(pc1->et) || apar > *(pc1->et + kn + kk - 1))
	goto err158;
     
     /* If inside the knot vector, but outside well define
	intervall, we shift the parameter value one period. */
     if (apar < *(pc1->et + kk - 1))
	apar += delta;
     if (apar > *(pc1->et + kn))
	apar -= delta;
     
     /* Now we create a new curve that is a copy of pc1,
	but with the period repeated once,
	this allows us to pick a whole period. */
     
     /* Get multiplisity at start of full basis interval */
     mu = s6knotmult(pc1->et, kk, kn, &kleft, pc1->et[kk-1], &kstat);
     if (kstat < 0) goto err153;
     if (mu >= kk) goto errinp;
     
     /* Copy ----------------------------------- */
     incr = kn - kk + mu;
     if ((scoef1 = newarray ((kn + incr) * kdim, double)) == SISL_NULL)
	goto err101;
     if ((st1 = newarray (kn + kk + incr, double)) == SISL_NULL)
	goto err101;
     
     memcopy (scoef1, pc1->ecoef, kn * kdim, double);
     memcopy (st1, pc1->et, kn + kk, double);
     memcopy (scoef1 + kn * kdim, pc1->ecoef + (kk - mu) * kdim, 
	      incr * kdim, double);
     
     
     for (ki = 0; ki < incr; ki++)
	st1[kn + kk + ki] = st1[kn + kk + ki - 1] +
	   (st1[2*kk - mu + ki] - st1[2*kk - mu + ki - 1]);
     if ((q1 = newCurve (kn + incr, kk, st1, scoef1,
			 newkind, pc1->idim, 2)) == SISL_NULL)
	goto err101;
     q1->cuopen = SISL_CRV_OPEN;
     
     /* Pick part (one period)------------------ */
     s1712 (q1, apar, apar + delta,
	    rcnew1, &kstat);
     if (kstat < 0)
	goto err153;
     freeCurve (q1);
     if (*rcnew1)
	(*rcnew1)->cuopen = SISL_CRV_CLOSED;
     
     /* Finished, exit */
     *jstat = 2;
     goto out;
     
  }
  
  /* End of periodic curve treatment, UJK jan 92------------- */

	    
  /* Check that the intersection point is an interior point. */
  /* Changed by UJK */
  /*if (apar <= *(pc1->et) || apar >= *(pc1->et+kn+kk-1)) goto err158; */
  if ((apar < pc1->et[kk - 1] || DEQUAL(apar, pc1->et[kk - 1]))||
      (apar > pc1->et[kn] || DEQUAL(apar, pc1->et[kn])))
    goto err158;
  
  
  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix.*/
  
  if (kk > 5)
  {
     if ((salfa = newarray (kk, double)) == SISL_NULL)	goto err101;
  }
  else salfa = salfa_local;
  
  
  /* Find the number of the knots which is smaller or like
     the intersection point, and how many knots we have to insert.*/
  
  s1 = pc1->et;
  kv = kk;     /* The maximum number of knots we may have to insert. */
  
  /* Using binear search*/
  kj1=0;
  kj2=kk+kn-1;
  knum = (kj1+kj2)/2;
  while (knum != kj1)
  {
     if ((s1[knum] < apar) && DNEQUAL (s1[knum], apar))
	kj1=knum; else kj2=knum;
     knum = (kj1+kj2)/2;
  }
  knum++;           /* The smaller knots. */
  
  while (DEQUAL (s1[knum], apar))
     /* The knots thats like the intersection point. */
  { 
     apar = s1[knum];
     knum++; 
     kv--;
  }
  
  
  /* Find the number of vertices in the two new curves. */
  
  kn1 = knum + kv - kk;
  kn2 = kn + kk - knum;
  
  
  
  /* Allocating the new arrays to the two new curves. */
  
  if (kn1>0)
  {
     if ((scoef1=newarray(kn1*kdim,double))==SISL_NULL) goto err101;
     if ((st1=newarray(kn1+kk,double))==SISL_NULL) goto err101;
  }
  if (kn2>0)
  {
     if ((scoef2=newarray(kn2*kdim,double))==SISL_NULL) goto err101;
     if ((st2=newarray(kn2+kk,double))==SISL_NULL) goto err101;
  }
  
  
  /* Copying the knotvectors, all but the intersection point from
     the old curve to the new curves */
  
  memcopy(st1,pc1->et,kn1,double);
  memcopy(st2+kk,pc1->et+knum,kn2,double);
  
  
  /* Updating the knotvectors by inserting a k-touple knot in
     the intersection point at each curve.*/
  
  for(s2=st1+kn1,s3=st2,s4=s3+kk; s3<s4; s2++,s3++) *s2 = *s3 = apar;
  
  
  /* Copying the coefisientvectors to the new curves.*/
  
  memcopy(scoef1,scoef,kdim*kn1,double);
  memcopy(scoef2,scoef+kdim*(knum-kk),kdim*kn2,double);
  
  
  /* Updating the coefisientvectors to the new curves.*/
  
  /* Updating the first curve. */
  knum -= kk - 1;
  for (ki=max(0,knum),kv1=max(0,-knum),s1=scoef1+ki*kdim; ki<kn1; ki++)
  {
     /* Initialising:
	knum = knum-kk+1, Index of the first vertice to change.
	ki = knum,        Index of the vertices we are going to
	change. Starting with knum, but if
	knum is negativ we start at zero.
	kv1 = 0,          Number if new knots between index ki
	and ki+kk. We are starting one below
	becase we are counting up before using
	it. If knum is negativ we are not
	starting at zero but at -knum.
	s1=scoef1+ki*kdim,SISLPointer at the first vertice to
	change. */
     
     
     /* Using the Oslo-algorithm to make a transformation-vector
	from the old vertices to one new vertice. */
     
     kmy=ki;
     s1700(kmy,kk,kn,++kv1,&kpl,&kfi,&kla,pc1->et,apar,salfa,&kstat);
     if (kstat) goto err153;
     
     
     /* Compute the kdim vertices with the same "index". */
     
     for (kj=0; kj<kdim; kj++,s1++)
	for (*s1=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla; kj1++,kj2++)
	   *s1 += salfa[kj2] * scoef[kj1*kdim+kj];
  }
  
  /* And the second curve. */
  
  for (ki1=min(kn1+kv-1,kn+kv),s1=scoef2; ki<ki1; ki++)
  {
     /* Initialising:
	ki1 = kn1+kv-1,   the index of the vertice next to the
	last vertice we have to change.
	If we do not have so many vertices,
	we have to use the index next to the
	last vertice we have, kn+kv.
	s1=scoef2         Pointer at the first vertice to
	change. */
     
     
     /* Using the Oslo-algorithm to make a transformation-vector
	from the old vertices to one new vertice. */
     
     s1700(kmy,kk,kn,kv1--,&kpl,&kfi,&kla,pc1->et,apar,salfa,&kstat);
     if (kstat) goto err153;
     
     
     /* Compute the kdim vertices with the same "index". */
     
     for (kj=0; kj<kdim; kj++,s1++)
	for (*s1=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla; kj1++,kj2++)
	   *s1 += salfa[kj2] * scoef[kj1*kdim+kj];
  }
  
  
  /* Allocating new curve-objects.*/
  /* use pc1->idim rather than kdim in case pc1 is rational  */
  
  if (kn1>0)
     if ((q1=newCurve(kn1,kk,st1,scoef1,newkind,pc1->idim,2)) == SISL_NULL)
								goto err101;
  if (kn2>0)
     if ((q2=newCurve(kn2,kk,st2,scoef2,newkind,pc1->idim,2)) == SISL_NULL)
								goto err101;
  
  
  /* Updating output. */
  
  *rcnew1 = q1;
  *rcnew2 = q2;
  *jstat = 0;
  goto out;
  
  
  /* Error. Error in low level routine. */
  
err153:
     *jstat = kstat;
  goto outfree;
  
  
  /* Error. Error in input */
errinp:
     *jstat = -154;
  goto outfree;
  
  /* Error. No curve to subdivide.  */
err150:
     *jstat = -150;
  s6err("s1231",*jstat,kpos);
  goto out;
  
  
  /* Error. The parameter value is outside the curve.  */
  
err158:
     *jstat = -158;
  s6err("s1231",*jstat,kpos);
  goto out;
  
  
  /* Error. Allocation error, not enough memory.  */
  
err101:
     *jstat = -101;
  s6err("s1231",*jstat,kpos);
  goto outfree;
  
  
outfree:
   if(q1) freeCurve(q1);
   if(q2) freeCurve(q2);   
   
   /* Free local used memory. */
   
out:
   if(!q1)
   {
      if (st1) freearray(st1);
      if (scoef1) freearray(scoef1);
   }
   
   if(!q2)
   {
      if (st2) freearray(st2);
      if (scoef2) freearray(scoef2);
   }
   
   if (kk > 5 && salfa)
      freearray (salfa);
   return;
}
