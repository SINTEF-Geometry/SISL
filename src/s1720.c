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
 * $Id: s1720.c,v 1.3 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1720

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1720(SISLCurve *pc,int ider,SISLCurve **rcnew,int *jstat)
#else
void s1720(pc,ider,rcnew,jstat)
     SISLCurve *pc;
     int   ider;
     SISLCurve **rcnew;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To express the ider-th derivative of a B-pline curve
*              as a B-spline curve with order ider less than the
*              orginal B-spline curve.
*
*
*
* INPUT      : pc       - SISLCurve to make a derivative from.
*              ider     - The derevative to be produced. 0 <= ider < pc->ik.
*
*
*
* OUTPUT     : rcnew    - The result of the ider differentiation of pc.
*              jstat    - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : For each ider differentiation we go through the vertices and
*              compute new vertices with use of a recursion formula.
*              At last the curve is checked removing redundant knots and
*              vertices.
*
*
* REFERENCES : Larry L. Schumaker, Spline Functions:Basic Theory.  p.195.
*              Carl de Boor, A Practial Guide to Spline.           p.139.
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              s1705.c   - Remove unnecesery knots and vertices from a curve.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, May 1992 (Added the NURBS handling.
*              Remark that the ider-th derivate of a NURBS curve of order
*              k has order (ider + 1)*k - ider. Remark also that we have to
*              use a completely different method for NURBS).
* REVISED BY : Christophe Birkeland, SI, 92-07 (remove unnecessary
*              knots in this routine (not in s1905), also testing if
*              ider=0).
* REVISED BY : Christophe Rene Birkeland, SINTEF, May 1993.
*              *jstat = 0   in start of routine.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994.  Changed icopy
*              flag from '1' to '2' for 'qc1' in the NURBS case - to reduce
*              memory needs and remove memory leak from 'scoef' and 'st'.
*              NOTE: this routine doesn't handle closed case correctly and
*              it also doesn't handle the periodic NURBS case correctly.
*
**********************************************************************/
{
  int kstat;             /* Local status variable.                  */
  int kpos=0;            /* Error posision.                         */
  int kk;                /* Order of the output curve.              */
  int kn;                /* Number of the vertices in output curves.*/
  int kdim=pc->idim;     /* Dimension of the space in which
			    the curves lies.                        */
  int kj;                /* Control variable in loop.               */
  double tdel;           /* Help variabel.                          */
  double *s1,*s2,*s3;    /* Pointers used in loop.                  */
  double *st=SISL_NULL;       /* The first new knot-vector.              */
  double *scoef=SISL_NULL;    /* The first new vertice.                  */
  SISLCurve *q1=SISL_NULL;    /* Pointer to new curve-object.            */

  /* NURBS variables: */
  SISLCurve *rat=SISL_NULL;   /* The denominator curve.                  */
  double *ratcoef=SISL_NULL;  /* The vertices of rat.                    */
  int rdim;              /* Rational dimension.                     */
  double eps;            /* Knot equality resolution.               */
  int multadd;           /* Added multiplicity of interior knots.   */
  int mult;              /* Original multiplicity of the knots.     */
  int oldprev;           /* Previous index in the original knots.   */
  int oldcurr;           /* Current index in the original knots.    */
  int newcurr;           /* Current index in the new knots.         */
  double *par = SISL_NULL;    /* Parameter values used for interpolation */
  int *der = SISL_NULL;       /* The derivative indicators (= 0).        */
  double *deriv = SISL_NULL;  /* The derivates returned by s1221.        */
  double *tau = SISL_NULL;    /* Interpolation points.                   */
  double denom;          /* The denominator.                        */
  int left = 0;          /* Interval indicator.                     */
  int ki;                /* Index in for loop.                      */

  /* Initialization */

  *jstat = 0;

  /* Check that we have a curve to differentiable. */

  if (!pc) goto err150;

  /* Check that the number of derivation is legal. */

  if (ider < 0) goto err156;

  /* If ider = 0; just copy input-curve to output and return */

  if (ider == 0)
  {
     if (pc->ikind == 2 || pc->ikind == 4)
     {
	if ((*rcnew=newCurve(pc->in,pc->ik,pc->et,pc->rcoef,
			    pc->ikind,pc->idim,1)) == SISL_NULL)
	   goto err101;
     }
     else
     {
	if ((*rcnew=newCurve(pc->in,pc->ik,pc->et,pc->ecoef,
			    pc->ikind,pc->idim,1)) == SISL_NULL)
	   goto err101;
     }
     goto out;
  }

  if (pc->ikind == 2 || pc->ikind == 4)
  {

     /* NURBS curve. */

     rdim = kdim + 1;

     /* Make the denominator curve. */

     ratcoef = newarray(pc->in, DOUBLE);
     for (kj=0; kj<pc->in; kj++)
	ratcoef[kj] = pc->rcoef[(kj + 1)*rdim - 1];
     rat = newCurve(pc->in, pc->ik, pc->et, ratcoef, 1, 1, 1);

     /* Make resolution for testing of knot equality. */

     eps = fabs(pc->et[pc->in] - pc->et[pc->ik - 1])*REL_PAR_RES;

     /* Make the new knot vector. */

     kn = 0;
     kk = (ider + 1)*pc->ik - ider;
     multadd = ider*pc->ik;
     st = newarray((2 + pc->in - pc->ik)*kk, DOUBLE);

     for (newcurr=0; newcurr<kk; newcurr++)
     {
	st[newcurr] = pc->et[pc->ik - 1];
        kn++;
     }

     oldcurr = pc->ik;
     oldprev = oldcurr;
     while ((pc->et[oldcurr] + eps) < pc->et[pc->in])
     {
	mult = 0;
	while ((pc->et[oldcurr] - pc->et[oldprev]) < eps)
	{
	   oldcurr++;
	   mult++;
	}
	mult += multadd;
	if (mult > kk) mult = kk;
	for (kj=0; kj<mult; kj++)
	{
	   st[newcurr + kj] = pc->et[oldprev];
	   kn++;
	}
	newcurr += mult;
	oldprev = oldcurr;
     }

     for (kj=0; kj<kk; kj++)
	st[newcurr + kj] = pc->et[pc->in];

     /* Calculate parameter values and derivate indicators. */

     s1890(st, kk, kn, &par, &der, &kstat);
     if (kstat < 0) goto error;

     /* Calculate interpolation points. */

     deriv = newarray((ider + 1)*kdim, DOUBLE);
     tau = newarray(kn*rdim, DOUBLE);

     for (kj=0; kj<kn; kj++)
     {
	s1221(rat, 0, par[kj], &left, &denom, &kstat);
	if (kstat < 0) goto error;
	denom = pow(denom, (ider + 1));
	s1221(pc, ider, par[kj], &left, deriv, &kstat);
	if (kstat < 0) goto error;
	for (ki=0; ki<kdim; ki++)
           tau[kj*rdim + ki] = deriv[ider*kdim + ki]*denom;
        tau[kj*rdim + kdim] = denom;
     }

     /* Make the new curve description. */

     s1891(par, tau, rdim, kn, 1, der, TRUE, st, &scoef, &kn,
	   kk, 0, 0, &kstat);
     if (kstat < 0) goto error;

     q1 = newCurve(kn, kk, st, scoef, pc->ikind, pc->idim, 2);
     /* (Note, copy == 2, so don't free 'st' or 'scoef' on exit). */

     /* Free allocated geometry. */

     if (rat != SISL_NULL) freeCurve(rat);
     rat = SISL_NULL;
     if (ratcoef != SISL_NULL) freearray(ratcoef);
     ratcoef = SISL_NULL;
     if (par != SISL_NULL) freearray(par);
     par = SISL_NULL;
     if (der != SISL_NULL) freearray(der);
     der = SISL_NULL;
     if (deriv != SISL_NULL) freearray(deriv);
     deriv = SISL_NULL;
     if (tau != SISL_NULL) freearray(tau);
     tau = SISL_NULL;

  }
  else
  {
     /* Not NURBS. */

     /* Find the number of vertices kn and the order kk for
        the derivative curve. */

     if (ider >= pc->ik)
       {
         kn = pc->in + pc->ik -1;
         kk = 1;
       }
     else
       {
         kn = pc->in + ider;
         kk = pc->ik - ider;
       }

     /* Allocating the new arrays to the new curve. */

     if ((st=newarray(kn+kk,double))==SISL_NULL) goto err101;
     if ((scoef=new0array(kn*kdim,double))==SISL_NULL) goto err101;

     /* Copying the knot vectors from the old curve to the new curve. */

     memcopy(st,pc->et,kn+kk,double);

     /* Copying the coeffisient vector from the old curve to the new curve.*/

     if (ider < pc->ik) memcopy(scoef,pc->ecoef,pc->in*kdim,double);

     /* Here we are computing a new coeffecient vector for each round. */

     if (ider < pc->ik)
       for (kj=1; kj<=ider; kj++)
         {
	   kk =pc->ik - kj;      /* The new order of the curve. */

	   s1 = scoef + (pc->in-1+kj)*kdim; /* The last new vertice. */

	   /* Here I just refere to the referenses.
	      The new vertices are computig from back to front,
	      and the "last" kdim vertices is beeing computed outside
	      the main loop because we did not have scoef[-1],
	      instead we use zeroes. */

	   for (s3=st+pc->in-1+kj; st<s3; s3--,s1-=2*kdim)
             {
	       tdel = s3[kk] - *s3;

	       if (DNEQUAL(tdel,DZERO))
	         for (s2=s1+kdim; s1<s2; s1++)
		   *s1=(*s1-s1[-kdim])*kk/tdel;
	       else
	         for (s2=s1+kdim; s1<s2; s1++) *s1 = DZERO;
             }

	   tdel = s3[kk] - *s3;

	   if (DNEQUAL(tdel,DZERO))
	     for (s2=s1+kdim; s1<s2; s1++) *s1 = *s1*kk/tdel;
	   else
             for (s2=s1+kdim; s1<s2; s1++) *s1 = DZERO;

         }

     /* Allocating new curve-object.*/

     if ((q1=newCurve(kn-2,kk,&st[1],&scoef[pc->idim],1,pc->idim,1))
	 == SISL_NULL) goto err101;
     freearray(st);
     freearray(scoef);
  }


  /* Remove unnessesary internal knots and vertises. */

  s1705(q1,&kstat);
  if (kstat < 0) goto error;

  /* Set periodicity flag if cyclic curve.  */

  test_cyclic_knots(q1->et,q1->in,q1->ik,&kstat);
  if (kstat<0) goto error;
  if (kstat == 2) q1->cuopen = SISL_CRV_PERIODIC;

  /* Updating output. */

  *rcnew = q1;
  goto out;

  /* Error. Error in lower level function. */

 error:
  *jstat = kstat;
  goto outfree;

  /* Error. No curve to subdevice.  */

 err150:
  *jstat = -150;
  s6err("s1720",*jstat,kpos);
  goto out;

  /* Error. Ileagal number of derevatives.  */

 err156:
  *jstat = -156;
  s6err("s1720",*jstat,kpos);
  goto out;

  /* Error. Allocation error, not enough memory.  */

 err101:
  *jstat = -101;
  s6err("s1720",*jstat,kpos);
  goto outfree;

 outfree:
  if(q1) freeCurve(q1);
  else
    {
      if (st) freearray(st);
      if (scoef) freearray(scoef);
    }

  /* Free local used memory. */

 out: return;
}
