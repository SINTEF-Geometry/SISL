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
 * $Id: s1711.c,v 1.4 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1711

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1711(SISLSurf *ps,int ipar,double apar,SISLSurf **rsnew1,SISLSurf **rsnew2,int *jstat)
#else
void s1711(ps,ipar,apar,rsnew1,rsnew2,jstat)
     SISLSurf   *ps;
     int    ipar;
     double apar;
     SISLSurf   **rsnew1;
     SISLSurf   **rsnew2;
     int    *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Subdivide a B-spline surface at a given parameter-value
*	       and a given direction.
*
*
*
* INPUT      : ps	- Surface to subdivide.
*	       ipar	- 1 means first direction, 2 means second direction.
*   	       apar	- Parameter-value at which to subdivide.
*
*
*
* OUTPUT     : rsnew1	- First part of the subdivided surface.
*              rsnew2	- Second part of the subdivided surface.
*              jstat	- status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Inserting knots at the subdividing-point until
*	       we have a ktuple-knot. Then we may separate the
*	       surface into two parts.
*	       The problem is to treat a three dimensional array that
*	       is locatet in a one dimensional array (coeffisients).
*	       To solve this problem we use constants to marsch in different
*	       direction in the array, and to mark end of lines or columns.
*
*
* REFERENCES :
*
*-
* CALLS      : newSurf  - Allocate space for a new curve-object.
*	       freeSurf - Free space occupied by given curve-object.
*	       S1700.C  - Making the knot-inserten-transformation matrix.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* MODIFIED BY : Mike Floater, SI, 90-12. Subdivide rational surfaces.
* MODIFIED BY : Arne Laksaa, SI, 92-09. Move apar to closest knot
*		if it is close to a knot. Using D(N)EQUAL().
* MODIFIED BY : Paal Fugelli, SINTEF, Oslo 15/07-1994. Changed from icopy==0
*               to icopy==2 in creation of output surface - memory leak.
*
**********************************************************************/
{
  int kstat;		/* Local status variable.		*/
  int kpos=0;		/* Position of error.			*/
  int kmy;		/* An index to the knot-vector.		*/
  int kv,kv1;		/* Number of knots we have to insert.	*/
  int kpl,kfi,kla;	/* To posisjon elements in trans.-matrix.*/
  int kk,kksec;		/* Order of the input surface.		*/
  int kn,knsec;		/* Number of the vertices in input curves.*/
  int kdim=ps->idim;	/* Dimensjon of the space in whice surface
			   lies.				*/
  int kind=ps->ikind;	/* Type of surface ps is.               */
  int kn1,kn2;		/* Number of vertices in the new surfaces.*/
  int knum;		/* Number of knots less and equal than
			   the intersection point.		*/
  int ki,ki1,ki2;	/* Control variable in loop.		*/
  int kj,kj1,kj2;	/* Control variable in loop.		*/
  int k1m,k2m,k3m,k4m;	/* Variables to mark directons in array.*/
  int newkind=1;	/* Type of surface subsurfaces are.     */
  double *s1,*s2,*s3,*s4;/* Pointers used in loop.		*/
  double *st,*stsec;	/* The old knot-vectors.		*/
  double *st1=SISL_NULL;	/* The first first new knot-vector.	*/
  double *st1sec=SISL_NULL;	/* The first second new knot-vector.	*/
  double *st2=SISL_NULL;	/* The second first new knot-vector.	*/
  double *st2sec=SISL_NULL;	/* The second second new knot-vector.	*/
  double *salfa=SISL_NULL;	/* A line of the trans.-matrix.		*/
  double *scoef;	/* Pointer to vertices.   		*/
  double *scoef1=SISL_NULL;	/* The first new vertice.		*/
  double *scoef2=SISL_NULL;	/* The second new vertice.		*/
  SISLSurf *q1=SISL_NULL;	/* Pointer to new surface-object.	*/
  SISLSurf *q2=SISL_NULL;	/* Pointer to new surface-object.	*/
  double salfa_local[5];/* Local help array.			*/

 /* if ps is rational, do subdivision in homogeneous coordinates */
 /* just need to set up correct dim and kind for the new surfaces at end of routine */
  if(kind == 2 || kind == 4)
  {
       scoef = ps->rcoef;
       kdim++;
       newkind++;
  }
  else
  {
       scoef = ps->ecoef;
  }

  /* Check that we have a surface to subdivide. */

  if (!ps) goto err150;

  /* Making constants and ponters to mark direction.  */

  if (ipar==1)
    {
      /* If ipar is 1 we have to split the "three" dimentional
	 coeffisient matrix along a column. In this case k4m is
	 the distance beetween each element in the clumn.
	 For each element in the column we have to treat a part
	 of a line, to march along the line we use k1m.*/

      st = ps->et1;
      stsec = ps->et2;
      kn = ps->in1;
      knsec = ps->in2;
      kk = ps->ik1;
      kksec = ps->ik2;
      k1m = kdim;
      k4m = kdim*kn;
    }
  else
    {
      /* If ipar is 2 we have to split the "three" dimentional
	 coeffisient matrix along a line. In this case k4m is
	 the distance beetween each element in the line.
	 For each element in the line we have to treat a part
	 of a column, to march along the column we use k1m.*/

      st = ps->et2;
      stsec = ps->et1;
      kn = ps->in2;
      knsec = ps->in1;
      kk = ps->ik2;
      kksec = ps->ik1;
      k1m = kdim*knsec;
      k4m = kdim;
    }

  /* Check that the intersection point is an interior point. */

  if ((apar < *st  && DNEQUAL(apar, *st)) ||
      (apar > st[kn+kk-1] && DNEQUAL(apar, st[kn+kk-1])))
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

  s1 = st;
  kv = kk;	/* The maximum number of knots we have to insert. */

  if ((apar > s1[0] && DNEQUAL(apar, s1[0])) &&
      (apar < s1[kn+kk-1] && DNEQUAL(apar, s1[kn+kk-1])))
  {
     /* Using binear search*/
     kj1=0;
     kj2=kk+kn-1;
     knum = (kj1+kj2)/2;
     while (knum != kj1)
     {
	if ((s1[knum] < apar ) && DNEQUAL(s1[knum], apar))
	   kj1=knum; else kj2=knum;
	knum = (kj1+kj2)/2;
     }
     knum++;           /* The smaller knots.*/

     while (DEQUAL(s1[knum], apar))
     {
	apar = s1[knum];
	knum++;
	kv--;
     }
     /* The knots thats like the */
     /*     intersection point.  */
  }
  else if (DEQUAL(apar,s1[0]))
  {
     apar = s1[0];
     knum = 0;
     while (s1[knum] == apar)
	/* The knots thats like the intersection point. */
	knum++;
  }
  else if (DEQUAL(apar,s1[kn+kk-1]))
  {
     apar = s1[kn+kk-1];
     knum = kn+kk-1;
     while (s1[knum-1] == apar)
	/* The knots thats like the intersection point. */
	knum--;
  }
  /* Find the number of vertices in the two new curves. */

  kn1 = knum + kv - kk;
  kn2 = kn + kk - knum;

  /* Allocating the new arrays to the two new curves. */

  if ((st1=newarray(kn1+kk,double))==SISL_NULL) goto err101;
  if ((st1sec=newarray(knsec+kksec,double))==SISL_NULL) goto err101;
  if ((st2=newarray(kn2+kk,double))==SISL_NULL) goto err101;
  if ((st2sec=newarray(knsec+kksec,double))==SISL_NULL) goto err101;
  if ((scoef1=newarray(kn1*kdim*knsec,double))==SISL_NULL) goto err101;
  if ((scoef2=newarray(kn2*kdim*knsec,double))==SISL_NULL) goto err101;

  /* Copying the knotvectors from the old curve to the new curves */

  memcopy(st1,st,kn1,double);
  memcopy(st2+kk,st+knum,kn2,double);
  memcopy(st1sec,stsec,knsec+kksec,double);
  memcopy(st2sec,stsec,knsec+kksec,double);

  /* Updating the knotvectors by inserting the new k-touple knot */

  for(s2=st1+kn1,s3=st2,s4=s3+kk;s3<s4;s2++,s3++) *s2 = *s3 = apar;

  /* Copying the coefisientvectors to the new curves.*/

  if (ipar == 1)
    for (ki=0; ki<knsec; ki++)
      {
	memcopy(scoef1+ki*kdim*kn1,scoef+ki*kdim*kn,
		kdim*kn1,double);
	memcopy(scoef2+ki*kdim*kn2,scoef+kdim*(ki*kn+knum-kk),
		kdim*kn2,double);
      }
  else
    {
      memcopy(scoef1,scoef,kdim*kn1*knsec,double);
      memcopy(scoef2,scoef+kdim*(knum-kk)*knsec,
	      kdim*kn2*knsec,double);
    }

  /* Updating the coefisientvectors to the new surfaces.*/

  /* Updating the first surface. */

  /* If we imagine that the matrix is turned in such a way that we are
     splitting it along a column, then for each element in the column
     we have to treat a par of a line, to march along the line
     in the first new matrix we use k1m, And we use k3m as a mark
     at the end of the column in this new matrix.*/

  if(ipar==1)
    {
      k2m=kdim*kn1;
      k3m=kdim*kn1*knsec;
    }
  else
    {
      k2m=kdim;
      k3m=kdim*knsec;
    }
  knum -= kk - 1;
  for (ki=max(0,knum),kv1=max(0,-knum),s1=scoef1+ki*k1m;ki<kn1;ki++,s1+=k1m)
    {
      /* Initialising:
	 knum = knum-kk+1, Index of the first vertice to change.
	 ki = knum, 	  Index of the vertices we are going to
	 change. Starting with knum, but if
	 knum is negativ we start at zero.
	 kv1 = 0,	  Number if new knots between index ki
	 and ki+kk. We are starting one below
	 becase we are counting up before using
	 it. If knum is negativ we are not
	 starting at zero but at -knum.
	 s1=scoef1+ki*k1m  Pointer at the first vertice to
	 change. */

      /* Using the Oslo-algorithm to make a transformation-vector
	 from the old vertices to one new vertice. */

      kmy=ki;
      s1700(kmy,kk,kn,++kv1,&kpl,&kfi,&kla,st,apar,salfa,&kstat);
      if (kstat) goto err153;

      /* Compute the knsec*kdim vertices with the "same index". */

      for (s2=s1,s3=s2+k3m,ki2=0; s2<s3; s2+=k2m,ki2+=k4m)
	for (kj=0,s4=s2; kj<kdim; kj++,s4++)
	  for (*s4=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla;kj1++,kj2++)
	    *s4 += salfa[kj2] * scoef[k1m*kj1+ki2+kj];
    }

  /* And the second surface. */

  /* If we imagine that the matrix is turned in such a way that we are
     splitting it along a column, then for each element in the column
     we have to treat a par of a line, to march along the line
     in the second new matrix we use k1m, And we use k3m as a mark
     at the end of the column in this new matrix.*/

  if(ipar==1)
    {
      k2m=kdim*kn2;
      k3m=kdim*kn2*knsec;
    }
  else
    {
      k2m=kdim;
      k3m=kdim*knsec;
    }

  for (ki1=min(kn1+kv-1,kn+kv),s1=scoef2; ki<ki1; ki++,s1+=k1m)
    {
      /* Initialising:
	 ki1 = kn1+kv-1,	  the index of the vertice next to the
	 last vertice we have to change.
	 If we do not have so many vertices,
	 we have to use the index next to the
	 last vertice we have, kn+kv.
	 s1=scoef2	  Pointer at the first vertice to
	 change. */


      /* Using the Oslo-algorithm to make a transformation-vector
	 from the old vertices to one new vertice. */

      s1700(kmy,kk,kn,kv1--,&kpl,&kfi,&kla,st,apar,salfa,&kstat);
      if (kstat) goto err153;


      /* Compute the knsec*kdim vertices with the "same index". */

      for (s2=s1,s3=s2+k3m,ki2=0; s2<s3; s2+=k2m,ki2+=k4m)
	for (kj=0,s4=s2; kj<kdim; kj++,s4++)
	  for (*s4=0,kj1=kfi,kj2=kfi+kpl; kj1<=kla;kj1++,kj2++)
	    *s4 += salfa[kj2] * scoef[k1m*kj1+ki2+kj];
    }


  /* Allocating new surface-objects.*/
 /* use ps->idim rather than kdim in case ps is rational  */


  if (ipar==1)
  {
    if ((q1=newSurf(kn1,knsec,kk,kksec,st1,st1sec,     /* PFU 15/07-94 */
                    scoef1,newkind,ps->idim,2)) == SISL_NULL) goto err101;
    if ((q2=newSurf(kn2,knsec,kk,kksec,st2,st2sec,     /* PFU 15/07-94 */
                    scoef2,newkind,ps->idim,2)) == SISL_NULL) goto err101;
  }
  else
  {
    if ((q1=newSurf(knsec,kn1,kksec,kk,st1sec,st1,     /* PFU 15/07-94 */
                    scoef1,newkind,ps->idim,2)) == SISL_NULL) goto err101;
    if ((q2=newSurf(knsec,kn2,kksec,kk,st2sec,st2,     /* PFU 15/07-94 */
                    scoef2,newkind,ps->idim,2)) == SISL_NULL) goto err101;
  }


  /* Updating output. */

  *rsnew1 = q1;
  *rsnew2 = q2;
  *jstat = 0;
  goto out;


  /* Error. Error in lower level function. */

 err153: *jstat = kstat;
  goto outfree;


  /* Error. No surface to subdevice.  */

 err150: *jstat = -150;
  s6err("s1711",*jstat,kpos);
  goto out;


  /* Error. The intersection-point is outside the surface.  */

 err158: *jstat = -158;
  s6err("s1711",*jstat,kpos);
  goto out;


  /* Error. Allocation error, not enough memory.  */

 err101: *jstat = -101;
  s6err("s1711",*jstat,kpos);
  goto outfree;


outfree:
   if(q1) freeSurf(q1);
   if(q2) freeSurf(q2);

   /* Free local used memory. */

out:
   if(!q1)
   {
      if (st1) freearray(st1);
      if (st1sec) freearray(st1sec);
      if (scoef1) freearray(scoef1);
   }

   if(!q2)
   {
      if (st2) freearray(st2);
      if (st2sec) freearray(st2sec);
      if (scoef2) freearray(scoef2);
   }

   if (kk > 5 && salfa)
      freearray (salfa);
   return;
}
