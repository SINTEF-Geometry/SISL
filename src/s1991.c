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
 * $Id: s1991.c,v 1.4 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1991

#include "sislP.h"
/*
#if defined(SISLNEEDPROTOTYPES)
static void s1991_s9smooth(double [],int,int,double,double [],int *);
#else
static void s1991_s9smooth();
#endif
*/

#if defined(SISLNEEDPROTOTYPES)
void
     s1991(SISLCurve *pc,double aepsge,int *jstat)
#else
void s1991(pc,aepsge,jstat)
     SISLCurve  *pc;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To make the orientation surface on the unit sphere to
*	       a b-spline curve, the surface is representated with
*	       a surrounding cone piced from the unit sphere.
*
*
*
* INPUT      : pc     - The orginal B-spline curve.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We are making a cone surrounding the orientating surface
*	       on the unit sphere. The cone is representated with senter
*	       coordinates and an angle. The orientation is computed
*	       from aproximation of the tangent to the curve.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* CORRECTED BY: Ulf J. Krystad, SI, 91-07
*               Problems in shevalc when clustering ceoeff's in s9smooth
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 09/09-1994. Commented out
*              call to s1991_s9smooth() and added call to memcopy() instead,
*              according to advice from Vibeke Skytt.
*********************************************************************
*/
{
  int kpos = 0;     /* Position of the error.                          */
  int kfirst = 1;   /* Flag to mark if the first tangent is treating.  */
  int kn;           /* Number of vertices of curve.                    */
  int kdim;	    /* Dimension of the space in which the objects lie.*/
  int kin;          /* The index to the vertice to treat.              */
  int k1,k2;        /* Control variables in loop.                      */
  double *t=SISL_NULL;   /* Tangent at each coeficient.                     */
  double tlen;      /* The length of a vector.                         */
  double tang;	    /* An angle between two vectors.		       */
  double t1,t2;     /* Help variables.				       */
  double *scoef;    /* Pointer to coefficients.                        */



  /* Test if the surfaces already have been treated.  */

  if (pc->pdir != SISL_NULL) goto out;


  /* Initialate dimentions. */

  kdim = pc -> idim;
  kn = pc -> in;


  /* Make a new direction cone. */

  if ((pc->pdir = newdir(kdim))==SISL_NULL) goto err101;

  /* UJK, Set default values in pdir. */
  pc->pdir->aang = DZERO;
  pc->pdir->igtpi = 0;
  pc->pdir->ecoef[0] = (double) 1.0;

  for (k2 = 1;k2<kdim;k2++)
    pc->pdir->ecoef[k2] = DZERO;


  /* Allocate local used array. */

  if ((t = newarray(kdim,double)) == SISL_NULL) goto err101;

  /* Allocate scratch for smoothed coefficients.  */

  if ((pc->pdir->esmooth = newarray(kn*kdim,DOUBLE)) == SISL_NULL) goto err101;
  scoef = pc->pdir->esmooth;

  /* Compute coefficients of smoothed curve.  */

   /* s1991_s9smooth(pc->ecoef,kn,kdim,aepsge,scoef,&kstat);
      if (kstat < 0) goto error; */
   /* (VSK 02-1994: no point in smoothing) */
   memcopy(scoef, pc->ecoef, kn*kdim, DOUBLE);

  /* Here we are treating each patch in the control polygon separately.*/

  for (k2=0,kin=0; kin < kn-1; kin++)
    {

      /* Here we make an aproximative tangents to the curve
	 using the control polygon. The tangents is also normalized
	 by deviding with its own length. */

      for (tlen=DZERO,k1=0; k1 < kdim; k1++,k2++)
	{
	  t[k1] = scoef[k2+kdim] - scoef[k2];
	  tlen += t[k1]*t[k1];
	}

      tlen = sqrt(tlen);

      if (tlen > aepsge)
	for (k1=0; k1 < kdim; k1++) t[k1] /= tlen;
      else
	{
	  /* UJK, whats wrong with colapsed polygons when computing directions? */
	  continue;

	  /* Vi have to be aware of colapsed polygon. */
	  /* pc->pdir->igtpi = 1;
	     goto out;             */

	}


      /* We are treating the first tangent. */

      if (kfirst)
	{

	  /* Computing the center coordinates of the cone.*/

	  for (k1=0; k1 < kdim; k1++)
	    pc->pdir->ecoef[k1]= t[k1];

	  /* Computing the angle of the cone. */

	  pc->pdir->aang = DZERO;

	  kfirst = 0;   /* The first tangent have been treated.*/
	}
      else
	{

	  /* Computing the angle beetween the senter of the cone
	     and the tangent. */

	  for (tang=DZERO,k1=0;k1<kdim;k1++)
	    tang += pc->pdir->ecoef[k1]*t[k1];

	  if (tang >= DZERO) tang = min((double)1.0,tang);
	  else               tang = max((double)-1.0,tang);

	  tang = acos(tang);

	  if (tang + pc->pdir->aang >= PI)
	    {
	      /* The angle is to great, give a meesage
		 to subdivied and exit this function. */

	      pc->pdir->igtpi = 1;
	      goto out;
	    }
	  else if (tang > pc->pdir->aang)
	    {
	      /* The tangent is not inside the cone, and we
		 have to compute a new cone. */

	      /* Computing the center coordinates.*/

	      t1 = (tang - pc->pdir->aang)/((double)2*tang);
	      t2 = (double)1 - t1;

	      for (tlen=DZERO,k1=0; k1<kdim; k1++)
		{
		  pc->pdir->ecoef[k1] =
		    pc->pdir->ecoef[k1]*t2 + t[k1]*t1;
		  tlen += pc->pdir->ecoef[k1]*
		    pc->pdir->ecoef[k1];
		}
	      tlen = sqrt(tlen);

	      if (tlen > DZERO)
		for (k1=0; k1 < kdim; k1++)
		  pc->pdir->ecoef[k1] /= tlen;
	      else
		{
		  /* Vi have to be aware of colapsed polyg.*/

		  pc->pdir->igtpi = 1;
		  goto out;
		}


	      /* Computing the angle of the cone. */

	      pc->pdir->aang = (tang + pc->pdir->aang)/
		(double)2;
	    }
	}
    }



  if (pc->pdir->aang >= SIMPLECASE)
    {
      /* The angle is to great, give a message
	 to subdivied and exit this function. */

      pc->pdir->igtpi = 3;
      goto out;
    }


  *jstat = 0;
  goto out;


  /* Error in space allacation.  */

 err101: *jstat = -101;
  s6err("s1991",*jstat,kpos);
  goto out;

 out:    if (t != SISL_NULL) freearray(t);

}

#if 0
#if defined(SISLNEEDPROTOTYPES)
static void
  s1991_s9smooth(double ecoef1[],int in,int idim,double aepsge,
		 double ecoef2[],int *jstat)
#else
static void s1991_s9smooth(ecoef1,in,idim,aepsge,ecoef2,jstat)
   double ecoef1[];
   int    in;
   int    idim;
   double aepsge;
   double ecoef2[];
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Perform noise filthering on the control polygon on a
*              B-spline curve with special emphasis to the ends of
*              the curve.
*
*
*
* INPUT      : ecoef1 - Original coefficients of curve.
*              in     - Number of coefficients.
*              idim   - Dimension of geometry space.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : ecoef2 - New coefficients after smoothing.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Start from both ends of the curve and traverse towards
*              the middle. Coefficients with distance less than the
*              tolerance from a line between nearby coefficients, are
*              projected down to this line.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s6dline  -  Distance between point and line.
*              s6dist   -  Distance between points.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-02.
* CORRECTED BY: Ulf J. Krystad, SI, 91-07
*               Problems in shevalc when clustering ceoeff's
*********************************************************************
*/
{
   int kstat = 0;      /* Local status variable.            */
   int kn2 = in/2;     /* Half the number of coefficients.  */
   int ki;             /* Loop control.                     */
   double *sdiff1=SISL_NULL;/* Diff vector                       */
   double *sdiff2=SISL_NULL;/* Diff vector                       */
   double alfa,dnum;   /* Factor and denominator in expr.   */
   double tdist;       /* Distance between point and line.  */
   double *s1,*s2,*s3; /* Pointers into coefficient array.  */
   double *st1,*st2;   /* Stop pointers in loop.            */

   /* Alloc scratch for locals */
   if ((sdiff1 = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   if ((sdiff2 = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;


   /* Copy coefficient array to output array.  */
   memcopy(ecoef2,ecoef1,in*idim,DOUBLE);

   /* Traverse and smooth first half of the coefficient array.  */

   for (s1=ecoef2, st1=s1+(kn2-1)*idim; s1<st1; s1=s2)
   {
      for (s2=s1+2*idim, st2=st1+idim; s2<=st2; s2+=idim)
      {
	 if (s6dist(s1,s2,idim) < aepsge) continue;

	 for (s3=s1+idim; s3<s2; s3+=idim)
	 {
	    /* Find distance between the point s3 and the line
	       segment between s1 and s2.   */

	    tdist = s6dline(s1,s2,s3,idim,&kstat);
	    if (kstat < 0) goto error;

	    /* Test if the point is close to a point within the
	       line segment, and the distance is less than
	       the tolerance.        */

	    if (kstat || tdist >= aepsge) break;  /* No smoothing
						     possible.  */
	 }
	 if (s3 < s2) break; /* No smoothing between s1 and s2. */
      }

      /* Project all coefficients between s1 and s2 down to the
	 line segment between s1 and s2.  */

      s2 -= idim;
      /*UJK Problems in shevalc when clustering ceoeff's,
	 let's really project! */
      /*for (s3=s1+idim; s3<s2; s3+=idim)
	 memcopy(s3,s1,idim,DOUBLE); */

      s6diff(s2,s1,idim,sdiff1);
      dnum = s6scpr(sdiff1,sdiff1,idim);


      for (s3=s1+idim; s3<s2; s3+=idim)
      {
	  if (dnum > DZERO)
	     {
		s6diff(s2,s3,idim,sdiff2);
		alfa = s6scpr(sdiff2,sdiff1,idim)/dnum;
		for (ki=0;ki<idim;ki++) s3[ki]=
		   alfa*s1[ki] + ((double)1.0 - alfa)*s2[ki];
	     }
      }

   }

   /* Traverse and smooth second half of the coefficient array.  */

   for (s1=ecoef2+(in-1)*idim, st1=s1-(kn2-1)*idim; s1>st1; s1=s2)
   {
      for (s2=s1-2*idim, st2=st1-idim; s2>=st2; s2-=idim)
      {
	 if (s6dist(s1,s2,idim) < aepsge) continue;

	 for (s3=s1-idim; s3>s2; s3-=idim)
	 {
	    /* Find distance between the point s3 and the line
	       segment between s1 and s2.   */

	    tdist = s6dline(s1,s2,s3,idim,&kstat);
	    if (kstat < 0) goto error;

	    /* Test if the point is close to a point within the
	       line segment, and the distance is less than
	       the tolerance.        */

	    if (kstat || tdist >= aepsge) break;  /* No smoothing
						     possible.  */
	 }
	 if (s3 > s2) break; /* No smoothing between s1 and s2. */
      }

      /* Project all coefficients between s1 and s2 down to the
	 line segment between s1 and s2.  */

      s2 += idim;
      /*UJK Problems in shevalc when clustering ceoeff's,
	 let's really project! */
      /*for (s3=s1-idim; s3>s2; s3-=idim)
	 memcopy(s3,s1,idim,DOUBLE); */

      s6diff(s2,s1,idim,sdiff1);
      dnum = s6scpr(sdiff1,sdiff1,idim);

      for (s3=s1-idim; s3>s2; s3-=idim)
      {
	  if (dnum > DZERO)
	     {
		s6diff(s2,s3,idim,sdiff2);
		alfa = s6scpr(sdiff2,sdiff1,idim)/dnum;
		for (ki=0;ki<idim;ki++) s3[ki]=
		   alfa*s1[ki] + ((double)1.0 - alfa)*s2[ki];
	     }
      }
   }

   /* Smoothing performed.  */

   *jstat = 0;
   goto out;

  /* Error in space allacation.  */

   err101: *jstat = -101;
   goto out;

   /* Error in lower level routine.  */

   error : *jstat = kstat;
   goto out;

   out :
      if (sdiff1) freearray(sdiff1);
      if (sdiff2) freearray(sdiff2);
      return;
}
#endif /* if 0 */
