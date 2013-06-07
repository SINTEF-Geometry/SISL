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
 * $Id: sh6idrcros.c,v 1.4 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6IDRMCROSS

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
      sh6idrmcross(SISLObject *po1, SISLObject *po2, SISLIntdat **pintdat,
		   SISLIntpt *vcross[], int incross, SISLIntpt *vpt[],
		   int inpt, int *jstat)
#else
void sh6idrmcross(po1, po2, pintdat, vcross, incross, vpt, inpt,jstat)
   SISLObject *po1;
   SISLObject *po2;
   SISLIntdat **pintdat;
   SISLIntpt  *vcross[];
   int        incross;
   SISLIntpt  *vpt[];
   int        inpt;
   int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an array of 4 cross intersection. These are
*              found by sh6idfcross. Check if the
*              intersections may be divided into two pairs where
*              the second pair consists of the cross intersections
*              of the first pair. In that case remove the cross
*              intersections.
*
*
* INPUT      : po1      - First object in the intersection.
*              po2      - Second object in the intersection.
*              pintdat  - Intersection data.
*              vcross   - Array of cross intersection point.
*              incross  - Number of points in vcross. incross = 4.
*              vpt      - Array of pointers into pintdat.
*              inpt     - Number of points in vpt.
*
* OUTPUT     : pintdat  - Updated intersection data.
*              vpt      - Updated array.
*              jstat    - Status
*                         jstat = 0   => No intersection points removed.
*                         jstat = 1   => Cross intersections removed.
*                         jstat < 0   => Error.
*
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* CALLS      : s1221, s1424, s6degnorm, s6scpr, s6length, sh6idkpt.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 12.92.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994. Init of
*              kleft1 and kleft2.
*
*********************************************************************
*/
{
   int kstat;          /* Status variable.                          */
   int ki,kj,kl;       /* Counters.                                 */
   int kdim;           /* Dimension of geometry space.              */
   int kleft1 = 0;     /* Parameters used in evaluator.             */
   int kleft2 = 0;
   int kdir1,kdir2;    /* Parameter directions.                     */
   int kpar;           /* Number of parameter directions.           */
   int k1par = po1->iobj; /* Number of par. dir. in first object.   */
   int kmin;           /* Number of minimum parameter point.        */
   double tmin;        /* Length of minimum parameter point as vector. */
   double thelp;       /* Help parameter. Length of vector.         */
   double sder1[27];   /* Position, derivative etc. of object 1.    */
   double sder2[27];   /* Position, derivative etc. of object 2.    */
   double stang1[3];   /* Tangent in first parameter dir., deg. surf. */
   double stang2[3];   /* Tangent in second parameter dir., deg. surf. */
   double snorm[3];    /* Normal of degenerated surface.            */

   *jstat = 0;

   /* Test input.  */

   if (incross != 4) goto err138;

   if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
   {
      *jstat = 0;
      goto out;
   }

   if (po1->iobj == SISLSURFACE)
   {
      /* Check if the intersection points have got one parameter in
	 common.      */

      for (kj=0; kj<2; kj++)
      {
	 for (ki=1; ki<incross; ki++)
	    if (DNEQUAL(vcross[ki]->epar[kj],vcross[0]->epar[kj])) break;
	 if (ki == incross) break;  /* Common parameter direction.  */
      }

      if (kj == 2)
      {
	 /* No common parameter direction, i.e. no cross intersection
	    to remove.  */

	 *jstat = 0;
	 goto out;
      }

      /* Set the parameter direction that is not constant.  */

      kdir1 = 1 - kj;
   }

   if (po2->iobj == SISLSURFACE)
   {
      /* Check if the intersection points have got one parameter in
	 common.      */

      for (kj=po1->iobj, kpar=vcross[0]->ipar; kj<kpar; kj++)
      {
	 for (ki=1; ki<incross; ki++)
	    if (DNEQUAL(vcross[ki]->epar[kj],vcross[0]->epar[kj])) break;
	 if (ki == incross) break;  /* Common parameter direction.  */
      }

      if (kj == kpar)
      {
	 /* No common parameter direction, i.e. no cross intersection
	    to remove.  */

	 *jstat = 0;
	 goto out;
      }

      /* Set the parameter direction that is not constant.  */

      kdir2 = kpar - 1 - kj;
   }

   /* Find the minimum parameter point in which to evaluate. */

   kmin = 0;
   tmin = s6length(vcross[0]->epar,vcross[0]->ipar,&kstat);

   for (kj=1; kj<incross; kj++)
   {
      thelp = s6length(vcross[kj]->epar,vcross[kj]->ipar,&kstat);
      if (thelp < tmin)
      {
	 tmin = thelp;
	 kmin = kj;
      }
   }

   /* Compute derivatives.  */

   if (po1->iobj == SISLCURVE)
   {
      kdir1 = 0;
      kdim = po1->c1->idim;
      s1221(po1->c1,1,vcross[kmin]->epar[kdir1],&kleft1,sder1,&kstat);
      if (kstat < 0) goto error;
   }
   else if (po1->iobj == SISLSURFACE)
   {
      kdim = po1->s1->idim;
      s1424(po1->s1,2,2,vcross[kmin]->epar,&kleft1,&kleft2,sder1,
	    &kstat);
      if (kstat < 0) goto error;
      s6crss(sder1+kdim, sder1+2*kdim, snorm);

      /* Check if the surface is degenerated in the wanted parameter
	 direction.   */

      if (s6length(sder1+(1+kdir1)*kdim,kdim,&kstat) <= REL_COMP_RES)
      {
	 /* Compute partial derivatives as a limit.  */

	 s6degnorm(po1->s1,2,vcross[kmin]->epar,sder1,stang1,stang2,
		   snorm,&kstat);
	 if (kstat < 0) goto error;

	 memcopy(sder1+kdim,stang1,kdim,DOUBLE);
	 memcopy(sder1+2*kdim,stang2,kdim,DOUBLE);
      }
   }

   /* Compute derivatives.  */

   if (po2->iobj == SISLCURVE)
   {
      kdir2 = 0;
      s1221(po2->c1,1,vcross[kmin]->epar[k1par+kdir2],&kleft1,sder2,&kstat);
      if (kstat < 0) goto error;
   }
   else if (po2->iobj == SISLSURFACE)
   {
      s1424(po2->s1,2,2,vcross[kmin]->epar+k1par,&kleft1,&kleft2,sder2,
	    &kstat);
      if (kstat < 0) goto error;
      s6crss(sder2+kdim, sder2+2*kdim, snorm);

      /* Check if the surface is degenerated in the wanted parameter
	 direction.   */

      if (s6length(sder2+(1+kdir2)*kdim,kdim,&kstat) <= REL_COMP_RES)
      {
	 /* Compute partial derivatives as a limit.  */

	 s6degnorm(po2->s1,2,vcross[kmin]->epar+k1par,sder2,stang1,stang2,
		   snorm,&kstat);
	 if (kstat < 0) goto error;

	 memcopy(sder2+kdim,stang1,kdim,DOUBLE);
	 memcopy(sder2+2*kdim,stang2,kdim,DOUBLE);
      }
   }

   if (s6ang(sder1+(kdir1+1)*kdim,sder2+(kdir2+1)*kdim,kdim) > ANGULAR_TOLERANCE
       && !(s6length(sder1+(kdir1+1)*kdim,kdim,&kstat) < REL_COMP_RES &&
	    s6length(sder2+(kdir2+1)*kdim,kdim,&kstat) < REL_COMP_RES))
   {
      *jstat = 0;
      goto out;
   }

   /* Check if the parameter directions of the objects are the same.  */

   if (s6scpr(sder1+(kdir1+1)*kdim,sder2+(kdir2+1)*kdim,kdim) >= 0)
   {
      /* Remove the pair of intersection point that do not have the
	 same "parameter direction" in both objects.  */

      for (ki=0; ki<incross; ki++)
	{
	  for (kj=1; kj<incross; kj++)
	    {
	      if (DNEQUAL(vcross[ki]->epar[kdir1],vcross[kj]->epar[kdir1]) &&
		  DNEQUAL(vcross[ki]->epar[k1par+kdir2],
			  vcross[kj]->epar[k1par+kdir2]))
		{
		  /* A pair is found. Check if the pair should be removed. */

		  if ((vcross[ki]->epar[kdir1] - vcross[kj]->epar[kdir1]) *
		      (vcross[ki]->epar[k1par+kdir2] -
		       vcross[kj]->epar[k1par+kdir2]) < 0)
		    {
		      /* Remove the points. First make sure that vpt will
			 not point to killed points.  */

		      for (kl=0; kl<inpt; kl++)
			if (vpt[kl] == vcross[ki] || vpt[kl] == vcross[kj])
			  vpt[kl] = SISL_NULL;

		      sh6idkpt(pintdat,&vcross[ki],1,&kstat);
		      if (kstat < 0) goto error;

		      sh6idkpt(pintdat,&vcross[kj],1,&kstat);
		      if (kstat < 0) goto error;

		      *jstat = 1;
		      break;
		    }
		}
	    }
	  if (kj < incross)
	    break;           /* The cross intersection is removed */
	}

      if (*jstat == 1) goto out;  /* Points removed.  */
   }
   else
   {
      /* Remove the pair of intersection point that have the
	 same "parameter direction" in both objects.  */

      for (ki=0; ki<incross; ki++)
	{
	  for (kj=1; kj<incross; kj++)
	    {
	      if (DNEQUAL(vcross[ki]->epar[kdir1],vcross[kj]->epar[kdir1]) &&
		  DNEQUAL(vcross[ki]->epar[k1par+kdir2],
			  vcross[kj]->epar[k1par+kdir2]))
		{
		  /* A pair is found. Check if the pair should be removed. */

		  if ((vcross[ki]->epar[kdir1] - vcross[kj]->epar[kdir1]) *
		      (vcross[ki]->epar[k1par+kdir2] -
		       vcross[kj]->epar[k1par+kdir2]) > 0)
		    {
		      /* Remove the points. First make sure that vpt will
			 not point to killed points. */

		      for (kl=0; kl<inpt; kl++)
			if (vpt[kl] == vcross[ki] || vpt[kl] == vcross[kj])
			  vpt[kl] = SISL_NULL;

		      sh6idkpt(pintdat,&vcross[ki],1,&kstat);
		      if (kstat < 0) goto error;

		      sh6idkpt(pintdat,&vcross[kj],1,&kstat);
		      if (kstat < 0) goto error;

		      *jstat = 1;
		      break;
		    }
		}
	    }
	  if (kj < incross)
	    break;           /* The cross intersection is removed */
	 }

      if (*jstat == 1) goto out;  /* Points removed.  */
   }

      /* No points are removed. Set status.  */

      *jstat = 0;
      goto out;


   /* Wrong number of intersection points.  */

   err138 : *jstat = -138;
   goto out;

   /* Error in lower level routine.  */

   error : *jstat = kstat;
   goto out;

   out :
      return;
}
