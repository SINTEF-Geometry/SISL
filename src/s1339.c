/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"


#define S1339

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1339(SISLSurf *ps, int ipow, double ea[], int inumprd, int incoef,
	 double aepsge, double ecimp[], double *cgrad, int *jstat)
#else
void s1339(ps, ipow, ea, inumprd, incoef, aepsge, ecimp, cgrad, jstat)
   SISLSurf *ps;
   int    ipow;
   double ea[];
   int    inumprd;
   int    incoef;
   double aepsge;
   double ecimp[];
   double *cgrad;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Solve an underdetermined equation system to find implicit
*              coefficents approximating a Bezier surface.
*
*
* INPUT      : ps       - Surface to approximate.
*              ipow     - Implicit degree.
*              ea       - Matrix of equation system of implicit coefficients
*              inumprd  - Number of products in the matrix ea, i.e. number
*                         of columns in the matrix.
*              incoef   - Number of coefficients in the matrix ea, i.e. number
*                         of rows in the matrix.
*              aepsge   - Geometry tolerance.
*
*
* OUTPUT     : ecimp    - Implicit coefficients.
*              cgrad    - Estimated medium gradient of the implicit surface.
*              jstat    - status messages
*                                = 1   : Implicit approximation not computed.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      :
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SINTEF Oslo, 94-09.
*
*********************************************************************
*/
{
   int kstat = 0;    /* Local status variable.                           */
   int knvec;        /* Number of vectors spanning the null space.       */
   int kntest = 10;  /* Number of sample points in each par. dir.        */
   int kr1, kr2, kl; /* Counters of sample points and solution vectors.  */
   int kl2, kl3;
   int ki, kj, kk;   /* Counters of multiplisity of factors.             */
   int kindex;       /* Index in snullvec.                               */
   int kleft1=0, kleft2=0;  /* Parameters used in surface evaluator.     */
   double tx, ty, tz;       /* Help variables used to compute gradient.  */
   double tx1, ty1, tz1;
   double tg1, tg2, tg3;    /* Help variables used to compute gradient.  */
   double tweight;          /* Weight belonging to a particular solution
			       to the approximation problem.             */
   double spar[2];          /* Parameter value of sampling point.        */
   double tint1, tint2;     /* Parametric intervals between sampling points. */
   double totaldot;         /* Sum of scalar products between gradients of
			       implicit surface and surface normal of
			       parametric surface.                         */
   double tepsge = MAX(aepsge/(double)1000, REL_COMP_RES);
                            /* Tolerance used to sort out negliable vectors
			       of the nullspace.                           */
   double *snullvec = NULL; /* Vectors spanning the null space.            */
   double sder[9];          /* Position and derivatives of spline surface. */
   double snorm[3];         /* Normal of spline surface.                   */
   double *sgrad = NULL;    /* Array storing gradients.                    */
   double *accdot = NULL;   /* Accumulated information about the scalar
			      product between the gradient of the implicit
			      surface and the normal of the parametric
			      surface in a number of points.              */
   double *accgrad = NULL; /* Accumulated information about the gradient
			      of the implicit surface.                    */

   /* Test input. */

   /*if (ps->in1 > ps->ik1 || ps->in2 > ps->ik2) 
     goto err122;*/

   /* Find a number of vectors spanning the null space of the equation
      system for finding the implicit coefficients.                                                    */

   s6nullspace(ea, incoef, inumprd, tepsge, &snullvec, &knvec, &kstat);
   if (kstat < 0) goto error;

   /* Allocate scratch for gradients. */

   if ((sgrad = new0array(5*knvec, DOUBLE)) == NULL) 
     goto err101;
   accgrad = sgrad + 3*knvec;
   accdot = accgrad + knvec;

   /* Accumulate the dot product between the normal of the surface and
      the gradient of the candidate implicit surfaces at a number of
      sample points. Also compute an estimate for the length of the
      medium gradient of the candidate implicit surfaces.              */

   /* Find parametric interval between sampling points. */

   tint1 = (ps->et1[ps->in1] - ps->et1[ps->ik1-1])/(kntest-1);
   tint2 = (ps->et2[ps->in2] - ps->et2[ps->ik2-1])/(kntest-1);

   for (kr1=0, spar[0]=ps->et1[ps->ik1-1]; kr1<kntest; kr1++, spar[0]+=tint1)
      for (kr2=0, spar[1]=ps->et2[ps->ik2-1]; kr2<kntest; kr2++, spar[1]+=tint2)
      {
	 /* Evaluate the surface to be approximated. */

	 s1421(ps, 1, spar, &kleft1, &kleft2, sder, snorm, &kstat);
	 if (kstat < 0) goto error;

	 /* Initialize vector for gradient storage to zero. */

	 memzero(sgrad, 3*knvec, DOUBLE);

	 /* Compute the gradient in the current sample point. */

	 for (kindex=0,tx1=tx=(double)1, ki=0; ki<=ipow; 
	      ki++, tx1=tx, tx*=sder[0])
	    for (ty1=ty=(double)1, kj=0; kj<=ipow-ki; 
		 kj++, ty1=ty, ty*=sder[1])
	       for (tz1=tz=(double)1, kk=0; kk<=ipow-ki-kj; 
		    kk++, tz1=tz, tz*=sder[2])
	       {
		 /* FIX start ------------->                               */
		 /* Fixed by HKE and UJK, Linear term not properly handled */
		 if (ki > 1 && fabs(sder[0]) > REL_COMP_RES)
		   tg1 = (double)ki*tx1*ty*tz;
		 else if (ki == 1)
		   tg1 = ty*tz;
		 else
		   tg1 = 0.0;

		 if (kj > 1 && fabs(sder[1]) > REL_COMP_RES)
		   tg2 = (double)kj*tx*ty1*tz;
		 else if (kj == 1)
		   tg2 = tx*tz;
		 else
		   tg2 = 0.0;

		 if (kk > 1 && fabs(sder[2]) > REL_COMP_RES)
		   tg3 = (double)kk*tx*ty*tz1;
		 else if (kk == 1)
		   tg3 = tx*ty;
		 else
		   tg3 = 0.0;

		 /* FIX end <-------------                                 */


		  /* Add contribution to the gradient. */
		  for (kl3=kl2=kl=0; kl<knvec; kl++, kl2+=inumprd, kl3+=3)
		  {
		     sgrad[kl3] += snullvec[kl2+kindex]*tg1;
		     sgrad[kl3+1] += snullvec[kl2+kindex]*tg2;
		     sgrad[kl3+2] += snullvec[kl2+kindex]*tg3;
		  }
		  kindex++;
	       }

	 /* Add information of the current sample point to
	    the accumulated scalar product between gradients and
	    normals, and to the accumulated lenghts of gradients. */

	 for (kl3=kl=0; kl<knvec; kl++, kl3+=3)
	 {
	    accdot[kl] += s6scpr(sgrad+kl3, snorm, 3);
	    accgrad[kl] += s6length(sgrad+kl3, 3, &kstat);
	 }
      }

   /* Compute sum of scalar products between gradients and normals,
      and compute medium of estimated gradient.     */

   for (totaldot=0.0, kl=0; kl<knvec; kl++)
   {
      totaldot += accdot[kl]*accdot[kl];
      accgrad[kl] /= kntest*kntest;
   }
   totaldot = sqrt(totaldot);

   /* Initiate output storage to zero. */

   memzero(ecimp, inumprd, DOUBLE);
   *cgrad = 0.0;

   /* Compute factors of each candidate solution to the set of implicit
      coefficients, and make the final sulution vector.                 */

   if (totaldot < REL_COMP_RES) goto warn1;  /* Implicitation not possible. */

   for (kl2=kl=0; kl<knvec; kl++, kl2+=inumprd)
   {
      tweight = fabs(accdot[kl])/totaldot;
      for (kindex=0; kindex<inumprd; kindex++)
	 ecimp[kindex] += tweight*snullvec[kl2+kindex];
      (*cgrad) += tweight*accgrad[kl];
   }

   /* Implicit approximation computed. */

   *jstat = 0;
   goto out;

   /* Implicit approximation not computed. */

   warn1 : *jstat = 1;
   goto out;

   err101 : *jstat = -101;   /* Error in scratch allocation.   */
   s6err("s1339",*jstat,0);
   goto out;

   error : *jstat = kstat;   /* Error in lower level routine. */
   s6err("s1339",*jstat,0);
   goto out;

   out:
      /* Free space occupied by local arrays. */

      if (snullvec != NULL) freearray(snullvec);
      if (sgrad != NULL) freearray(sgrad);

      return;
}
