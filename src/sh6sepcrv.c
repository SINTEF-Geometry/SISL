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
 * $Id: sh6sepcrv.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6SEPCRV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh6sepcrv_s9circle(double [], double[], double[], double, 
				  double[], double[], double *, int *);
#else
static void sh6sepcrv_s9circle();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh6sepcrv (SISLCurve *pc1, SISLCurve *pc2, double aepsge, double ecentre[],
	     double *crad, int *jstat)
#else
void
   sh6sepcrv (pc1, pc2, aepsge, ecentre, crad, jstat)
   SISLCurve *pc1;
   SISLCurve *pc2;
   double aepsge;
   double ecentre[];
   double *crad;
   int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find splitting geometry between two B-spline curves
*              The geometry object is a sphere,
*              and is used to put between two input objects to check if
*              an intersection is possible.
*
*
* INPUT      : pc1      - First curve.
*              pc2      - Second curve.
*              aepsge   - Geometry tolerance. 
*
*
* OUTPUT     : ecentre  - Centre of sphere.
*              crad     - Radius of sphere.
*              jstat    - status messages
*                                = 1   : A sphere is found.
*                                = 0   : Ok. No splitting geometry is found.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : s1421   -  Surface evaluator.
*              s1773   -  Closest point between point and surface.
*              s6ang   -  Angle between vectors.
*              s6scpr  -  Scalar product between vectors.
*              s6norm  -  Normalize vector.
*              s6diff  -  Difference vector.
*              s6dist  -  Distance between points. 
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, 93-02.
* CHANGED BY : Vibeke Skytt, SI, 93-05. New splits introduced.
*
*********************************************************************
*/
{
   int kstat = 0;         /* Status variable.       */
   int ki,kj;             /* Counters.              */
   int kleft1=0; /* Parameters to surface evaluator. */
   int kdim=pc1->idim;    /* Dimension of geometry space.      */
   double tpar;       /* Midpoint of surface.              */
   double sparc1[3];      /* Corner parameters of first surface. */
   double sparc2[3];      /* Parameters of closest points on second surface. */
   double scorn1[9];     /* Corners of first surface.           */
   double scorn2[9];     /* Closest points in the other surface. */
   SISLPoint *qp = SISL_NULL;  /* Representing a surface corner as a point. */
   double tstart;      /* Start parameters of second surface.       */
   double tend;       /* End parameters of second surface.       */
   double saxis[3];       /* Normal to circle between edges.         */
   double tdot;
   double tsign;
   double tang;
   double tpi4 = PI/(double)4.0;
   
   /* Test dimension.  */
   
   if (kdim != 3)
   {
      *jstat = 0;
      goto out;
   }
   
   /* Test if the cones of the surfaces is less than pi, otherwise
      no attempt to find splitting geometry is made.  */
   
   if (pc1->pdir->igtpi != 0 || pc2->pdir->igtpi != 0)
   {
      *jstat = 0;
      goto out;
   }
   
   
   /* Check that the objects are not too large, i.e. contain to many
      vertices to be put into a sphere equation effectively. */
   
   if (pc1->in > 4*pc1->ik || pc2->in > 4*pc2->ik)
   {
      *jstat = 0;
      goto out;
   }
   
   /* Make sure that the cones lies in the same area, otherwise
      return.   */

   tdot = s6scpr(pc1->pdir->ecoef,pc2->pdir->ecoef,kdim);
   tsign = (tdot >= DZERO) ? (double)1.0 : -(double)1.0;

   tang = s6ang(pc1->pdir->ecoef,pc2->pdir->ecoef,kdim);
   if (tang > tpi4)
   {
      *jstat = 0;
      goto out;
   }
 
   /* Try to find a circle splitting the edge curves of the surface and the
      cyrve, and extend this circle to a sphere.  */
   
   sparc1[0] = *(pc1->et+pc1->ik-1);
   sparc1[2] = *(pc1->et+pc1->in);
   sparc1[1] = (double)0.5*(sparc1[0] + sparc1[2]);
   
   tstart = *(pc2->et + pc2->ik - 1);
   tend = *(pc2->et + pc2->in);
   tpar= (double)0.5*(tstart + tend);
   
   for (ki=0; ki<3; ki++)
   {
      /* Evaluate curve.  */
      
      s1221(pc1, 0, sparc1[ki], &kleft1, scorn1+ki*kdim, &kstat);
      if (kstat < 0) goto error;
      
      /* Find the closest point in the other surface. First express
	 the corner as a SISLPoint. */
      
      if ((qp = newPoint(scorn1+ki*kdim, kdim, 1)) == SISL_NULL) goto err101;
      s1771(qp, pc2, aepsge, tstart, tend, tpar, sparc2+ki, &kstat);
      if (kstat < 0) goto error;
      
      /* Evaluate second curve. */
      
      s1221(pc2, 0, sparc2[ki], &kleft1, scorn2+ki*kdim, &kstat);
      if (kstat < 0) goto error;
      
      if (qp != SISL_NULL) freePoint(qp); qp = SISL_NULL;
   }
   
   
   /* Find middle points between the sets of closest points. */
   
   for (kj=0; kj<3; kj++)
      for (ki=0; ki<kdim; ki++)
	 scorn1[kj*kdim+ki] = (double)0.5*(scorn1[kj*kdim+ki] + 
					   scorn2[kj*kdim+ki]);
   
   /* Compute splitting cylinder. */
   
   sh6sepcrv_s9circle(scorn1, scorn1+kdim, scorn1+2*kdim,
		       aepsge, ecentre, saxis, crad, &kstat);
   if (kstat < 0) goto error;
   if (kstat > 0)
   {
      *jstat = 0;
      goto out;
   }
   
   /* Output sphere. */
   
   *jstat = 1;
   goto out;

   err101 : *jstat = -101;
   goto out;
   
   error : *jstat = kstat;
   goto out;
   
   out :
      return;
}

   
#if defined(SISLNEEDPROTOTYPES)
static void
  sh6sepcrv_s9circle(double apt1[], double apt2[], double apt3[],
		     double aepsge, double ecentre[], double eaxis[],
		     double *crad, int *jstat)
#else
static void sh6sepcrv_s9circle(apt1, apt2, apt3, aepsge, ecentre, eaxis,
			       crad, jstat)
      double apt1[];
      double apt2[];
      double apt3[];
      double aepsge;
      double ecentre[];
      double eaxis[];
      double *crad;
      int *jstat;
#endif      
/*
*********************************************************************
*
* PURPOSE    : Compute the circle passing through 3 points in 3D, and
*              the normal to the plane in which the cirle lies.
*
*
* INPUT      : apt1   - First point on the circle.
*              apt2   - Second point on the circle.
*              apt3   - Third point on the circle.
*              aepsge - Geometry tolerance.
*              
*
* OUTPUT     : ecentre - Centre of the circle.
*              eaxis   - Normal to the plane in which the circle lies.
*              crad    - Radius of the circle.
*              jstat   - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* CALLS      : s6diff, s6crss, s6length, s6scpr, s6dist, s6lufacp, s6lusolp 
*
*********************************************************************
*/
{
   int kstat = 0;
   int ki;
   int kdim = 3;
   int lpiv[3];
   double snorm[3];
   double smid1[3];
   double smid2[3];
   double sdiff1[3];
   double sdiff2[3];
   double smat[9];
   double sright[3];
   
   /* Compute difference vectors between the 1. and 2. and 2. and 3. point. */
   
   s6diff(apt1, apt2, kdim, sdiff1);
   s6diff(apt3, apt2, kdim, sdiff2);
   
   /* Compute the normal of the plane in which the circle lies. */
   
   s6crss(sdiff1, sdiff2, snorm);
   
   /* Compute the normals to the planes normal to the first plane and
      perpendicular to the difference vectors. */
   
   /* s6crss(sdiff1, snorm, snorm1);
   s6crss(sdiff2, snorm, snorm3); */
   
   /* Check normals.  */
   
   if (s6norm(sdiff1, kdim, sdiff1, &kstat) < aepsge) goto warn1;
   if (s6norm(snorm, kdim, snorm, &kstat) < aepsge) goto warn1;
   if (s6norm(sdiff2, kdim, sdiff2, &kstat) < aepsge) goto warn1; 
   
   /* Compute the midpoints of the difference vectors. */
   
   for (ki=0; ki<kdim; ki++)
   {
      smid1[ki] = (double)0.5*(apt1[ki] + apt2[ki]);
      smid2[ki] = (double)0.5*(apt2[ki] + apt3[ki]);
   }
   
   /* Set up equation system.  */

   memcopy(smat, snorm, kdim, DOUBLE);
   memcopy(smat+kdim, sdiff1, kdim, DOUBLE);
   memcopy(smat+2*kdim, sdiff2, kdim, DOUBLE);
   
   sright[0] = s6scpr(apt2, snorm, kdim);
   sright[1] = s6scpr(smid1, sdiff1, kdim);
   sright[2] = s6scpr(smid2, sdiff2, kdim);
   
   /* Solve equation system.  */
   
   s6lufacp(smat, lpiv, 3, &kstat);
   if (kstat < 0) goto error;
   
   s6lusolp(smat, sright, lpiv, 3, &kstat);
   if (kstat < 0) goto error;
   
   /* Prepare output.  */
   
   memcopy(eaxis, snorm, kdim, DOUBLE);
   memcopy(ecentre, sright, kdim, DOUBLE);
   *crad = s6dist(ecentre, apt2, kdim);
   
   *jstat = 0; 
   goto out;
   
   /* Almost singular equation system.  */
   
   warn1 :
      *jstat = 1;
   goto out;
   
   /* Error in lower level routine.  */
   
   error :
      *jstat = kstat;
   goto out;
   
   out :
      return;
}
      
