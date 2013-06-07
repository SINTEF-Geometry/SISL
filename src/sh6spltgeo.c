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
 * $Id: sh6spltgeo.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH6SPLITGEOM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh6splitgeom_s9circle(double [], double[], double[], double, 
				  double[], double[], double *, int *);
#else
static void sh6splitgeom_s9circle();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh6splitgeom (SISLSurf *ps1, SISLSurf *ps2, double aepsge, double ecentre[],
		double eaxis[], double *cdist, double *crad, int *jstat)
#else
void
sh6splitgeom (ps1, ps2, aepsge, ecentre, eaxis, cdist, crad, jstat)
   SISLSurf *ps1;
   SISLSurf *ps2;
   double aepsge;
   double ecentre[];
   double eaxis[];
   double *cdist;
   double *crad;
   int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find splitting geometry between two B-spline surfaces.
*              The geometry object is either a plane, a sphere, a cylinder
*              or a torus and is used to put between two surfaces to check 
*              if an intersection is possible.
*
*
* INPUT      : ps1      - First surface.
*              ps2      - Second surface.
*              aepsge   - Geometry tolerance. 
*
*
* OUTPUT     : ecentre  - Centre of sphere/cylinder/torus or point in plane.
*              eaxis    - Axis of cylinder/torus or plane normal.
*              cdist    - Large radius of torus.
*              crad     - Radius of sphere/cylinder or small radius of torus.
*              jstat    - status messages
*                                = 4   : A torus is found.
*                                = 3   : A cylinder is found.
*                                = 2   : A sphere is found.
*                                = 1   : A plane is found.
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
* CHANGED BY : Vibeke Skytt, SI, 93-05. New geometry introduced.
*
*********************************************************************
*/
{
   int kstat = 0;         /* Status variable.       */
   int ki,kj,k1,k2;       /* Counters.              */
   int kleft1=0, kleft2=0; /* Parameters to surface evaluator. */
   int kder=0;            /* Evaluate only position.           */
   int kdim=ps1->idim;    /* Dimension of geometry space.      */
   double tpi6=PI/(double)6;
   double tsign;          /* Sign of vector.        */
   double tdot;           /* Scalar product.        */
   double tang;           /* Angle between vectors. */
   double trad1, trad2;   /* Curvature radius.          */
   double tmaxrad;        /* Maximum radius of sphere/torus/cylinder. */
   double tminfac = (double)0.9; /* Minimum factor between radiuses 
				    for a sphere. */
   double spar1[2],spar2[2]; /* Paramater value in which to evaluate
				surfaces.                         */
   double sder1[18];      /* Value of 1. surface. */
   double snorm1[3];      /* Normal of 1. surface.  */
   double sder2[18];      /* Value of 2. surface. */
   double snorm2[3];      /* Normal of 2. surface.  */
   double scentre1[3];    /* Centre of first circle. */
   double scentre2[3];    /* Centre of second circle. */
   double svec[3];        /* Vector used to find midpoint of 
			     splitting geometry.    */
   double sdiff[3];       /* Difference vector between midpoints. */
   double sparc1[10];      /* Corner parameters of first surface. */
   double sparc2[10];      /* Parameters of closest points on second surface. */
   double scorn1[15];     /* Corners of first surface.           */
   double scorn2[15];     /* Closest points in the other surface. */
   SISLPoint *qp = SISL_NULL;  /* Representing a surface corner as a point. */
   double start2[2];      /* Start parameters of second surface.       */
   double send2[2];       /* End parameters of second surface.       */
   double sdist[4];       /* Distance between closest points.        */
   
   /* Test if the cones of the surfaces is less than pi, otherwise
      no attempt to find splitting geometry is made.  */
   
   if (ps1->pdir->igtpi != 0 || ps2->pdir->igtpi != 0)
   {
      *jstat = 0;
      goto out; 
   }
   
   /*   if (ps1->pdir->aang > tpi4 || ps2->pdir->aang > tpi4)
   {
    *jstat = 0;
      goto out; 
   } */
   
   /* Make sure that the cones lies in the same area, otherwise
      return.   */
   
   tdot = s6scpr(ps1->pdir->ecoef,ps2->pdir->ecoef,kdim);
   tsign = (tdot >= DZERO) ? (double)1.0 : -(double)1.0;
   
   tang = s6ang(ps1->pdir->ecoef,ps2->pdir->ecoef,kdim);
   if (tang > tpi6)
   {
      *jstat = 0;
      goto out; 
   } 
   
   /* Check that the surfaces is not too large, i.e. contain to many
      vertices to be put into a sphere- or cylinder equation effectively. */
   
   if (ps1->in1 > 2*ps1->ik1 || ps1->in2 > 2*ps1->ik2 ||
       ps2->in1 > 2*ps2->ik1 || ps2->in2 > 2*ps2->ik2)
   {
      *jstat = 0;
      goto out; 
   } 
   
   /* Compute the midvector between the axises of the surface cones. */
   
   for (ki=0; ki<kdim; ki++)
      svec[ki] = (double)0.5*(tsign*ps1->pdir->ecoef[ki] + ps2->pdir->ecoef[ki]);
   (void)s6norm(svec,kdim,svec,&kstat);
   if (!kstat)
   {
      *jstat = 0;
      goto out; 
   }
   
   /* Set maximum radius. */
   
   tmaxrad = ps1->pbox->e2max[2][0] - ps1->pbox->e2min[2][0];
   tmaxrad = MAX(tmaxrad, ps1->pbox->e2max[2][1]-ps1->pbox->e2min[2][1]);
   tmaxrad = MAX(tmaxrad, ps1->pbox->e2max[2][2]-ps1->pbox->e2min[2][2]);
   tmaxrad *= (double)10.0;
   
   /* Set parameter bourders of second surface. */
   
   start2[0] = *(ps2->et1 + ps2->ik1 - 1);
   start2[1] = *(ps2->et2 + ps2->ik2 - 1);
   send2[0] = *(ps2->et1 + ps2->in1);
   send2[1] = *(ps2->et2 + ps2->in2);
   
   /* Evaluate the surfaces in their midpoints up to 2. order
      derivatives.                                             */
   
   spar1[0] = (double)0.5*(ps1->et1[ps1->ik1-1] + ps1->et1[ps1->in1]);
   spar1[1] = (double)0.5*(ps1->et2[ps1->ik2-1] + ps1->et2[ps1->in2]);
      
   s1421(ps1,kder,spar1,&kleft1,&kleft2,sder1,snorm1,&kstat);
   if (kstat < 0) goto error;
   
   spar2[0] = (double)0.5*(ps2->et1[ps2->ik1-1] + ps2->et1[ps2->in1]);
   spar2[1] = (double)0.5*(ps2->et2[ps2->ik2-1] + ps2->et2[ps2->in2]);
      
   s1421(ps2,kder,spar2,&kleft1,&kleft2,sder2,snorm2,&kstat);
   if (kstat < 0) goto error;
   
   /* Check if the difference vector between the midpoints point in 
      about the same direction as the vector svec.  */
   
   s6diff(sder1, sder2, kdim, sdiff);
   tang = s6ang(sdiff, svec, kdim);
   if (tang < tpi6 || tang > (double)5.0*tpi6)
   {
      /* Set up parameter values for evaluation of first surface in
	 the midpoint and in the midpoints of each edge curve.       */
      
      memcopy(sparc1, spar1, 2, DOUBLE);
      sparc1[3] = *(ps1->et2+ps1->ik2-1);
      sparc1[4] = *(ps1->et1+ps1->in1);
      sparc1[7] = *(ps1->et2+ps1->in2);
      sparc1[8] = *(ps1->et1+ps1->ik1-1);
      sparc1[2] = sparc1[6] = spar1[0];
      sparc1[5] = sparc1[9] = spar1[1];
      
      for (ki=0; ki<5; ki++)
      {
	 /* Evaluate point.  */
	 
	 if (ki == 0)
	    memcopy(scorn1, sder1, kdim, DOUBLE);
	 else
	 {
	    s1421(ps1, 0, sparc1+2*ki, &kleft1, &kleft2, scorn1+ki*kdim,
		  snorm1, &kstat);
	    if (kstat < 0) goto error;
	 }
	 
	 /* Find the closest point in the other surface. First express
	    the corner as a SISLPoint. */
	 
	 if ((qp = newPoint(scorn1+ki*kdim, kdim, 1)) == SISL_NULL) goto err101;
	 s1773(qp, ps2, aepsge, start2, send2, spar2, sparc2+2*ki, &kstat);
	 if (kstat < 0) goto error;
	 
	 /* Evaluate surface. */
	 
	 s1421(ps2, 0, sparc2+2*ki, &kleft1, &kleft2, scorn2+ki*kdim,
	       snorm2, &kstat);
	 if (kstat < 0) goto error;
	 
	 /* Compute midpoint. */

	 for (kj=0; kj<kdim; kj++)
	    scorn1[ki*kdim+kj] = (double)0.5*(scorn1[ki*kdim+kj] + scorn2[ki*kdim+kj]);
	    
	 if (qp != SISL_NULL) freePoint(qp); qp = SISL_NULL;
      }
      
      /* Estimate circles.  */
      
      sh6splitgeom_s9circle(scorn1+kdim, scorn1, scorn1+3*kdim,
			    aepsge, scentre1, snorm1, &trad1, &kstat);
      if (kstat < 0) goto error;
      if (kstat > 0)
	 *jstat = 1;  /* Find plane. */

      sh6splitgeom_s9circle(scorn1+4*kdim, scorn1, scorn1+2*kdim,
			    aepsge, scentre2, snorm2, &trad2, &kstat);
      if (kstat < 0) goto error;
      if (kstat > 0)
	 *jstat = 1;  /* Find plane. */
      
      /* Find kind of splitting geometry.  */
      
      if (*jstat == 1 || (trad1 > tmaxrad && trad2 > tmaxrad))
      {
	 /* Set plane geometry. */
	 
	 *jstat = 1;
	 memcopy(ecentre, scorn1, kdim, DOUBLE);
	 s6diff(scorn1+2*kdim, scorn1, kdim, scorn1+2*kdim);
	 s6diff(scorn1+3*kdim, scorn1, kdim, scorn1+3*kdim);
	 s6crss(scorn1+2*kdim, scorn1+3*kdim, eaxis);
      }
      else if (MAX(trad1,trad2) > tmaxrad)
      {
	 /* Set cylinder geometry.  */
	 
	 *jstat = 3;
	 *crad = MIN(trad1, trad2);
	 if (trad1 < trad2)
	 {
	    memcopy(ecentre, scentre1, kdim, DOUBLE);
	    memcopy(eaxis, snorm1, kdim, DOUBLE);
	 }
	 else
	 {
	    memcopy(ecentre, scentre2, kdim, DOUBLE);
	    memcopy(eaxis, snorm2, kdim, DOUBLE);
	 }
      }
      else if (MIN(trad1,trad2)/MAX(trad1,trad2) > tminfac)
      {
	 /* Set sphere geometry. */
	 
	 *jstat = 2;
	 *crad = (double)0.5*(trad1 + trad2);
	 for (kj=0; kj<kdim; kj++) 
	    ecentre[kj] = (double)0.5*(scentre1[kj] + scentre2[kj]);
      }
      else if (MAX(trad1,trad2)/MIN(trad1,trad2) > (double)25.0)
      {
	 /* Little chance of success in interception. */
	 
	 *jstat = 0;
	 goto out;
      }
      else
      {
	 /* Set torus geometry.  */
	 
	 *jstat = 4;
	 *crad = MIN(trad1, trad2);
	 *cdist = MAX(trad1, trad2) - (*crad);
	 *crad = MIN(trad1, trad2);
	 if (trad1 < trad2)
	 {
	    memcopy(ecentre, scentre2, kdim, DOUBLE);
	    memcopy(eaxis, snorm2, kdim, DOUBLE);
	 }
	 else
	 {
	    memcopy(ecentre, scentre1, kdim, DOUBLE);
	    memcopy(eaxis, snorm1, kdim, DOUBLE);
	 }
	     
      }
   }
   else if (tang > (double)2.0*tpi6 && tang < (double)4.0*tpi6)
   {
      /* Try to find a circle splitting the edge curves of the surfaces,
	 and extend this circle to a cylinder. First find closest edgecurves
	 by feching the corners of the first surface and finding the closest
	 points in the other surface. */
      
      sparc1[6] = sparc1[0] = *(ps1->et1+ps1->ik1-1);
      sparc1[2] = sparc1[4] = *(ps1->et1+ps1->in1);
      sparc1[1] = sparc1[3] = *(ps1->et2+ps1->ik2-1);
      sparc1[5] = sparc1[7] = *(ps1->et2+ps1->in2);
      
      for (ki=0; ki<4; ki++)
      {
	 /* Evaluate corner.  */
	 
	 s1421(ps1, 0, sparc1+2*ki, &kleft1, &kleft2, scorn1+ki*kdim,
	       snorm1, &kstat);
	 if (kstat < 0) goto error;
	 
	 /* Find the closest point in the other surface. First express
	    the corner as a SISLPoint. */
	 
	 if ((qp = newPoint(scorn1+ki*kdim, kdim, 1)) == SISL_NULL) goto err101;
	 s1773(qp, ps2, aepsge, start2, send2, spar2, sparc2+2*ki, &kstat);
	 if (kstat < 0) goto error;
	 
	 /* Evaluate surface. */
	 
	 s1421(ps2, 0, sparc2+2*ki, &kleft1, &kleft2, scorn2+ki*kdim,
	       snorm2, &kstat);
	 if (kstat < 0) goto error;
	 
	 /* Compute distance. */
	 
	 sdist[ki] = s6dist(scorn1+ki*kdim, scorn2+ki*kdim, kdim);
	 
	 if (qp != SISL_NULL) freePoint(qp); qp = SISL_NULL;
      }
      
      /* Check if the two closest points lies on a common edge. */
      
      if (sdist[0] < MIN(sdist[1],sdist[2]) && 
	  sdist[3] < MIN(sdist[1],sdist[2]))
      {
	 k1 = 0; k2 = 3;
      }
      else if (sdist[1] < MIN(sdist[0],sdist[3]) && 
	       sdist[2] < MIN(sdist[0],sdist[3]))
      {
	 k1 = 1; k2 = 2;
      }
      else if (sdist[0] < MIN(sdist[2],sdist[3]) && 
	       sdist[1] < MIN(sdist[2],sdist[3]))
      {
	 k1 = 0; k2 = 1;
      }
      else if (sdist[2] < MIN(sdist[0],sdist[1]) && 
	       sdist[3] < MIN(sdist[0],sdist[1]))
      {
	 k1 = 2; k2 = 3;
      }
      else
      {
	 *jstat = 0;
	 goto out;
      }
      
      /* Compute closest point to the midpoint between the two closest
	 corners.                                                       */
      
      sparc1[8] = (double)0.5*(sparc1[2*k1] + sparc1[2*k2]);
      sparc1[9] = (double)0.5*(sparc1[2*k1+1] + sparc1[2*k2+1]);
      
      /* Evaluate point.  */
      
      s1421(ps1, 0, sparc1+8, &kleft1, &kleft2, scorn1+4*kdim,
	    snorm1, &kstat);
      if (kstat < 0) goto error;
      
      /* Find the closest point in the other surface. First express
	 the corner as a SISLPoint. */
      
      if ((qp = newPoint(scorn1+4*kdim, kdim, 1)) == SISL_NULL) goto err101;
      s1773(qp, ps2, aepsge, start2, send2, spar2, sparc2+8, &kstat);
      if (kstat < 0) goto error;
      
      if (qp != SISL_NULL) freePoint(qp); qp = SISL_NULL;
      
      /* Evaluate surface. */
      
      s1421(ps2, 0, sparc2+8, &kleft1, &kleft2, scorn2+4*kdim,
	    snorm2, &kstat);
      if (kstat < 0) goto error;

      /* Find middle points between the sets of closest points. */
      
      for (ki=0; ki<kdim; ki++)
      {
	 scorn1[k1*kdim+ki] = (double)0.5*(scorn1[k1*kdim+ki] + scorn2[k1*kdim+ki]);
	 scorn1[k2*kdim+ki] = (double)0.5*(scorn1[k2*kdim+ki] + scorn2[k2*kdim+ki]);
	 scorn1[4*kdim+ki] = (double)0.5*(scorn1[4*kdim+ki] + scorn2[4*kdim+ki]);
      }
      
      /* Compute splitting cylinder. */
      
      sh6splitgeom_s9circle(scorn1+k1*kdim, scorn1+4*kdim, scorn1+k2*kdim,
			    aepsge, ecentre, eaxis, crad, &kstat);
      if (kstat < 0) goto error;
      if (kstat > 0 || *crad > tmaxrad)
      {
	 /* Make plane. */
	    
	 *jstat = 1;
	 memcopy(ecentre, scorn1+4*kdim, kdim, DOUBLE);
	 s6diff(scorn1+k1*kdim, scorn1+4*kdim, kdim, scorn1+k1*kdim);
	 s6diff(scorn1+k2*kdim, scorn1+4*kdim, kdim, scorn1+k2*kdim);
	 s6crss(scorn1+k1*kdim, scorn1+k2*kdim, eaxis);
      }
      else
	 
	 /* Output cylinder. */
	 
	 *jstat = 3;
   }
   else *jstat = 0;
   
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
  sh6splitgeom_s9circle(double apt1[], double apt2[], double apt3[],
			double aepsge, double ecentre[], double eaxis[],
			double *crad, int *jstat)
#else
   static void sh6splitgeom_s9circle(apt1, apt2, apt3, aepsge, ecentre, eaxis,
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
   double snorm2[3];
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
   
   s6crss(sdiff1, sdiff2, snorm2);
   
   /* Compute the normals to the planes normal to the first plane and
      perpendicular to the difference vectors. */
   
   /* s6crss(sdiff1, snorm2, snorm1);
   s6crss(sdiff2, snorm2, snorm3); */
   
   /* Check normals.  */
   
   if (s6norm(snorm2, kdim, snorm2, &kstat) < aepsge) goto warn1;
   
   /* Compute the midpoints of the difference vectors. */
   
   for (ki=0; ki<kdim; ki++)
   {
      smid1[ki] = (double)0.5*(apt1[ki] + apt2[ki]);
      smid2[ki] = (double)0.5*(apt2[ki] + apt3[ki]);
   }
   
   /* Set up equation system.  */

   memcopy(smat, snorm2, kdim, DOUBLE);
   memcopy(smat+kdim, sdiff1, kdim, DOUBLE);
   memcopy(smat+2*kdim, sdiff2, kdim, DOUBLE);
   
   sright[0] = s6scpr(apt2, snorm2, kdim);
   sright[1] = s6scpr(smid1, sdiff1, kdim);
   sright[2] = s6scpr(smid2, sdiff2, kdim);
   
   /* Solve equation system.  */
   
   s6lufacp(smat, lpiv, 3, &kstat);
   if (kstat < 0) goto error;
   
   s6lusolp(smat, sright, lpiv, 3, &kstat);
   if (kstat < 0) goto error;
   
   /* Prepare output.  */
   
   memcopy(eaxis, snorm2, kdim, DOUBLE);
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
      
