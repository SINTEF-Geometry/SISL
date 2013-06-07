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

#define S1789

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
static void s1789_s9eval(double [], double [], double [],double [], int, int *);
static int  s1789_s9knot(double [], int, int, double, double, int *, int *);
#else
static void s1789_s9eval();
static int  s1789_s9knot();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
      s1789(SISLPoint *ppoint,SISLSurf *psurf,double aepsge,
	   double epar1[],double epar2[],int *jstat)
#else
   void s1789(ppoint,psurf,aepsge,epar1,epar2,jstat)
      SISLPoint  *ppoint;
      SISLSurf   *psurf;
      double     aepsge;
      double     epar1[];
      double     epar2[];
      int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Test if a point and a surface coincide beetween
*              two intersection points.
*              This function is used when the partial derivatives
*              of the surface are matching for both
*              intersection points.
*
*
* INPUT      : ppoint   - Pointer to the point.
*              psurf    - Pointer to the surface.
*              aepsge   - Geometry resolution.
*              epar1[3] - Parameter values for the first intersection point.
*              epar2[3] - Parameter values for the second intersection point.
*
*
*
* OUTPUT     :  jstat   - status messages
*                                = 1   : Coincidence.
*                                = 0   : No coincidence.
*                                < 0   : Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SINTEF Oslo, Norway. 10.94
* REVISED BY : Vibeke Skytt, 03.98.  Use test also in 3D case.
*
*********************************************************************
*/
{
   int kstat;          /* Status variable                                 */
   int ki;             /* Counter.                                        */
   int kleft1=0;       /* Left indicator for point calculation in 1. par.
			  direction of surface.                           */
   int kleft2=0;       /* Left indicator for point calculation in 2. par dir.*/
   int kknot1, kknot2; /* Indicates whether there is a knot between the
			  input points in 1. and 2. parameter direction.  */
   int kmy1, kmy2;     /* Index of an eventual knot.                      */
   int kk1,kk2,kn1,kn2;/* Orders and nu,ber of vertices of surface        */
   int kdims;          /* Dimension of space where the surface lies       */
   int kpos=0;         /* Position of error                               */
   int kders=2;        /* Number of derivatives to be calculated on surface
			  If step lenght is to be generated from surface,
			  kders must be equal to 2.                       */
   int kpar;           /* Parameter value of constant parameter curve.    */
   double snorm[3];    /* Normal vector of surface                        */
   double *st1;        /* First knot direction of surface                 */
   double *st2;        /* Second knot direction of surface                */
   double sders[18];   /* Position, first and second derivatives of surface */
   double tstep;       /* Final step length     */
   double tlengthend;  /* Length of 1st derivative at end of segment */
   double tincre;      /* Parameter value increment */
   double tsmax;       /* Local maximal step length based of boxsizes of objects */
   double tdist;       /* Distance */
   double tref;        /* Referance value in equality test.               */
   double sstart[2];   /* Lower boundary of parameter intervals */
   double send[2];     /* Upper bounadry of parameter intervals */
   double spos[2];     /* New iteration  point on surface                 */
   double spos1[2];    /* New iteration  point on surface                 */
   double spos2[2];    /* New iteration  point on surface                 */
   double sint[2];     /* Interval between test points in par. space.     */
   double snext[2];    /* Save previous intersection point.               */
   double sdiff[2];    /* Difference vector between input int. pts.       */
   double spardir[2];  /* Direction of coincidence curve in parameter area. */
   double tbeta;       /* Scaling factor between partial derivatives of sf. */
   double stanc[2];    /* Direction of coincidence curve in surface.        */
   double sder2[10];   /* Information about curve in surface.               */
   double tdot;        /* Scalar product to test direction of vectors.      */
   double td;          /* Distance between current and last point.          */
   double s3dinf2[10]; /* Marching information to decide step length.       */
   SISLCurve *qc = SISL_NULL;   /* Constant parameter curve.                       */

   *jstat = 0;

   /* Make maximal step length based on box-size of surface */

   sh1992su(psurf,0,aepsge,&kstat);
   if (kstat < 0) goto error;

   tsmax = MAX(psurf->pbox->e2max[0][0] - psurf->pbox->e2min[0][0],
	       psurf->pbox->e2max[0][1] - psurf->pbox->e2min[0][1]);

   /* Copy surface attributes to local parameters.  */

   kdims = psurf -> idim;
   kk1   = psurf -> ik1;
   kk2   = psurf -> ik2;
   kn1   = psurf -> in1;
   kn2   = psurf -> in2;
   st1   = psurf -> et1;
   st2   = psurf -> et2;

   /* Set reference value.  */

   tref = MAX(st1[kn1]-st1[kk1-1],st2[kn2]-st2[kk2-1]);

   /* Check dimension  */

   if (ppoint->idim != kdims || (kdims != 2 && kdims != 3))
     goto err105;

   sstart[0] = st1[kk1-1];
   sstart[1] = st2[kk2-1];
   send[0] = st1[kn1];
   send[1] = st2[kn2];

   /* Set start point for marching on surface */

   spos1[0] = epar1[0];
   spos1[1] = epar1[1];

   /* Set difference vector between input points. */

   s6diff(epar2, epar1, 2, sdiff);

   /* Evaluate start point of surface.  */

   s1421(psurf,kders,spos1,&kleft1,&kleft2,sders,snorm,&kstat);
   if (kstat < 0) goto error;

   /* While end not reached */

   td = s6dist(spos1, epar2, 2);
   while (td > REL_PAR_RES)
   {
      /* Compute direction of marching. The partial derivatives of the
	 surface in this point must be almost parallel. Find the factor
	 that makes the partial derivatives sum up to zero (approximately). */

     if (kdims == 2)
       {
	 if (DEQUAL(sders[kdims]+tref,tref) && 
	     DEQUAL(sders[kdims+1]+tref,tref) &&
	     DEQUAL(sders[kdims+2]+tref,tref)) break;

	 if (sders[2] >= sders[3])
	   {
	     if (DEQUAL(sders[4]+tref,sders[2]+tref))
	       tbeta = (double)0.5;
	     else
	       tbeta = (double)1/((double)1 - (sders[4]/sders[2]));
	   }
	 else
	   {
	     if (DEQUAL(sders[5]+tref,sders[3]+tref))
	       tbeta = (double)0.5;
	     else
	       tbeta = (double)1/((double)1 - (sders[5]/sders[3]));
	   } 


	 spardir[0] = (double)1-tbeta;
	 spardir[1] = tbeta;
       }
     else
       {
	 spardir[0] = epar2[0]-epar1[0];
	 spardir[1] = epar2[1]-epar1[1];
       }

      tdot = s6norm(spardir, 2, spardir,&kstat);
      if (tdot < REL_PAR_RES)
      {
	 *jstat = 0;
	 goto out;
      }

      for (ki=0; ki<kdims; ki++)
	 stanc[ki] = spardir[0]*sders[kdims+ki] + spardir[1]*sders[2*kdims+ki];

      tdot = s6scpr(stanc, sdiff, kdims);
      if (tdot < DZERO)
      {
	 stanc[0] *= -(double)1;
	 stanc[1] *= -(double)1;
      }

      /* Compute position, first and second derivative of the curve in the
	 surface going through the evaluated point in this point. */

      s1789_s9eval(sders,snorm,stanc,sder2,kdims,&kstat);
      if (kstat < 0) goto error;

      /* Calculate unit tangent and radius of curvature of curve in surface.*/

      s1307(sder2,kdims,s3dinf2,&kstat);
      if (kstat<0) goto error;

      /* Calculate step length based on curvature */

      tstep = s1311(s3dinf2[3*kdims],aepsge,tsmax,&kstat);
      if (kstat<0) goto error;

      tlengthend = s6length(sder2+kdims,kdims,&kstat);
      if (kstat<0) goto error;

      /* Find candidate end point, make sure that no breaks in tangent or
	 curvature exists between start and endpoints of the segment      */

      /* Make step length equal to resolution if the length is zero */

      /* Find parameter value of candidate end point of segment */

      if (DEQUAL(tlengthend+tref,tref))
	 tincre = REL_PAR_RES;
      else
	 tincre = tstep/tlengthend;

      spos2[0] = spos1[0] + tincre*spardir[0];
      spos2[1] = spos1[1] + tincre*spardir[1];

     /* Make sure not to jump out of the surface */
     if ((epar2[0] > epar1[0] && spos2[0] >= epar2[0]) ||
	 (epar2[0] < epar1[0] && spos2[0] <= epar2[0]) ||
	 (epar2[1] > epar1[1] && spos2[1] >= epar2[1]) ||
	 (epar2[1] < epar1[1] && spos2[1] <= epar2[1]))
       {
	 spos2[0] = epar2[0];
	 spos2[1] = epar2[1];
       }

      if (s6dist(spos1, spos2, kdims) > s6dist(spos1, epar2, kdims))
	 memcopy(spos2, epar2, 2, DOUBLE);

      /* Check if any knot line exist within the step. */

      kknot1 = s1789_s9knot(st1, kk1, kn1, spos1[0], spos2[0], &kmy1, &kstat);
      if (kstat < 0) goto error;

      kknot2 = s1789_s9knot(st2, kk2, kn2, spos1[1], spos2[1], &kmy2, &kstat);
      if (kstat < 0) goto error;

      if ((kknot1 && !kknot2) ||
	  (kknot1 && kknot2 && spardir[1]*(st1[kmy1]-spos1[0]) <
	   spardir[0]*(st2[kmy2]-spos1[1])))
      {
	 /* Pull back to knotline in first parameter direction. */

	 spos2[0] = psurf->et1[kmy1];   /* Parameter value of knotline. */
	 spos2[1] = spos1[1] + (spos2[0]-spos1[0])*spardir[1]/spardir[0];
	 kpar = 1;
      }
      else if (kknot2)
      {
	 /* Pull back to knot line in second parameter direction. */

	 spos2[1] = psurf->et2[kmy2];
	 spos2[0] = spos1[0] + (spos2[1] - spos1[1])*spardir[0]/spardir[1];
	 kpar = 2;
      }
      else
      {
	 /* No knot line. Decide in which parameter direction to iterate. */

	 if (spardir[1]*fabs(st1[kmy1]-spos1[0]) <
	     spardir[0]*fabs(st2[kmy2]-spos1[1]))
	    kpar = 1;
	 else
	    kpar = 2;
      }

      sint[0] = (spos2[0]-spos1[0])/(double)3;
      sint[1] = (spos2[1]-spos1[1])/(double)3;

      for (ki=0, spos[0]=spos1[0]+sint[0], spos[1]=spos1[1]+sint[1];
       ki<3; ki++, spos[0]+=sint[0], spos[1]+=sint[1])
      {

	 if (kpar == 1)
	 {
	    /* Pick constant parameter curve in 1. par. dir. */

	    s1437(psurf, spos[0], &qc, &kstat);
	    if (kstat < 0) goto error;

	    /* Iterate down to the curve. */

	    s1771(ppoint, qc, aepsge, qc->et[qc->ik-1], qc->et[qc->in],
		  spos[1], &spos[1], &kstat);
	    if (kstat < 0) goto error;
	 }
	 else
	 {
	    /* Pick constant parameter curve in 2. par. dir. */

	    s1436(psurf, spos[1], &qc, &kstat);
	    if (kstat < 0) goto error;

	    /* Iterate down to the curve. */

	    s1771(ppoint, qc, aepsge, qc->et[qc->ik-1], qc->et[qc->in],
		  spos[0], &spos[0], &kstat);
	    if (kstat < 0) goto error;
	 }

	 memcopy(snext, spos, 2, DOUBLE);

	 /* Calculate point and derivatives in surface */

	 s1421(psurf,kders,spos,&kleft1,&kleft2,sders,snorm,&kstat);
	 if (kstat<0) goto error;

	 /* Check if the input point and surface point are within positional
	    tolerance. */

	 tdist = s6dist(ppoint->ecoef,sders,kdims);

	 if (tdist>aepsge)
	 {
	    /* Points not within tolerances, no coincide. */

	    goto war01;
	 }

	 /* Test whether the marching has advanced. */

	 if (s6dist(spos1, spos, 2) < REL_PAR_RES) goto war01;

         /* Free memory occupied by local curve. */

        if (qc != SISL_NULL) freeCurve(qc);
        qc = SISL_NULL;
      }

      /* Update start parameter of step. */

      spos1[kpar-1] = spos2[kpar-1];
      spos1[2-kpar] = snext[2-kpar];
      td = s6dist(spos1, epar2, 2);
   }

   if (td > REL_PAR_RES) *jstat = 0;
   else *jstat = 1;

   goto out;

   /* Point and surface not within tolerance */
   war01: *jstat = 0;
   goto out;

   /* Error in input, dimension not equal to 2 or 3 */

   err105: *jstat = -105;
   s6err("s1789",*jstat,kpos);
   goto out;

   /* Error in lower level function */

   error:  *jstat = kstat;
   s6err("s1789",*jstat,kpos);
   goto out;


   out:
      if (qc != SISL_NULL) freeCurve(qc);

      return;
}


#if defined(SISLNEEDPROTOTYPES)
   static void
            s1789_s9eval(double eders[],double enorms[],double etanc[],
			 double ederc[],int idim, int *jstat)
#else
	       static void s1789_s9eval(eders,enorms,etanc,ederc,idim,jstat)
		  double eders[];
		  double enorms[];
		  double etanc[];
		  double ederc[];
		  int idim;
		  int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Compute the position, first and second derivative of
*              a curve going through a given point of the surface when
*              the 0-2'th derivatives of the surface is given. The
*              tangent of the wanted curve is parallel to the projection
*              of a given vector into the tangent plane of the surface.
*
*
* INPUT      : eders    - 0-2'th derivatives of the surface. Dimension
*                         is 6*idim.
*              enorms   - Normal vector of the surface. Dimension is idim.
*              etanc    - Vector to be projected into the tangent plane
*                         of the surface. Dimension is idim.
*              idim     - Dimension of geometry space.
*
*
*
* OUTPUT     : ederc   - 0-2'th derivative of the curve in the surface.
*              jstat   - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, Oslo, Norway. Oct. 1990
*
*********************************************************************
*/
{
   int kstat = 0;         /* Status variable.  */
   int ki;                /* Counter.          */
   int ksign = 1;         /* Parameter used in s6findfac.     */
   double tfac1,tfac2,tfac3;  /* Factors found by s6findfac.  */

   /* Copy position of surface to output array.   */

   memcopy(ederc,eders,idim,DOUBLE);

   /* Compute the factors used to express etanc by the derivatives and normal
      of the surface.  */

   s6findfac(eders+idim,eders+2*idim,enorms,etanc,idim,ksign,&tfac1,&tfac2,
	     &tfac3,&kstat);
   if (kstat < 0) goto error;

   /* Compute first and second derivative of the curve in the surface.  */

   for (ki=0; ki<idim; ki++)
   {
      ederc[idim+ki] = tfac1*eders[idim+ki] + tfac2*eders[2*idim+ki];
      ederc[2*idim+ki] = tfac1*tfac1*eders[3*idim+ki]
	 + (double)2.0*tfac1*tfac2*eders[4*idim+ki] + tfac2*tfac2*eders[5*idim+ki];
   }

   *jstat = 0;
   goto out;

   /* Error in lower level routine.  */

   error:
      *jstat = kstat;
   goto out;

   out:
      return;
}



#if defined(SISLNEEDPROTOTYPES)
   static int
    s1789_s9knot(double et[], int ik, int in, double ax1, double ax2,
		 int *jmy, int *jstat)
#else
    static int s1789_s9knot(et, ik, in, ax1, ax2, jmy, jstat)
       double et[];
       int    ik;
       int    in;
       double ax1;
       double ax2;
       int    *jmy;
       int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if there is any knots between two parameter values
*              on a given knot vector. In that case, return the knot index.
*
*
* INPUT      : et       - Knot vector.
*              ik       - Order of spline space.
*              in       - Number of coefficients in spline space.
*              ax1      - First parameter value.
*              ax2      - Second parameter value.
*
*
*
* OUTPUT     : s9knot  - 1 if such a knot exist, 0 otherwise.
*              jmy     - Index of knot. If s9knot == 0, not reliable.
*              jstat   - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, Oslo, Norway. Nov. 1994
*
*********************************************************************
*/
{
   int kstat = 0;         /* Status variable.  */
   int kleft1 = 0;
   int kleft2 = 0;
   int kknot;
   double tref = et[in] - et[ik-1];

   /* Initialize input. */

   *jmy = 0;

   /* Find position of the input parameter values in the given knot vector. */

   s1219(et, ik, in, &kleft1, ax1, &kstat);
   if (kstat < 0) goto error;

   s1219(et, ik, in, &kleft2, ax2, &kstat);
   if (kstat < 0) goto error;

   if (kleft1 != kleft2)
   {
      /* Not the same knot interval. */

      if (ax1 < ax2) (*jmy) = kleft1 + 1;
      else
      {
	 (*jmy) = kleft1 - 1;
	 while (DEQUAL(et[*jmy], et[kleft1])) (*jmy)--;
      }
   }

   if (kleft1 == kleft2 ||
       DEQUAL(et[*jmy]+tref, ax2+tref) ||
       (DEQUAL(et[kleft1]+tref, ax1+tref) && kleft2 == (*jmy) &&
	DEQUAL(et[kleft2]+tref, ax2+tref)))
      kknot = 0;     /* No knot found between the parameter values. */
   else kknot = 1;   /* Knot with index (*jmy) found.               */

   *jstat = 0;
   goto out;

   /* Error in lower level routine. */

   error : *jstat = kstat;
   goto out;

   out:
      return kknot;
}
