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
 * $Id: s1015.c,v 1.4 2001-03-19 15:58:40 afr Exp $
 *
 */
#define S1015

#include "sislP.h"

/*
*  Forward declarations.
*  ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static
void
c_c_f_s9corr(double [],double,double,double,double,double,double);
static
void
c_c_f_s9dir(double *,double *,double *,double [],double [],double [],int);
#else
static void c_c_f_s9corr();
static void c_c_f_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
      s1015(SISLCurve *pc1,SISLCurve *pc2,double aepsge,
	   double eps1[],double eps2[],double aradius,double *parpt1,
	   double *parpt2,double center[],int *jstat)
#else
void s1015(pc1,pc2,aepsge,eps1,eps2,aradius,parpt1,parpt2,
           center,jstat)
     SISLCurve  *pc1;
     SISLCurve  *pc2;
     double aepsge;
     double eps1[];
     double eps2[];
     double aradius;
     double *parpt1;
     double *parpt2;
     double center[];
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To iterate to the start and end point of a fillet between
*              two 2D curves. The center of the (circular) fillet is also
*              calculated.
*
* INPUT      : pc1    - The first 2D input curve.
*              pc2    - The second 2D input curve.
*              aepsge - Geometry resolution.
*              eps1   - 2D point telling that the fillet should be put on
*                       the side of curve 1 where eps1 is situated.
*              eps2   - 2D point telling that the fillet should be put on the
*                       side of curve 2 where eps2 is situated.
*              aradius - The radius to be used on the fillet.
*
* INPUT/OUTPUT :
*              parpt1 - Parameter value of the point on curve 1 where the
*                       fillet starts. Input is a guess value for the iteration.
*              parpt2 - Parameter value of the point on curve 2 where the
*                       fillet ends. Input is a guess value for the iteration.
*
* OUTPUT     : center - 2D center of the (circular) fillet.  Space must be
*                       allocated outside the function.
*              status - Status:
*                              = 1      : converged
*                              = 2      : diverged
*                              < 0      : error
*
* METHOD     :
*
* USE        :
*
* REFERENCES :
*
*-
* CALLS      :
*
*
* WRITTEN BY : Johannes Kaasa, SI, April 1992.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Clearified
*              header and changed the function interface declaration of 'center'
*              - must now be allocated outside the function (was memory overwrite).
*              Initialized 'p1', 'p2' and 'kstat' and fixed memory leak from
*              'p1' and 'p2'.
*
*********************************************************************
*/
{

  int kstat=0;        /* Status variable                                  */
  int kpos=0;         /* Position of error                                */

  SISLPoint *p1=SISL_NULL; /* SISLPoint form of eps1.                          */
  SISLPoint *p2=SISL_NULL; /* SISLPoint form of eps2.                          */
  double tpar1;       /* Test parameter on curve 1.                       */
  double tpar2;       /* Test parameter on curve 2.                       */
  int kder = 1;       /* Number of derivatives.                           */
  int kleft1=0;       /* Pointer to the interval in the knot vector.      */
  int kleft2=0;       /* Pointer to the interval in the knot vector.      */
  double derive[4];   /* Position and derivative on a curve.              */
  double tdot;        /* Dot product of two vectors.                      */
  double rad1;        /* Offset radius of curve 1.                        */
  double rad2;        /* Offset radius of curve 2.                        */

  int kdim = 2;       /* Dimension of the geometry.                       */
  int knbit;          /* Number of iterations                             */
  int kdir;           /* Changing direction.                              */
  double tdelta1;     /* Parameter interval of the curves.                */
  double tdelta2;     /* Parameter interval of the curves.                */
  double tdist;       /* Distance between position and origo.             */
  double td[2];       /* Distances between old and new parameter
			 value in the two parameter directions.           */
  double t1[2];       /* Distances between old and new parameter
			 value in the two parameter directions.           */
  double tdn[2];      /* Distances between old and new parameter
			 value in the two parameter directions.           */
  double tprev;       /* Previous difference between the curves.          */
  double *sval1=SISL_NULL; /* Value ,first and second derivatie on curve 1     */
  double *sval2;      /* Value ,first and second derivatie on curve 2     */
  double *sdiff;      /* Difference between the curves                    */
  double start1;      /* Start parameter of curve 1.                      */
  double start2;      /* Start parameter of curve 2.                      */
  double end1;        /* End parameter of curve 1.                        */
  double end2;        /* End parameter of curve 2.                        */

  /* Check dimensions */

  if (pc1->idim != 2 || pc2->idim != 2) goto err105;

  start1 = pc1->et[pc1->ik - 1];
  start2 = pc2->et[pc2->ik - 1];
  end1 = pc1->et[pc1->in];
  end2 = pc2->et[pc2->in];

  /* Calculate the offset radius of curve 1. */

  if (!(p1 = newPoint(eps1, 2, 0))) goto err101;
  s1771(p1, pc1, REL_COMP_RES, start1, end1, *parpt1,
	&tpar1, &kstat);
  if (kstat < 0) goto error;

  s1221(pc1, kder, tpar1, &kleft1, derive, &kstat);
  if (kstat < 0) goto error;

  tdot = (eps1[0] - derive[0])*(-derive[3])
                              + (eps1[1] - derive[1])*derive[2];
  if (tdot > (double)0.)
    rad1 = aradius;
  else
    rad1 = -aradius;

  /* Calculate the offset radius of curve 2. */

  if (!(p2 = newPoint(eps2, 2, 0))) goto err101;
  s1771(p2, pc2, REL_COMP_RES, start2, end2, *parpt2,
	&tpar2, &kstat);
  if (kstat < 0) goto error;

  s1221(pc2, kder, tpar2, &kleft2, derive, &kstat);
  if (kstat < 0) goto error;

  tdot = (eps2[0] - derive[0])*(-derive[3])
                              + (eps2[1] - derive[1])*derive[2];
  if (tdot > (double)0.)
    rad2 = aradius;
  else
    rad2 = -aradius;

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  tdelta1 = end1 - start1;
  tdelta2 = end2 - start2;

  /* Allocate local used memory */

  sval1 = newarray((2*kder+5)*kdim,double);
  if (sval1 == SISL_NULL) goto err101;

  sval2 = sval1 + (kder+2)*kdim;
  sdiff = sval2 + (kder+2)*kdim;

  /* Initiate variables.  */

  tprev = (double)HUGE;

  /* Evaluate 0-1.st derivatives of both curves */

  ev_cv_off(pc1, kder, *parpt1, &kleft1, rad1, sval1, &kstat);
  if (kstat < 0) goto error;

  ev_cv_off(pc2, kder, *parpt2, &kleft2, rad2, sval2, &kstat);
  if (kstat < 0) goto error;

  /* Compute the distanse vector and value and the new step. */

  c_c_f_s9dir(&tdist,td,td+1,sdiff,sval1,sval2,kdim);

  /* Correct if we are not inside the parameter intervall. */

  t1[0] = td[0];
  t1[1] = td[1];
  c_c_f_s9corr(t1,*parpt1,*parpt2,start1,end1,start2,end2);

  /* Iterate to find the intersection point.  */

  for (knbit = 0; knbit < 50; knbit++)
    {
      /* Evaluate 0-1.st derivatives of both curves */

       ev_cv_off(pc1, kder, *parpt1+t1[0], &kleft1, rad1, sval1, &kstat);
      if (kstat < 0) goto error;

      ev_cv_off(pc2, kder, *parpt2+t1[1], &kleft2, rad2, sval2, &kstat);
      if (kstat < 0) goto error;

      /* Compute the distanse vector and value and the new step. */

      c_c_f_s9dir(&tdist,tdn,tdn+1,sdiff,sval1,sval2,kdim);

      /* Check if the direction of the step have change. */

      kdir = (s6scpr(td,tdn,2) >= DZERO);     /* 0 if changed. */

      /* Ordinary converging. */

      if (tdist < tprev*(double)0.9 || kdir)
	{
	  *parpt1 += t1[0];
	  *parpt2 += t1[1];

          td[0] = tdn[0];
          td[1] = tdn[1];

	  /* Correct if we are not inside the parameter intervall. */

	  t1[0] = td[0];
	  t1[1] = td[1];
	  c_c_f_s9corr(t1,*parpt1,*parpt2,start1,end1,start2,end2);

          if ( (fabs(t1[0]/tdelta1) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta2) <= REL_COMP_RES) ) break;

          tprev = tdist;
	}

      /* Not converging, corrigate and try again. */

      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          knbit--;
	}
    }

  /* Calculate the centerpoint. */

  center[0] = (sval1[0] + sval2[0])/2.;
  center[1] = (sval1[1] + sval2[1])/2.;

  /* Iteration stopped, test if point founds found is within resolution */

  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;

  /* Iteration completed.  */


  goto out;

  /* Error in allocation */

 err101: *jstat = -101;
  s6err("s1015",*jstat,kpos);
  goto out;

  /* Dimension nmot equal to 2 */

 err105:
  *jstat = -105;
  s6err("s1015",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1015",*jstat,kpos);
  goto out;

 out:
  if (sval1 != SISL_NULL) freearray(sval1);
  if ( p1 != SISL_NULL )  freePoint(p1);
  if ( p2 != SISL_NULL )  freePoint(p2);
}

#if defined(SISLNEEDPROTOTYPES)
static
void
   c_c_f_s9corr(double gdn[],double acoef1,double acoef2,
		   double astart1,double aend1,double astart2,double aend2)
#else
static void c_c_f_s9corr(gdn,acoef1,acoef2,astart1,aend1,astart2,aend2)
     double gdn[];
     double acoef1;
     double acoef2;
     double astart1;
     double aend1;
     double astart2;
     double aend2;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To be sure that we are inside the boorder of the
*              parameter plan. If we are outside clipping is used
*	       to corrigate the step value.
*
*
* INPUT      : acoef1  - Coeffisient in the first direction.
*              acoef2  - Coeffisient in the second direction.
*              astart1 - The lower boorder in first direction.
*              aend1   - The higher boorder in first direction.
*              astart2 - The lower boorder in second direction.
*              aend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gdn   - Old and new step value.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/
{
  if (acoef1 + gdn[0] < astart1)
    gdn[0] = astart1 - acoef1;
  else if (acoef1 + gdn[0] > aend1)
    gdn[0] = aend1 - acoef1;

  if (acoef2 + gdn[1] < astart2)
    gdn[1] = astart2 - acoef2;
  else if (acoef2 + gdn[1] > aend2)
    gdn[1] = aend2 - acoef2;
}

#if defined(SISLNEEDPROTOTYPES)
static
void
   c_c_f_s9dir(double *cdist,double *cdiff1,double *cdiff2,
		  double gdiff[],double eval1[],double eval2[],int idim)
#else
static void c_c_f_s9dir(cdist,cdiff1,cdiff2,gdiff,eval1,eval2,idim)
     double *cdist;
     double *cdiff1;
     double *cdiff2;
     double eval1[];
     double eval2[];
     double gdiff[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the distance vector and value beetween
*	       a points on the first curve and a point on the second
*	       curve. And to compute a next step on both curves.
*	       This is equivalent to the nearest way to the
*	       parameter plan in the tangent plan from a point in the
*	       distance surface between two curves.
*
*
* INPUT      : eval1 - Value and derevative in first parameterdirection.
*              eval2 - Value and derevative in second parameterdirection.
*	       idim  - Dimension of space the curves lie in.
*
*
* OUTPUT     : gdiff   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to each other.
*	       The differens vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the two
*	       points on the curve is "diff" and the two derevativ vectors
*	       are "der1" and "der2", and the two wanted parameter distance
*	       are "dt1" and "dt2", then we get the following system of
*	       equations:
*		 <dt1*der1+dist-dt2*der2,der2> = 0
*		 <dt1*der1+dist-dt2*der2,der1> = 0
*	       This is futher:
*
*		 | -<der1,der2>   <der2,der2> |  | dt1 |   | <dist,der2> |
*		 |                            |  |     | = |             |
*		 | -<der1,der1>   <der1,der2> |  | dt2 |   | <dist,der1> |
*
*	       The solution of this matrix equation is the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/
{
  int kstat;		           /* Local status variable. */
  register double tdet;		   /* Determinant */
  register double t1,t2,t3,t4,t5;  /* Variables in equation system */

  /* Computing the different vector */

  s6diff(eval1,eval2,idim,gdiff);

  /* Computing the length of the different vector. */

  *cdist = s6length(gdiff,idim,&kstat);

  t1 =  s6scpr(eval1+idim,eval1+idim,idim);
  t2 =  s6scpr(eval1+idim,eval2+idim,idim);
  t3 =  s6scpr(eval2+idim,eval2+idim,idim);
  t4 =  s6scpr(gdiff,eval1+idim,idim);
  t5 =  s6scpr(gdiff,eval2+idim,idim);

  /* Computing the determinant. */

  tdet = t2*t2 - t1*t3;

  if (DEQUAL(tdet,DZERO))
    {
      *cdiff1 = DZERO;
      *cdiff2 = DZERO;
    }
  else
    {
      /* Using Cramer's rule to find the solution of the system. */

      *cdiff1 =  (t4*t3 - t5*t2)/tdet;
      *cdiff2 =  (t2*t4 - t1*t5)/tdet;
    }
}
