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
 * $Id: s1771.c,v 1.5 2005-02-28 09:04:49 afr Exp $
 *
 */
#define S1771

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static
void
   s1771_s9point(SISLCurve *,double [],double [],double [],double,double,
		 int,double *,double *,double,double *,double,int,int *);
/* static
void
 s1771_s9corr(double[],double,double,double); */
static
double
s1771_s9del(double *,double *,double *,int);
#else
static void s1771_s9point();
/* static void s1771_s9corr(); */
static double s1771_s9del();
#endif


#if defined(SISLNEEDPROTOTYPES)
void
s1771(SISLPoint *ppoint,SISLCurve *pcurve,double aepsge,
	   double astart,double aend,double anext,double *cpos,int *jstat)
#else
void s1771(ppoint,pcurve,aepsge,astart,aend,anext,cpos,jstat)
     SISLPoint  *ppoint;
     SISLCurve  *pcurve;
     double aepsge;
     double astart;
     double aend;
     double anext;
     double *cpos;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameter is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint  - The point in the closest point problem.
*              pcurve  - The curve in the closest point problem.
*              aepsge  - Geometrical resolution.
*              astart  - Curve parameter giving the start of the search
*                        interval.
*              aend    - Curve parameter giving the end of the search
*                        interval.
*              anext   - Curve guess parameter for the closest point
*                        iteration.
*
*
*
* OUTPUT     : cpos    - Resulting curve parameter from the iteration.
*              jstat   - status messages
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/
{
  int kstat = 0;         /* Local status variable.                           */
  int kpos = 0;          /* Position of error.                               */
  int kleft=0;           /* Variables used in the evaluator.                 */
  int kdim;              /* Dimension of space the curves lie in             */
  double tdelta;         /* Parameter interval of the curve.                 */
  double tdist;          /* Distance between position and origo.             */
  double td;             /* Distances between old and new parameter value in
			    the two parameter directions.                    */
  double tprev;          /* Previous difference between the curves.          */
  double *sval=SISL_NULL;     /* Value ,first and second derivatie on curve 1    */
  double *sdiff;         /* Difference between the curves                    */
  int quick = (*jstat);  /* Indicates if the exactness requirement is
                            relaxed.                                         */
  int max_it = 20;       /* Maximum number of iterations.                    */

  if (quick) max_it = 10;

  /* Test input.  */

  if (ppoint->idim != pcurve->idim) goto err106;

  kdim = pcurve -> idim;

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  tdelta = pcurve->et[pcurve->in] - pcurve->et[pcurve->ik - 1];

  /* Allocate local used memory */

  sval = newarray(4*kdim,double);
  if (sval == SISL_NULL) goto err101;

  sdiff = sval + 3*kdim;

  /* Initiate variable.  */

  tprev = (double)HUGE;

  /* Evaluate 0-2.st derivatives of the curve. */

  s1221(pcurve,2,anext,&kleft,sval,&kstat);
  if (kstat < 0) goto error;

  s6diff(ppoint->ecoef,sval,kdim,sdiff);

  tdist = s6length(sdiff,kdim,&kstat);

  td = s1771_s9del(sdiff,sval+kdim,sval+kdim+kdim,kdim);

  /* Correct if we are not inside the parameter intervall. */

  if (anext + td < astart) td = astart - anext;
  else if (anext + td > aend) td = aend - anext;

  s1771_s9point(pcurve,ppoint->ecoef,sval,sdiff,astart,aend,max_it,&anext,
	  &td,tdelta,&tdist,tprev,kleft,&kstat);
  if (kstat < 0) goto error;

  /* Iteration stopped, test if point found is within resolution */

  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;

  *cpos = anext;

  /* Iteration completed.  */

  goto out;

  /* Error in allocation */

 err101: *jstat = -101;
  s6err("s1771",*jstat,kpos);
  goto out;

  /* Error in input. Conflicting dimensions.  */

 err106: *jstat = -106;
  s6err("s1771",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1771",*jstat,kpos);
  goto out;

 out:    if (sval != SISL_NULL) freearray(sval);
}

#if defined(SISLNEEDPROTOTYPES)
static
void
s1771_s9point(SISLCurve *pcurve,double eval1[],double eval2[],double ediff[],
	      double astart,double aend,int max_it,double *cnext,double *ad,
	      double adel,double *cdist,double aprev,int ileft,int *jstat)
#else
static void s1771_s9point(pcurve,eval1,eval2,ediff,astart,aend,max_it,cnext,ad,adel,
			  cdist,aprev,ileft,jstat)
     SISLCurve  *pcurve;
     double eval1[];
     double eval2[];
     double ediff[];
     double astart;
     double aend;
     int    max_it;
     double *cnext;
     double *ad;
     double adel;
     double *cdist;
     double aprev;
     int    ileft;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a point and a curve to find a closest point or an
*              intersection point.
*
*
* INPUT      : pcurve  - Pointer to the curve.
*              eval1[] - Array containing the coefisients to the point.
*              eval2[] - Array containing the coefisients to a start
*                        iteration point on the curve.
*              ediff[] - The vector between the point and the start
*                        point on the curve.
*              astart  - Start parameter value on the curve.
*              aend    - End parameter value on the curve.
*              max_it  - Maximum number of iterations.
*              ad      - A first computed parametric step.
*              adel    - The differrence between start and end
*                        parameter value.
*              aprev    -A prevous distance between the curve and
*                        the point.
*              ileft   - An estimate of the knot number whice is
*                        closest to the point on the curve.
*
*
* OUTPUT     : cnext   - Parameter value of curve in the closest point.
*              cdist    -The current distance between the curve and
*                        the point.
*              jstat   - status messages
*                                >= 0   : OK.
*                                <  0   : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, May 1989
* MODIFIED BY : Vibeke Skytt, SINTEF SI, 10.93. Stop earlier when divergence.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int kdim = pcurve->idim;  /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdiv = 0;             /* Counts number of diverging steps.           */
  int kdiv2 = 0;            /* Counts number of almost divergence.         */

  /* Correct if we are not inside the parameter intervall. */

  if (*cnext + *ad < astart)    *ad = astart - *cnext;
  else if (*cnext + *ad > aend) *ad = aend - *cnext;

  for (knbit=0;knbit<max_it;knbit++)
    {
      /* Evaluate 0-2.st derivatives of the curve. */

      s1221(pcurve,2,*cnext + *ad,&ileft,eval2,&kstat);
      if (kstat < 0) goto error;

      s6diff(eval1,eval2,kdim,ediff);

      *cdist = s6length(ediff,kdim,&kstat);

      if (*cdist -aprev  <= REL_COMP_RES)
	{
	   if (kdiv2 > 4) break;
	   if (*cdist -aprev >= DZERO) kdiv2++;

	   kdiv = 0;
	  aprev = *cdist;
	  *cnext += *ad;

	  *ad = s1771_s9del(ediff,eval2+kdim,eval2+kdim+kdim,kdim);

	  /* Correct if we are not inside the parameter intervall. */

	  if (*cnext + *ad < astart)    *ad = astart - *cnext;
	  else if (*cnext + *ad > aend) *ad = aend - *cnext;
	}
      else
	{
	   kdiv++;
	   if (kdiv > 3) break;
	  (*ad) /= (double)2;
	  /** knbit--; */
	}
      if (fabs((*ad)/MAX(fabs(*cnext),adel)) <= REL_COMP_RES) break;
    }

  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1771_s9point",*jstat,0);
  goto out;

 out:   ;
}
/*

#if defined(SISLNEEDPROTOTYPES)
static
void
s1771_s9corr(double gdn[],double acoef,double astart,double aend)
#else
static void s1771_s9corr(gdn,acoef,astart,aend)
     double gdn[];
     double acoef;
     double astart;
     double aend;
#endif
*/
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
* INPUT      : acoef  - Coeffisient in the first direction.
*              astart - The lower boorder in first direction.
*              aend   - The higher boorder in first direction.
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
/*
{
  if (acoef + gdn[0] < astart)
    gdn[0] = astart - acoef;
  else if (acoef + gdn[0] > aend)
    gdn[0] = aend - acoef;
}
*/

#if defined(SISLNEEDPROTOTYPES)
static
double
s1771_s9del(double *eco,double *eco1,double *eco2,int idim)
#else
static double s1771_s9del(eco,eco1,eco2,idim)
     double *eco;
     double *eco1;
     double *eco2;
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the distance on the parameter line to a point
*            on the curve on which the tangent is orthogonal to the
*            difference vector from this point on the curve to the
*            point in the space.
*
*
* INPUT      : eco   - The differens vector beetween the point and the
*                      current posision on the curve.
*              eco1  - The first derevative vector in the  current posision
*                      on the curve.
*              eco2  - The second derevative vector in the  current posision
*                      on the curve.
*              idim  - Dimension of space the vectors lie in.
*
*
* OUTPUT     : s1771_s9del - The computed parameter distance.
*
*
* METHOD     : We have to find the parameter distance "dt" from
*              the equation:
*                <ecoef-dt*ecoef1,ecoef1+dt*ecoef2> = 0.
*              This may be written:
*                  <ecoef,ecoef1> + <ecoef,ecoef2>*dt
*                - <ecoef1,ecoef1>*dt + <ecoef1,ecoef2>*dt*dt = 0
*              The following function is the solution of this second
*              degree equation. We are always using the solution
*              with the smallest absolute value.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Mar 1989
*
*********************************************************************
*/
{
  double t1,t2,t3,t4,t5,t6;   /* Constants in equation.                    */
  double tmax,tmax1;          /* Max values in equation.                   */
  double ttol=(double)1e-10;  /* Relative tolerance in equation.           */

  t1 =  s6scpr(eco,eco1,idim);
  t3 =  s6scpr(eco1,eco1,idim);
  t2 =  t3 - s6scpr(eco,eco2,idim);
  t4 =  -(double)2 * s6scpr(eco1,eco2,idim);

  tmax  = max(fabs(t1),fabs(t2));
  tmax1 = max(fabs(t3),fabs(t4));
  tmax  = max(tmax1,tmax);

  if (DEQUAL(tmax,DZERO))                    return DZERO;

  else if (fabs(t4)/tmax < ttol) /* The second degree part is degenerated. */
    {
      if (fabs(t2)/tmax < ttol)
	{
          if (fabs(t3)/tmax < ttol)        return DZERO;
          else                             return (t1/t3);
	}
      else                                  return (t1/t2);
    }
  else  /* An ordinary second degree equation.    */
    {
      t5 = t2*t2 - (double)2*t4*t1;
      if (t5 < DZERO)                       return (t1/t3);
      else
	{
          t6 = sqrt(t5);
          t5 = (t2 + t6)/t4;
          t6 = (t2 - t6)/t4;
	  t1 *= t3;


          /* We have two solutions and we want to use the one
	     with the same sign as we get while using an other
	     metode t1/t3. If both solutions have the same
	     sign we use the one with smallest value. */

          if (t1 < DZERO)
	    {
	      if (t5 <= DZERO && t6 <= DZERO)
		{
		  if (t5 > t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 <= DZERO)        return t5;
	      else if (t6 <= DZERO)        return t6;
              else                         return min(t5,t6);
	    }
	  else if (t1 > DZERO)
	    {
	      if (t5 >= DZERO && t6 >= DZERO)
		{
		  if (t5 < t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 >= DZERO)        return t5;
	      else if (t6 >= DZERO)        return t6;
              else                         return max(t5,t6);
	    }
	  else                             return min(fabs(t5),fabs(t6));
	}
    }
}
