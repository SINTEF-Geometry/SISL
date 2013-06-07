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
 * $Id: s1252.c,v 1.3 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1252

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined (SISLNEEDPROTOTYPES)
static void s1252_s6corr(double *,double,double [],int,int,int *,int *);
static void s1252_s6dir(double *,double,double [],double,double);
#else
static void s1252_s6corr();
static void s1252_s6dir();
#endif

#if defined (SISLNEEDPROTOTYPES)
void

     s1252(SISLCurve *pcurve,double aepsge,double astart,double *cpos,int *jstat)
#else
void s1252(pcurve,aepsge,astart,cpos,jstat)
     SISLCurve  *pcurve;
     double aepsge;
     double astart;
     double *cpos;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration to a local maximum on a function
*              in one variable.
*
* INPUT      : pcurve  - Pointer to the first curve in the intersection.
*              aepsge  - Geometry resolution.
*              astart  - Start value of the first curve to the iteration.
*
*
*
* OUTPUT     : cpos    - Parameter value of of first curve in intersection
*                        point.
*              jstat   - status messages
*                                = 2   : Divergence or approximative
*                                        intersection found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in one parameter direction.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221 - Evaluate expression of curve in given
*                         parameter values.
*
* WRITTEN BY : Tor Dokken, SI, Mars 1989
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable.                          */
  int kpos = 0;         /* Position of error.                              */
  int kleft=0;          /* Variables used in the evaluator.                */
  int kder=3;           /* Order of derivatives to be calulated            */
  int kdim;             /* Dimension of space the curves lie in            */
  int knbit;            /* Number of iterations                            */
  int kn,kk;            /* Number of vertices and order                    */
  int kdir=1;           /* Direction of derivative to be calculated        */
  double tstart,tend;   /* Ends of parameter interval of first curve.      */
  double tdelta;        /* Parameter interval of the curves.               */
  double tdist=DZERO;   /* Distance between position and origo.            */
  double td;        	/* Distances between old and new parameter value   */
  double tnext;         /* Parameter-value of expression in first curve.   */
  double tprev;         /* Previous difference between the curves.         */
  double sval[4];       /* Value ,first and second derivative on function  */
  double *st;           /* Knot vector                                     */
  double ref;           /* Refferance value for equality test.             */

  /* Test input.  */

  if (pcurve->idim != 1) goto err106;

  kdim = pcurve -> idim;

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  st = pcurve->et;
  kn = pcurve->in;
  kk = pcurve->ik;

  tstart = *(pcurve->et + pcurve->ik - 1);
  tend   = *(pcurve->et + pcurve->in);
  tdelta = tend - tstart;
  if (tdelta == DZERO) tdelta = fabs(tend);
  if (tdelta == DZERO) tdelta = (double)1.0;

  /* Initiate variables.  */

  tnext = astart;

  /* Evaluate 0-1.st derivatives of function */

  s1221(pcurve,kder,tnext,&kleft,sval,&kstat);
  if (kstat < 0) goto error;

  tprev = sval[0];

  /* Evaluate step */

  s1252_s6dir(&td,tnext,sval,tstart,tend);

  /* Correct if we not are inside the parameter intervall. */

  s1252_s6corr(&td,tnext,st,kn,kk,&kleft,&kdir);

  /* Iterate to find the intersection point.  */

  for (knbit = 0; knbit < 20; knbit++)
    {

      /* If the tnext is a break point test if it is a local maximum */

      if (kdir == -2 || kdir == 2)
	{
	  double tder1,tder2;
	  /* Break point, test if local maximum */

	  s1221(pcurve,kder,tnext,&kleft,sval,&kstat);
	  if (kstat < 0) goto error;
	  tder2 = sval[1];

	  s1227(pcurve,kder,tnext,&kleft,sval,&kstat);
	  if (kstat < 0) goto error;
	  tder1 = sval[1];

	  /*    Test if top point */

	  if (tder1>=DZERO && tder2<=DZERO) break;

	  /*    Not a top point */
	}


      /* Evaluate 0-1.st derivatives of both curves, dependent of the
	 sign of td we calculate derivatives from the right or the left */

      if (kdir>=1)
	{
	  s1221(pcurve,kder,tnext+td,&kleft,sval,&kstat);
	  if (kstat < 0) goto error;
	}
      else
	{
	  s1227(pcurve,kder,tnext+td,&kleft,sval,&kstat);
	  if (kstat < 0) goto error;
	}

        tdist = sval[0];
        if (fabs(tdist) < (double)1.0) ref = (double)2.0;
	else                           ref = DZERO;

        if (tdist >= tprev || DEQUAL(ref+tdist,ref+tprev))
	{
	   tnext += td;

	   /* Evaluate step */
	   s1252_s6dir(&td,tnext,sval,tstart,tend);
	   s1252_s6corr(&td,tnext,st,kn,kk,&kleft,&kdir);

	   if (fabs(td/tdelta) <= REL_COMP_RES) break;

	   tprev = tdist;

	}

      /* Not converging, correct and try again. */

      else
	{

	  td /= (double)2;
	  if (fabs(td/tdelta) <= REL_COMP_RES) break;
	}


    }


  /* Iteration stopped, test if point founds found is within resolution */

  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;

  /*ujk,july 89:*/
  /* Test if the iteration is close to a knot */
  if (DEQUAL(tnext,pcurve->et[kleft]))
    *cpos = pcurve->et[kleft];
  else if (DEQUAL(tnext,pcurve->et[kleft+1]))
    *cpos = pcurve->et[kleft+1];
  else
    *cpos = tnext;

  /* Iteration completed.  */

  goto out;


  /* Error in input. Conflicting dimensions.  */

 err106: *jstat = -106;
  s6err("S1252",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("S1252",*jstat,kpos);
  goto out;

 out:;
}

#if defined (SISLNEEDPROTOTYPES)
static void
   s1252_s6corr(double *gdn,double acoef,double et[],
		int in,int ik,int *ileft,int *jdir)
#else
static void s1252_s6corr(gdn,acoef,et,in,ik,ileft,jdir)
     double *gdn;
     double acoef;
     double et[];
     int    in;
     int    ik;
     int    *ileft;
     int    *jdir;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Adjust the step to not cross knot values or out
*              of the curve
*
*
* INPUT      : acoef   - Current parameter value
*              st      - knots
*              in      - number of vertices
*              ik      - polynomial order
*
* INPUT/OUTPUT :
*              ileft - Pointer to the interval in the knot vector
*                       where ax is located, check the relations above.
*
*
* OUTPUT     : gdn     - Old and new step value.
*              jdir    - Direction of derivative to be calculated
*                         -2 Negative direction acoef at break point
*                         -1 Negative direction
*                         +1 Positive direction
*                         +2 Positive direction acoef at break point
*
* METHOD     : We are making the step keep inside the parameter interval.
*
* REFERENCES :
*
*-
* CALLS      :
*
*
* WRITTEN BY : Tor Dokken, SI, Mars 1989
*
*********************************************************************
*/
{
  int kmult,kstat;

  /* Make sure the point is inside the interval */

  *gdn = MAX(et[ik-1]-acoef,*gdn);
  *gdn = MIN(et[in]  -acoef,*gdn);

  if (acoef+*gdn<et[*ileft] && acoef>et[*ileft])
    {
      *gdn = MAX(et[*ileft]-acoef,*gdn);
    }

  else if(acoef<et[*ileft+1] && acoef+*gdn>et[*ileft+1])
    {
      /*  We cross a knot value */

      *gdn = MIN(et[*ileft+1]-acoef,*gdn);
    }

  /* Make sure that we calculate the left or right handed derivatives */

  if (*gdn>=0)
    {
      *jdir = 1;
    }
  else
    {
      *jdir = -1;
    }

  kmult = s6knotmult(et,ik,in,ileft,acoef,&kstat);

  if (acoef==et[*ileft])
    {

      if(kmult>ik-2)
        {
	  if (*jdir == -1)
            {
	      *jdir = -2;
            }
	  else
            {
	      *jdir =  2;
            }
        }
    }
}

#if defined (SISLNEEDPROTOTYPES)
static void
  s1252_s6dir(double *cdiff,double acoef,double eval[],double astart,
	      double aend)
#else
static void s1252_s6dir(cdiff,acoef,eval,astart,aend)
     double *cdiff;
     double acoef;
     double eval[];
     double astart;
     double aend;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the next iteration step
*
* INPUT      : eval    - Value and derivative
*              astart  - The lower interval end
*              aend    - The higher interval end
*
*
* OUTPUT     : cdiff   - Iteration step
*-
*
* WRITTEN BY : Tor Dokken, SI, Mars 1989
*
*********************************************************************
*/
{
  double t1,t2,t3,t4,t5,t6;   /* Constants in equation.                    */
  double tmax;                /* Max values in equation.                   */
  double ttol=(double)1e-10;  /* Relative tolerance in equation.           */

  /* Dummy statements to avoid warning. */
  t1=acoef;
  t2=astart;
  t3=aend;


  t1 =  eval[1];
  t2 =  eval[2];
  t3 =  eval[3]/(double)2.0;

  tmax  = max(fabs(t1),fabs(t2));
  tmax  = max(fabs(t3),tmax);

  if (DEQUAL(tmax,DZERO))                    *cdiff = DZERO;
  else if (fabs(t3)/tmax < ttol) /* The second degree part is degenerated. */
	{
          if (fabs(t2) == DZERO )      *cdiff = DZERO;
	  else                        *cdiff = (-t1/t2);
	}
  else
	{
          /* An ordinary second degree equation.    */
	   t4 = t2*t2 - (double)4*t3*t1;
	   if (t4 < DZERO)
	    {
	      /* Use linear equation. */
	      if (fabs(t2) == DZERO )      *cdiff = DZERO;
              else                        *cdiff = (-t1/t2);
      	    }

           else
	    {
	       t6 = sqrt(t4);
	       t5 = (-t2 + t6)/((double)2*t3);
	       t6 = (-t2 - t6)/((double)2*t3);
	       t4 = min(fabs(t5),fabs(t6));

               /* We have two solutions and we want to use the one
	          with the one with smallest value. */

               if (t4 == DZERO)
                {
	          /* Use linear equation. */
	          if (fabs(t2) == DZERO )      *cdiff = DZERO;
                  else                        *cdiff = (-t1/t2);
	        }
               else if (fabs(t5) <= fabs(t6))  *cdiff = t5;
               else                            *cdiff = t6;
             }
	}
}
