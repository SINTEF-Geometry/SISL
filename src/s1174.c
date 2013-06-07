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
 * $Id: s1174.c,v 1.5 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1174

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void s1174_s9corr(double [],double,double,double,double,double,double);
static void s1174_s9dir(double *,double *,double []);
#else
static void s1174_s9corr();
static void s1174_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
s1174(SISLSurf *psurf,double estart[],
     double eend[], double enext[], double gpos[],int *jstat)
#else
void s1174(psurf,estart,eend,enext,gpos,jstat)
     SISLSurf         *psurf;
     double       estart[];
     double       eend[];
     double       enext[];
     double       gpos[];
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on a onedimensional surface.
*              The function finds a local extremum.
*
*
* INPUT      : psurface - Pointer to the surface.
*              estart   - Start values of parameter intervalls.
*              eend     - End value of parameter intervalls.
*              enext    - Parameter start value for iteration.
*
*
* OUTPUT     : gpos    - Parameter values of the found extremum.
*              jstat   - status messages
*                                = 1   : Extremum found.
*                                = 0   : Extremum NOT found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, OCTOBER 1990
* Changed by : Per OEyvind Hvidsten, SINTEF, 11-94
*              Added an initialization of tdist.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int kleft1=0;             /* Variables used in the evaluator.            */
  int kleft2=0;             /* Variables used in the evaluator.            */
  int kder=2;               /* Order of derivatives to be calulated        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist = 0.0;       /* Euclidian norm of derivative vector         */
  double tprev;             /* Previous Euclidian norm of derivative vector*/
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the two parameter directions.      */
  double sval[7];           /* Value ,first and second derivatiev of surf. */
  double *snorm=sval+7;     /* Normal vector of the surface, dummy.        */
  double snext[2];          /* Parameter values                            */
  double tol = (double)10000.0*REL_COMP_RES; /* Singularity tolerance        */
  /* --------------------------------------------------------------------- */

  /* Test input.  */
  if (psurf->idim != 1) goto err106;

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  tdelta[0] = psurf->et1[psurf->in1] - psurf->et1[psurf->ik1 - 1];
  tdelta[1] = psurf->et2[psurf->in2] - psurf->et2[psurf->ik2 - 1];


  /* Initiate variables.  */
  gpos[0] = enext[0];
  gpos[1] = enext[1];

  /* Evaluate 0-2.st derivatives of surface */
  s1421(psurf,kder,gpos,&kleft1,&kleft2,sval,snorm,&kstat);
  if (kstat < 0) goto error;

  /* Get Euclidian norm of derivative vector */
  tprev = sqrt(sval[1]*sval[1] + sval[2]*sval[2]);

  /* Compute the Newton stepdistanse vector. */
  s1174_s9dir(td,td+1,sval);

  if ( (fabs(td[0]/tdelta[0]) <= REL_COMP_RES) &&
      (fabs(td[1]/tdelta[1]) <= REL_COMP_RES))
     goto stop_it;

  /* Adjust if we are not inside the parameter intervall. */
  t1[0] = td[0];
  t1[1] = td[1];
  s1174_s9corr(t1,gpos[0],gpos[1],estart[0],eend[0],estart[1],eend[1]);

  /* Iterate to find the intersection point.  */

  for (knbit = 0; knbit < 50; knbit++)
    {
      /* Evaluate 0-2.st derivatives of surface */

      snext[0] = gpos[0] + t1[0];
      snext[1] = gpos[1] + t1[1];

      s1421(psurf,kder,snext,&kleft1,&kleft2,sval,snorm,&kstat);
      if (kstat < 0) goto error;

      /* Get Euclidian norm of derivative vector */
      tdist = sqrt(sval[1]*sval[1] + sval[2]*sval[2]);

      /* Compute the Newton stepdistanse vector. */
      s1174_s9dir(tdn,tdn+1,sval);

      /* Check if the direction of the step have change. */

      kdir = (s6scpr(td,tdn,2) >= DZERO);     /* 0 if changed. */

      if (tdist <= tprev || kdir)
	{
	  /* Ordinary converging. */

          gpos[0] += t1[0];
          gpos[1] += t1[1];

          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];

	  /* Adjust if we are not inside the parameter intervall. */
	  s1174_s9corr(t1,gpos[0],gpos[1],estart[0],eend[0],estart[1],eend[1]);


          if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES))
	    {
	      gpos[0] += t1[0];
	      gpos[1] += t1[1];

	      break;
	    }

          tprev = tdist;
	}

      else
	{
	  /* Not converging, half step length try again. */

          t1[0] /= (double)2;
          t1[1] /= (double)2;
	  /*         knbit--;  */
	}
    }

  /* Iteration stopped, test if point is extremum */

  stop_it:

  if (tdist <= tol)
    *jstat = 1;
  else
    *jstat = 0;


  /* Test if the iteration is close to a knot */
  if (fabs(gpos[0] - psurf->et1[kleft1])/tdelta[0] < tol)
    gpos[0] = psurf->et1[kleft1];
  else if (fabs(gpos[0] - psurf->et1[kleft1+1])/tdelta[0] < tol)
    gpos[0] = psurf->et1[kleft1+1];

  if (fabs(gpos[1] - psurf->et2[kleft2])/tdelta[1] < tol)
    gpos[1] = psurf->et2[kleft2];
  else if (fabs(gpos[1] - psurf->et2[kleft2+1])/tdelta[1] < tol)
    gpos[1] = psurf->et2[kleft2+1];

  /* Iteration completed.  */
  goto out;

 /* --------------------------------------------------------------------- */
  /* Error in input. Dimension not equal to 1 */
 err106: *jstat = -106;
  s6err("s1174",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("s1174",*jstat,kpos);
  goto out;

 out:;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1174_s9corr(double gd[], double acoef1,double acoef2,double astart1,
		   double aend1,double astart2, double aend2)
#else
static void s1174_s9corr(gd,acoef1,acoef2,astart1,aend1,astart2,aend2)
     double gd[];
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
*              estart2 - The lower boorder in second direction.
*              eend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gd    - Old and new step value.
*
*
* METHOD     : We are cutting a line inside a rectangle.
*	       In this case we always know that the startpoint of
*	       the line is inside the rectangel, and we may therfor
*	       use a simple kind of clipping.
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
  if (acoef1 + gd[0] < astart1)  gd[0] = astart1 - acoef1;
  else if (acoef1 + gd[0] > aend1) gd[0] = aend1 - acoef1;

  if (acoef2 + gd[1] < astart2)  gd[1] = astart2 - acoef2;
  else if (acoef2 + gd[1] > aend2) gd[1] = aend2 - acoef2;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1174_s9dir(double *cdiff1, double *cdiff2,double evals[])
#else
static void s1174_s9dir(cdiff1,cdiff2,evals)
     double *cdiff1;
     double *cdiff2;
     double evals[];
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the step according to a Newton scheme
*              for finding an extremal to a one dimensional surface.
*
*
* INPUT      : evals - Value and derivatives in point on surface.
*
*
* OUTPUT     : cdiff1  - Parameter increment in first direction.
*              cdiff2  - Parameter increment in second direction.
*
*
*
* METHOD     : This is a one dimensional case. We want to find x,y such that
*
*                    x,y : Sx(x,y) = 0
*                          Sy(x,y) = 0
*
*              Using Taylor we get:
*
*                    x,y : Sx+DxSxx+DySxy =0
*                          Sy+DySyy+DxSxy =0
*
*                    x,y: SxxDx + SxyDy = - Sx
*                         SyyDy + SxyDx = - Sy
*
*	       The solution of this matrix equation is the
*	       following function.
*              No effort is done if the matrix is singular,
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, OCTOBER 1990
*
*********************************************************************
*/
{
  double tdiv;		      /* Determinant                               */
  double ta11,ta12,ta21,ta22; /* The matrix                  		   */
  double tmax;                /* The largest value in matrix               */
  double tb1,tb2;             /* The right hand side.                      */
  double tderx,tderxx;        /* Derivatives                               */
  double tdery,tderyy;
  double tderxy;
  double tdeltax,tdeltay;   /* Locals for the step value to be determined. */
  /* --------------------------------------------------------------------- */

  /* Init */
  tderx  = evals[1];
  tdery  = evals[2];
  tderxx = evals[3];
  tderxy = evals[4];
  tderyy = evals[5];
  tdeltax = DZERO;
  tdeltay = DZERO;
  *cdiff1  = DZERO;
  *cdiff2  = DZERO;


  /* Building the matrix. */

  ta11 = tderxx;
  ta12 = tderxy;
  ta21 = tderxy;
  ta22 = tderyy;
  tb1  = -tderx;
  tb2  = -tdery;

  tmax = max(fabs(ta11),max(fabs(ta12),max(fabs(ta21),fabs(ta22))));

  if (DEQUAL(tb1+tmax,tmax) && DEQUAL(tb2+tmax,tmax))
    {
      /* Finished, we have found a max. */
    }
  else
    {
      tdiv    = ta11*ta22 - ta21*ta12;
      if (fabs(tdiv) > MAX(tmax*REL_COMP_RES,REL_COMP_RES))
	{
	  /* The matrix is ok, solve the system using Cramers rule. */
	  tdeltax = tb1*ta22 - tb2*ta12;
	  tdeltay = ta11*tb2 - ta21*tb1;
	  tdeltax /= tdiv;
	  tdeltay /= tdiv;
	}
      else if (max (fabs(ta11),fabs(ta22)) > REL_COMP_RES)
	{
	   if (fabs(ta11) > fabs(ta22))
	     tdeltax = tb1/ta11;
	   else
	     tdeltay = tb2/ta22;
	}

    }

  *cdiff1  = tdeltax;
  *cdiff2  = tdeltay;

}
