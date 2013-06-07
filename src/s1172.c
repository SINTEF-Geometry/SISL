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
 * $Id: s1172.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1172

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void s1172_s9corr(double *,double,double,double);
static void s1172_s9dir(double *,double []);
#else
static void s1172_s9corr();
static void s1172_s9dir();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1172(SISLCurve *pcurve,double astart,
     double aend, double anext, double *cpos,int *jstat)
#else
void s1172(pcurve,astart,aend,anext,cpos,jstat)
     SISLCurve        *pcurve;
     double       astart;
     double       aend;
     double       anext;
     double       *cpos;
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on a onedimensional curve.
*              The function finds a local extremum.
*
*
* INPUT      : pcurve   - Pointer to the curve.
*              astart   - Start values of parameter intervall.
*              aend     - End value of parameter intervall.
*              anext    - Parameter start value for iteration.
*
*
* OUTPUT     : cpos    - Parameter value of the found extremum.
*              jstat   - status messages  
*                                = 1   : Extremum found.
*                                = 0   : Extremum NOT found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in one parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, OCTOBER 1990
* CORRECTED BY:  Ulf J. Krystad, SI, AUGUST 1991
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int kleft=0;              /* Variables used in the evaluator.            */
  int kder=3;               /* Order of derivatives to be calulated        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  double tdelta;            /* Parameter intervals of the Curve.        */
  double tdist;             /* Euclidian norm of derivative vector         */
  double tprev;             /* Previous Euclidian norm of derivative vector*/
  double td,t1,tdn;         /* Distances between old and new parameter
			       value in the two parameter directions.      */
  double sval[4];           /* Value ,first and second derivatiev of Curve.*/ 
  double tnext;             /* Parameter values                            */
  double tol = (double)1000.0*REL_COMP_RES; /* Singularity tolerance      */
  /* --------------------------------------------------------------------- */
  
  /* Test input.  */
  if (pcurve->idim != 1) goto err106;
  
  /* Fetch endpoints and the interval of parameter interval of curves.  */
  
  tdelta = pcurve->et[pcurve->in] - pcurve->et[pcurve->ik - 1];
  
  /* Evaluate 0-2.st derivatives of curve */
  s1221(pcurve,kder,anext,&kleft,sval,&kstat);
  if (kstat < 0) goto error;

  /* Get Euclidian norm of derivative */
  tprev = fabs(sval[1]);
  
  /* Compute the Newton stepdistanse vector. */
  s1172_s9dir(&td,sval);
  
  /* Adjust if we are not inside the parameter intervall. */
  t1 = td;
  s1172_s9corr(&t1,anext,astart,aend);
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 50; knbit++)
    {
      /* Evaluate 0-3.st derivatives of curve */
      
      tnext = anext + t1;
      
      s1221(pcurve,kder,tnext,&kleft,sval,&kstat);
      if (kstat < 0) goto error;

      /* Get Euclidian norm of derivative */
      tdist = fabs(sval[1]);
  
      /* Compute the Newton stepdistanse vector. */
      s1172_s9dir(&tdn,sval);
      
      /* Check if the direction of the step have change. */
      
      kdir = (td*tdn >= DZERO);     /* 0 if changed. */
      
      if (tdist <= tprev || kdir)
	{
	  /* Ordinary converging. */
      
          anext += t1;

          td = t1 = tdn;
	  
	  /* Adjust if we are not inside the parameter intervall. */
	  s1172_s9corr(&t1,anext,astart,aend);
	  
	  
          if (fabs(t1/tdelta) <= REL_COMP_RES)
	    {
	      anext += t1;
	      break;
	    }
	  
          tprev = tdist;
	}
      
      else
	{
	  /* Not converging, half step length try again. */
      
          t1 /= (double)2;
	  /*         knbit--;  */
	}
    }
  
  /* Iteration stopped, test if point is extremum */
  
  if (tdist <= tol)
    *jstat = 1;
  else
    *jstat = 0;

 
  /* Test if the iteration is close to a knot */
  if (fabs(anext - pcurve->et[kleft])/tdelta < tol)
    anext = pcurve->et[kleft];
  else if (fabs(anext - pcurve->et[kleft+1])/tdelta < tol)
    anext = pcurve->et[kleft+1];

  /* Uppdate output.  */
  *cpos = anext;
  
  /* Iteration completed.  */
  goto out;
  
 /* --------------------------------------------------------------------- */ 
  /* Error in input. Dimension not equal to 1 */
 err106: *jstat = -106;
  s6err("s1172",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("s1172",*jstat,kpos);
  goto out;                  
  
 out:;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1172_s9corr(double *cd, double acoef,double astart,double aend)
#else
static void s1172_s9corr(cd,acoef,astart,aend)
     double *cd;
     double acoef;
     double astart;
     double aend;
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
* INPUT      : acoef   - Parameter value to clipp
*              astart  - The lower boorder
*              aend    - The higher boorder
*
*
*
* INPUT/OUTPUT : cd    - Old and new step value.
*
*
* METHOD     : We are cutting a line.
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
  if (acoef + *cd < astart)  *cd = astart - acoef;
  else if (acoef + *cd > aend) *cd = aend - acoef;  
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1172_s9dir(double *cdiff,double evals[])
#else
static void s1172_s9dir(cdiff,evals)
     double *cdiff;
     double evals[];
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the step according to a Newton scheme
*              for finding an extremal to a one dimensional curve.
*
*
* INPUT      : evals - Value and derivatives in point on curve.
*              
*
* OUTPUT     : cdiff  - Parameter increment in first direction.
*            
*
*
* METHOD     : This is a one dimensional case. We want to find x such that
*
*                    x : S'(x) = 0      
*
*              Using Taylor we get:
*                  
*                    x : S'+dxS''+0.5 dxdxS''' = 0
*
*	       The solution of this equation is the
*	       following function.
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, OCTOBER 1990
*
*********************************************************************
*/                       
{                        
   double a,b,c,d,d1,d2;

   a = evals[3];
   b = evals[2];
   c = b*b - 2.0*a*evals[1];
   
   if (fabs(b) > DZERO)  d = -evals[1]/b;
   else                  d = 0.0;
   
   
   if (c < DZERO)                    *cdiff = d;
   else if (fabs(a) > DZERO)
   {
      c = sqrt(c);
      d1 = (-b + c)/a;
      d2 = (-b - c)/a;
      if (DEQUAL(b,c))               *cdiff = d;
      else
	if (fabs(d1-d) < fabs(d2-d)) *cdiff = d1;
      else                           *cdiff = d2;
   }
   else                              *cdiff = d;
}

