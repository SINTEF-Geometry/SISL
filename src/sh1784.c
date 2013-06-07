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
 * $Id: sh1784.c,v 1.4 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1784

#include "sislP.h"


typedef  void (*fevalcProc)(
#if defined(SISLNEEDPROTOTYPES)
                        SISLCurve *,
                        int,
                        double ,
                        int *,
                        double [],
                        int *
#endif
);


#if defined(SISLNEEDPROTOTYPES)
void
sh1784 (SISLCurve * pcurve, SISLSurf * psurf, double aepsge,
	double epar[], int icur, int idirc, double elast[],
	double enext[], int *jstat)
#else
void
sh1784 (pcurve, psurf, aepsge, epar, icur, idirc, elast, enext, jstat)
     SISLCurve *pcurve;
     SISLSurf *psurf;
     double aepsge;
     double epar[];
     int icur;
     int idirc;
     double elast[];
     double enext[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : March along a curve as long as it coincides with a given
*              surface. Start marching in an intersection point. Return
*              the last found point of coincidence, and the first found
*              point outside the coincidence interval.
*
*
* INPUT      : pcurve   - Pointer to the curve.
*              psurf    - Pointer to the surface.
*              aepsge   - Geometry resolution.
*              epar[3]  - Parameter values for the given intersection point.
*              icur     -         = 1 , epar[0] is the parameter value
*                                       for the curve and epar[1],epar[2] is
*                                       the parameter values for the surface.
*                                 = 0 , epar[2] is the parameter value
*                                       for the curve and epar[0],epar[1] is
*                                       the parameter values for the surface.
*              idirc    - Direction to march the curve.
*                         =  1 : March in the parameter direction.
*                         = -1 : March in the opposite direction.
*
*
*
* OUTPUT     : elast[3] - Parameter values of the last point within the
*                         interval of coincidence.
*              enext[3] - Parameter values of the first point found outside
*                         this interval.
*              jstat   -  status messages
*                                = 2   : Coincidence until an edge of the surface.
*                                = 1   : Coincidence of complete curve.
*                                = 0   : Ok. End of interval found.
*                                < 0   : Error.
*
*
* METHOD     : March along the curve and iterate down to the surface for
*              each midpoint and endpoint of each step. The steplenght
*              is computed from the curvature and the knot vector
*              of the curve. If any knot lines of the surface are crossed,
*              the marching are drawn back to the knot line, and a new
*              step is initiated.
*              If the closest point on the surface is at an edge of
*              the surface then we stop with jstat=2.
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, Oslo, Norway. 04.91.
* REVISED BY : Michael Floater,  SI, Oslo, Norway.  08.91.
*                 Stopping condition at edge of surface.
* CORRECTED BY: UJK,  SI, Oslo, Norway. Oct. 91
* REVISED BY : Vibeke Skytt, SI, 02.93. Reduce step length if curve
*                                       tangent changes direction.
*********************************************************************
*/
{
  int kstat;			/* Status variable                                 */
  int ki;			/* Counter.                                        */
  int kleftc = 0;		/* Left indicator for point calculation            */
  int kleft1 = 0;		/* Left indicator for point calculation in 1. par.
			           direction of surface.                           */
  int kleft2 = 0;		/* Left indicator for point calculation in 2. par dir.*/
  int kleft1prev, kleft2prev;	/* Previous left indicators of surface.    */
  int kn;			/* The number of B-splines, i.e., the dimension of
			           the spline space associated with the knot
			           vector.                                         */
  int kk;			/* The polynomial order of the curve.              */
  int kk1, kk2, kn1, kn2;	/* Orders and nu,ber of vertices of surface        */
  int kdimc;			/* The dimension of the space in which the curve
			           lies. Equivalently, the number of components
			           of each B-spline coefficient.                   */
  int kdims;			/* Dimension of space where the surface lies       */
  int kpos = 0;			/* Position of error                               */
  int kderc = 2;		/* Number of derivatives to be claculated on curve */
  int kders = 1;		/* Number of derivatives to be calculated on surface
			           If step lenght is to be generated from surface,
			           kders must be equal to 2.                       */
  int kdum;			/* Temporary variable                              */
  int kpar;			/* Parameter value of constant parameter curve.    */
  int kiterate;                 /* Indicates if further iteration is necessary
				   after curve-curve iteration.                    */
  double tref;                  /* Referance value in equality test.               */
  double tclose1, tclose2;	/* Parameter values of closest point between curves. */
  double tangdot;               /* Scalar product between curve tangents.          */
  double snorm[3];		/* Normal vector of surface                        */
  double s3dinf1[10];		/* Pointer to storage for point info of curve
			           (10 dobules prpoint when idim=3, 7 when idim=3) */
  double *st;			/* Pointer to the first element of the knot vector
			           of the curve. The knot vector has [kn+kk]
			           elements.                                       */
  double *st1;			/* First knot direction of surface                 */
  double *st2;			/* Second knot direction of surface                */
  double sfirst[2];		/* Start parameter par in surface                  */
  double tfirst;		/* Fist parameter on curve                         */
  double tend;			/* Last parameter on curve                         */
  double sderc[9];		/* Position, first and second derivative of curve  */
  double stangprev[3];          /* Previous tangent of curve.                      */
  double sders[18];		/* Position, first and second derivatives of surface */
  double tx, tx1, tx2;		/* Parameter value */
  double tstep;			/* Final step length     */
  double tmaxinc;		/* Maximal increment in parameter value along curve*/
  double tlengthend;		/* Length of 1st derivative at end of segment */
  double tincre;		/* Parameter value increment */
  double tsmax, tcmax;		/* Local maximal step length based of boxsizes of objects */
  double tdist = DZERO;		/* Distance */
  double sstart[2];		/* Lower boundary of parameter intervals */
  double send[2];		/* Upper bounadry of parameter intervals */
  double snext[3];		/* Existing iteration point on  surface            */
  double spos[3];		/* New iteration  point on surface                 */
  double snext2[2];		/* Help parameter values.                          */
  SISLPoint *qpoint = SISL_NULL;
  SISLCurve *qc = SISL_NULL;		/* Constant parameter curve.                       */

  /* Pointer to curve evaluator routine of the curve.  */

  fevalcProc fevalc;
/*
 #if defined(SISLNEEDPROTOTYPES)
   void (*fevalc) (SISLCurve *, int, double, int *, double[], int *);
 #else
   void (*fevalc) ();
 #endif
 */
  /* Make maximal step length based on box-size of surface */

  sh1992su (psurf, 0, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  tsmax = MAX (psurf->pbox->e2max[0][0] - psurf->pbox->e2min[0][0],
	       psurf->pbox->e2max[0][1] - psurf->pbox->e2min[0][1]);
  tsmax = MAX (tsmax, psurf->pbox->e2max[0][2] - psurf->pbox->e2min[0][2]);

  /* Make maximal step length based on box-size of curve */

  sh1992cu (pcurve, 0, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  tcmax = MAX (pcurve->pbox->e2max[0][0] - pcurve->pbox->e2min[0][0],
	       pcurve->pbox->e2max[0][1] - pcurve->pbox->e2min[0][1]);
  tcmax = MAX (tcmax, pcurve->pbox->e2max[0][2] - pcurve->pbox->e2min[0][2]);

  /* Copy curve attributes to local parameters.  */

  kdimc = pcurve->idim;
  kk = pcurve->ik;
  kn = pcurve->in;
  st = pcurve->et;

  /* Copy surface attributes to local parameters.  */

  kdims = psurf->idim;
  kk1 = psurf->ik1;
  kk2 = psurf->ik2;
  kn1 = psurf->in1;
  kn2 = psurf->in2;
  st1 = psurf->et1;
  st2 = psurf->et2;

  /* Set reference value.  */

  tref = MAX(st[kn]-st[kk-1],MAX(st1[kn1]-st1[kk1-1],st2[kn2]-st2[kk2-1]));

  /* Check that dimensions are 3 */

  if (kdimc != 3 || kdims != 3)
    goto err105;

  sstart[0] = st1[kk1 - 1];
  sstart[1] = st2[kk2 - 1];
  send[0] = st1[kn1];
  send[1] = st2[kn2];

  /* Copy interval description into local variables */

  if (icur == 1)
    {
      sfirst[0] = epar[1];
      sfirst[1] = epar[2];
      tfirst = epar[0];
      tend = (idirc == 1) ? st[kn] : st[kk - 1];
    }
  else
    {
      sfirst[0] = epar[0];
      sfirst[1] = epar[1];
      tfirst = epar[2];
      tend = (idirc == 1) ? st[kn] : st[kk - 1];
    }

  /* To make sure we do not start outside or end outside the curve we
     truncate tfirst to the knot interval of the curve */

  tfirst = (idirc == 1) ? MAX (tfirst, st[kk - 1]) : MIN (tfirst, st[kn]);

  /* Set start point of iteration on surface */

  spos[0] = sfirst[0];
  spos[1] = sfirst[1];

  /* Set curve evaluator of the curve.  */

  fevalc = (idirc == 1) ? s1221 : s1227;

  /* Store knot values at start of curve */

  tx2 = tfirst;
  kdum = MAX (kk1, kk2);
  kdum = MAX (kdum, kk);
  tmaxinc = fabs (tend - tfirst) / (kdum * kdum);

  /* Make start point of curve  */

  fevalc (pcurve, kderc, tx2, &kleftc, sderc, &kstat);
  if (kstat < 0) goto error;

  /* Make start point of surface.  */

  s1421 (psurf, kders, spos, &kleft1, &kleft2, sders, snorm, &kstat);
  if (kstat < 0) goto error;

  /* While end not reached */

  while (idirc * tx2 < idirc * tend)
    {
      /* Save parameters of previous step.   */

      tx1 = tx2;
      snext[0] = spos[0];
      snext[1] = spos[1];
      kleft1prev = kleft1;
      kleft2prev = kleft2;

      /* Calculate unit tangent and radius of curvature of curve. */

      s1307 (sderc, kdimc, s3dinf1, &kstat);
      if (kstat < 0)
	goto error;

      /* Calculate step length based on curvature */

      tstep = s1311 (s3dinf1[3 * kdimc], aepsge, tsmax, &kstat);
      if (kstat < 0)
	goto error;

      /* Remember length of start tangent, end of zero segment */

      tlengthend = s6length (sderc + kdimc, kdimc, &kstat);
      if (kstat < 0)
	goto error;


      /* Find candidate end point, make sure that no breaks in tangent or
         curvature exists between start and endpoints of the segment     */
      /* Make step length equal to resolution if the length is zero */

      if (DEQUAL (tlengthend, DZERO))
	tincre = REL_PAR_RES;
      else
	tincre = tstep / tlengthend;

      tincre = MIN (tincre, tmaxinc);

      /*  Make sure that we don't pass any knots of the curve. */

      if (idirc * (tx1 + tincre) > idirc * (st[kleftc + idirc] + REL_PAR_RES))
	tincre = idirc * (st[kleftc + idirc] - tx1);

      if (idirc < 0 && (tx1 - tincre < st[kleftc] - REL_PAR_RES))
	tincre = idirc * (st[kleftc] - tx1);

      /* Find parameter value of candidate end point of segment */

      tx2 = tx1 + idirc * tincre;

      for (ki = 0, tx = (tx1 + tx2) / (double) 2.0; ki < 2; ki++, tx = tx2)
	{
	  if (idirc * tx >= idirc * tend)
	    break;

	  /* Make point sderc at curve at tx */

	  fevalc (pcurve, kderc, tx, &kleftc, sderc, &kstat);
	  if (kstat < 0) goto error;

	  /* Test if the step is legal.  */

	  if (DNEQUAL(tx1,tfirst) || ki>0)
	  {
	     tangdot = s6scpr(stangprev, sderc+kdimc, kdimc);
	     while (tangdot < DZERO)
	     {
		/* The step is not legal. Reduce step length. */

		if (ki == 0)
		{
		   tx2 = tx;
		   tx = (tx1 + tx2)/(double)2.0;
		}
		else
		{
		   tx2 = tx1 + (double)0.75*(tx2-tx1);
		   tx = tx2;
		}

		/* Make point sderc at curve at tx */

		fevalc (pcurve, kderc, tx, &kleftc, sderc, &kstat);
		if (kstat < 0) goto error;

		tangdot = s6scpr(stangprev, sderc+kdimc, kdimc);
	     }
	  }
	  /* Find closest point on surface to sderc */

	  qpoint = newPoint (sderc, kdimc, 0);
	  if (qpoint == SISL_NULL)
	    goto err101;

	  snext2[0] = snext[0];
	  snext2[1] = snext[1];
	  s1773 (qpoint, psurf, aepsge, sstart, send, snext2, spos, &kstat);
	  if (kstat < 0)
	    goto error;

	  freePoint (qpoint);
	  qpoint = SISL_NULL;

	  /* Check to see if we have crossed an edge of the
	     surface, i.e. we have gone outside the parameter
	     area for psurf. */

          if(spos[0] <= st1[kk1-1] || spos[0] >= st1[kn1] ||
             spos[1] <= st2[kk2-1] || spos[1] >= st2[kn2])
          {
	      /* Coincidence! Finish with a message. */
	      goto edge_of_surf;
	  }

	  /* Calculate point and derivatives in surface */

	  s1421 (psurf, kders, spos, &kleft1, &kleft2, sders, snorm, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* Check if point on curve and surface are within positional and
             angular tolerances */

	  tdist = s6dist (sderc, sders, kdimc);

	  if (tdist > aepsge)
	    {
	      /* Points not within tolerances, curve and surface do not
	         coincide */
	      goto no_coin;
	    }

	  /* Check if any parameter lines of the surface is crossed in the 1.
             parameter direction.  */

	  if (kleft1 != kleft1prev &&
	      ((DNEQUAL(spos[0]+tref,st1[kleft1]+tref) &&
		DNEQUAL(snext[0]+tref,st1[kleft1]+tref)) ||
	       kleft1 != kleft1prev+1) &&
	      ((DNEQUAL(snext[0]+tref,st1[kleft1prev]+tref) &&
		DNEQUAL(spos[0]+tref,st1[kleft1prev]+tref)) ||
	       kleft1 != kleft1prev - 1))
	    {
	      /* At least one parameter line is crossed. Fetch the constant parameter
	         curve at the closest parameter line in the direction of the marching. */

	      if (kleft1 > kleft1prev)
		kpar = kleft1prev + 1;
	      else if (snext[0] != st1[kleft1prev])
		kpar = kleft1prev;
	      else
		kpar = kleft1prev - 1;

	      /* Pick constant parameter curve.   */

	      s1437 (psurf, st1[kpar], &qc, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Find the closest point between the input curve and the constant
	         parameter curve.    */

		/* UJK Oct 91, Nice trap ! tx1 > tx */
		/*  s1770 (pcurve, qc, aepsge, tx1, st2[kk2 - 1], tx, st2[kn2], (tx1 + tx) / (double) 2.0,
		     st2[kleft2], &tclose1, &tclose2, &kstat); */
	       s1770 (pcurve, qc, aepsge, min(tx1,tx), st2[kk2 - 1], max(tx1,tx),
		      st2[kn2], (tx1 + tx) / (double) 2.0,
		     (double)0.5*(st2[kleft2]+st2[kleft2+1]),
		     &tclose1, &tclose2, &kstat);
	      if (kstat < 0)
		goto error;

	      if (kstat == 2 || fabs(tclose1-tx1) < REL_PAR_RES)
		 /* No intersection point is found. Mark that surface-point
		    iteration is necessary.  */

		 kiterate = 1;
	      else kiterate = 0;

	      /* Set new parameter values to the iteration.  */

	      spos[0] = st1[kpar];
	      spos[1] = tclose2;
	      if (fabs(tclose1-tx1) > REL_PAR_RES) tx2 = tclose1;

	      /* Test midpoint of reduced step. First evaluate curve in midpoint. */

	      tx = (tx1 + tx2) / (double) 2.0;

	      fevalc (pcurve, kderc, tx, &kleftc, sderc, &kstat);
	      if (kstat < 0) goto error;

	      /* Find closest point on surface to sderc */

	      qpoint = newPoint (sderc, kdimc, 0);
	      if (qpoint == SISL_NULL)
		goto err101;

	      snext2[0] = snext[0];
	      snext2[1] = snext[1];
	      s1773 (qpoint, psurf, aepsge, sstart, send, snext2, snext2, &kstat);
	      if (kstat < 0)
		goto error;

	      freePoint (qpoint);
	      qpoint = SISL_NULL;

	      /* Calculate point and derivatives in surface */

	      s1421 (psurf, kders, snext2, &kleft1, &kleft2, sders, snorm, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Check if point on curve and surface are within positional and
	         angular tolerances */

	      tdist = s6dist (sderc, sders, kdimc);

	      if (tdist > aepsge)
		{
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto no_coin;
		}

	      /* Calculate point and derivatives in the curve in the endpoint of the step. */

	      fevalc (pcurve, kderc, tx2, &kleftc, sderc, &kstat);
	      if (kstat < 0) goto error;

	      if (kiterate)
	      {
		 /* Relax the point on the curve down to the surface. */

		 qpoint = newPoint (sderc, kdimc, 0);
		 if (qpoint == SISL_NULL)
		    goto err101;

		 spos[0] = snext2[0];
		 spos[1] = snext2[1];
		 s1773 (qpoint, psurf, aepsge, sstart, send, spos, spos, &kstat);
		 if (kstat < 0)
		    goto error;

		 freePoint (qpoint);
		 qpoint = SISL_NULL;
	      }

	      /* Calculate point and derivatives in the surface.  */

	      s1421 (psurf, kders, spos, &kleft1, &kleft2, sders, snorm, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Check if point on curve and surface are within positional and
	         angular tolerances */

	      tdist = s6dist (sderc, sders, kdimc);

	      if (tdist > aepsge)
		{
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto no_coin;
		}

	      /* Mark that a new step is to be initiated.  */

	      ki = 2;

	      /* Free constant parameter curve.  */

	      if (qc != SISL_NULL)
		freeCurve (qc);
	      qc = SISL_NULL;
	    }

	  /* Check if any parameter lines of the surface is crossed in the 2.
             parameter direction.  */

	  if (kleft2 != kleft2prev &&
	      ((DNEQUAL(spos[1]+tref,st2[kleft2]+tref) &&
		DNEQUAL(snext[1]+tref,st2[kleft2]+tref)) ||
	       kleft2 != kleft2prev+1) &&
	      ((DNEQUAL(snext[1]+tref,st2[kleft2prev]+tref) &&
		DNEQUAL(spos[1]+tref,st2[kleft2prev]+tref)) ||
	       kleft2 != kleft2prev - 1))
	    {
	      /* At least one parameter line is crossed. Fetch the constant parameter
	         curve at the closest parameter line in the direction of the marching. */

	      if (kleft2 > kleft2prev)
		kpar = kleft2prev + 1;
	      else if (snext[1] != st2[kleft2prev])
		kpar = kleft2prev;
	      else
		kpar = kleft2prev - 1;

	      /* Pick constant parameter curve.   */

	      s1436 (psurf, st2[kpar], &qc, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Find the closest point between the input curve and the constant
	         parameter curve.    */

		/* UJK Oct 91, Nice trap ! tx1 > tx */
		s1770 (pcurve, qc, aepsge, min(tx1,tx), st1[kk1 - 1], max(tx,tx1),
		       st1[kn1], (tx1 + tx) / (double) 2.0,
		       (double)0.5*(st1[kleft1]+st1[kleft1+1]),
		       &tclose1, &tclose2, &kstat);
	      if (kstat < 0)
		goto error;

	      if (kstat == 2 || fabs(tclose1-tx1) < REL_PAR_RES)
		 /* No intersection point is found. Mark that surface-point
		    iteration is necessary.  */

		 kiterate = 1;
	      else kiterate = 0;

	      /* Set new parameter values to the iteration.  */

	      spos[0] = tclose2;
	      spos[1] = st2[kpar];
	      if (fabs(tclose1-tx1) > REL_PAR_RES) tx2 = tclose1;

	      /* Test midpoint of reduced step. First evaluate curve in midpoint. */

	      tx = (tx1 + tx2) / (double) 2.0;

	      fevalc (pcurve, kderc, tx, &kleftc, sderc, &kstat);
	      if (kstat < 0) goto error;

	      /* Find closest point on surface to sderc */

	      qpoint = newPoint (sderc, kdimc, 0);
	      if (qpoint == SISL_NULL)
		goto err101;

	      snext2[0] = snext[0];
	      snext2[1] = snext[1];
	      s1773 (qpoint, psurf, aepsge, sstart, send, snext2, snext2, &kstat);
	      if (kstat < 0)
		goto error;

	      freePoint (qpoint);
	      qpoint = SISL_NULL;

	      /* Calculate point and derivatives in surface */

	      s1421 (psurf, kders, snext2, &kleft1, &kleft2, sders, snorm, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Check if point on curve and surface are within positional and
	         angular tolerances */

	      tdist = s6dist (sderc, sders, kdimc);

	      if (tdist > aepsge)
		{
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto no_coin;
		}

	      /* Calculate point and derivatives in the curve.    */

	      fevalc (pcurve, kderc, tx2, &kleftc, sderc, &kstat);
	      if (kstat < 0) goto error;

	      if (kiterate)
	      {
		 /* Relax the point on the curve down to the surface. */

		 qpoint = newPoint (sderc, kdimc, 0);
		 if (qpoint == SISL_NULL)
		    goto err101;

		 spos[0] = snext2[0];
		 spos[1] = snext2[1];
		 s1773 (qpoint, psurf, aepsge, sstart, send, spos, spos, &kstat);
		 if (kstat < 0)
		    goto error;

		 freePoint (qpoint);
		 qpoint = SISL_NULL;
	      }


	      /* Calculate point and derivatives in the surface.  */

	      s1421 (psurf, kders, spos, &kleft1, &kleft2, sders, snorm, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Check if point on curve and surface are within positional and
	         angular tolerances */

	      tdist = s6dist (sderc, sders, kdimc);

	      if (tdist > aepsge)
		{
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto no_coin;
		}

	      /* Mark that a new step is to be initiated.  */

	      ki = 2;

	      /* Free constant parameter curve.  */

	      if (qc != SISL_NULL)
		freeCurve (qc);
	      qc = SISL_NULL;
	    }

	  /* Save tangent of curve.  */

	  memcopy(stangprev, sderc+kdimc, kdimc, DOUBLE);
	}
    }

  /* Coincidence interval along complete curve. */

  *jstat = 1;
  if (icur == 1)
    {
      elast[0] = tx1;
      elast[1] = snext[0];
      elast[2] = snext[1];
    }
  else
    {
      elast[0] = snext[0];
      elast[1] = snext[1];
      elast[2] = tx1;
    }
  goto out;

  /* Curve and surface not within tolerance */
no_coin:*jstat = 0;
  if (icur == 1)
    {
      elast[0] = tx1;
      elast[1] = snext[0];
      elast[2] = snext[1];
      enext[0] = tx2;
      enext[1] = spos[0];
      enext[2] = spos[1];
    }
  else
    {
      elast[0] = snext[0];
      elast[1] = snext[1];
      elast[2] = tx1;
      enext[0] = spos[0];
      enext[1] = spos[1];
      enext[2] = tx2;
    }
  goto out;

  /* Curve and surface are within tolerance up to an edge
     of the surface. */
edge_of_surf:
  *jstat = 2;
  if (icur == 1)
    {
      elast[0] = tx1;
      elast[1] = snext[0];
      elast[2] = snext[1];
      enext[0] = tx2;
      enext[1] = spos[0];
      enext[2] = spos[1];
    }
  else
    {
      elast[0] = snext[0];
      elast[1] = snext[1];
      elast[2] = tx1;
      enext[0] = spos[0];
      enext[1] = spos[1];
      enext[2] = tx2;
    }
  goto out;

  /* Error in memory allocation */

err101:*jstat = -101;
  s6err ("sh1784", *jstat, kpos);
  goto out;

  /* Error in input, dimension not equal to 2 or 3 */

err105:*jstat = -105;
  s6err ("sh1784", *jstat, kpos);
  goto out;

  /* Error in lower level function */

error:*jstat = kstat;
  s6err ("sh1784", *jstat, kpos);
  goto out;


out:

  return;
}
