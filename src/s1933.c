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
 * $Id: s1933.c,v 1.4 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1933

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1933 (int inbcrv, SISLCurve * crvarr[], double start, double stop,
       double **it, int *in, int *iordr, int *jstat)
#else
void
s1933 (inbcrv, crvarr, start, stop, it, in, iordr, jstat)
     int inbcrv;
     SISLCurve *crvarr[];
     double start;
     double stop;
     double **it;
     int *in;
     int *iordr;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To map a number of knot vectors of different lengths
*		into the same interval, and to produce a union of the
*		knot vectors taking care to use the greatest of the
*		orders of the input knot vectors and to reflect the
*		continuity requirements reflected in all input knot
*		vectors. To avoid too many new knots, knots lying
*		closer than 0.0001 times the parameter interval
*               (relative distance), might be moved.
*
*
* INPUT      :	inbcrv	- Number of curves to be interpolated by a
*			  spline lofted surface.
*		crvarr	- Array (length inbcrv) of pointers to the curves
*			  in the curve-set.
*		start	- Start parameter value of B-spline basis to
*			  be made.
*		stop	- End parameter value of B-spline basis to
*			  be made.
*
*
* OUTPUT     :  it	- The new knot vector (union of the old ones).
*		in	- The number of B-spline basis functions in basis
*			  produced.
*		iordr	- The order of the B-spline basis produced.
*               jstat   - Output status:
*                          < 0: Error.
*                          = 0: Ok.
*                          > 0: Warning.
*
* METHOD     : 	Phase 1:
*                   Find highest order of all curves.
*                   Map all knot-vectors into [astart,astop].
*                   Lift all knot-vectors to the highest order and make union.
*           	Phase 2:
*                   Run through the union vector and detect knots lying
*                   closer to each other than 10E-4 times the parameter
*                   interval.
*                   Find values to be used for the close knot-values.
*                   Move knots in original knot-vectors to match these values.
*           	Phase 3:
*                   Lift all knot-vectors (now adjusted) to the highest
*                   order and make union.
*
* REFERENCES :  Fortran version:
*               T.Dokken, SI, 1981-12
*
* CALLS      :  s1934,s1754,s1935,s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07
* REVISED BY :  Paal Fugelli, SINTEF, Oslo 18/07-1994.  Removed memory
*               leaks.
*
*********************************************************************
*/
{
  int ki, kj, kl, kr, kp;	/* Loop control variables		*/
  int ktell;
  int kfirst;			/* Pointer to first knot value in
				 * interval [tinmin,tinmax]		*/
  int klast;			/* Pointer to last knot value in
				 * interval [tinmin,tinmax]		*/
  int kn;			/* Number of vertices in curve		*/
  int kbgn;			/* Pointer to first knot value
				 * satisfying kt[kbgn] > kt[*iordr-1]	*/
  int kend;			/* Pointer to last knot value
				 * satisfying kt[kend] < kt[kn]		*/
  int kmax;			/* Number of distinct knot values to be
				 * kept in the interval [tinmin,tinmax]	*/
  int kdum1;			/* Dummy variables for use in algorithm */
  int kdum2;
  int kdiff;
  int kdf;			/* Equals (kdiff div 2)			*/
  int knumb;			/* Number of vertices in a curve	*/
  int kpos = 0;			/* Error position indicator 		*/
  int kstat;
  double tepsco;		/* Greatest distance for two knots to be
				 * regarded as the same knot		*/
  double dum;
  double tlast;			/* Value of knot vector			*/
  double tval;
  double tincre;
  double tincr2;		/* Equals   tincre / 2.0		*/
  double tmin, tmax;
  double tinmin, tinmax;
  double tbgn, tend;		/* tbgn=kt[*iordr-1] & tend = kt[kn]
				 * Used to find kbgn and kend		*/
  double *kt = SISL_NULL;		/* New knot vector			*/
  double *knot = SISL_NULL;		/* Knot vector				*/
  double *incknt = SISL_NULL;	/* Used to store the union of two
				 * knot vectors				*/
  SISLCurve *curve = SISL_NULL;

  *jstat = 0;


  /* Initailzation of variables */

  *in = 0;
  *iordr = 0;


  /* Test if legal input */

  if (inbcrv < 2)
    goto err179;
  for (ki = 0; ki < inbcrv; ki++)
    if (crvarr[ki]->in <crvarr[ki]->ik || crvarr[ki]->ik < 1)
      goto err112;


  /* Find highest order used in description */

  for (ki = 0; ki < inbcrv; ki++)
    *iordr = MAX (*iordr, crvarr[ki]->ik);

  kn = *iordr;


  /* Normalize the knot vectors */

  for (ki = 0; ki < inbcrv; ki++)
    {
      curve = crvarr[ki];
      s1934 (curve->et, curve->in, curve->ik, start, stop, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* Make start knot vector */

  kt = newarray (*iordr * 2, double);
  if (kt == SISL_NULL)
    goto err101;
  for (ki = 0; ki < *iordr; ki++)
    {
      kt[ki] = start;
      kt[*iordr + ki] = stop;
    }

  /* Run through all knot vectors, lift the order, map into right
   * interval and make the union with the candidates already existing */

  for (ki = 0; ki < inbcrv; ki++)
    {
      /* Increase order of basis */

      curve = crvarr[ki];
      s1754 (curve->et, curve->in, curve->ik, *iordr,
	     &incknt, &knumb, &kstat);
      if (kstat < 0)
	goto error;


      /* Make union of old union and new knot vector */

      s1935 (kt, kn, incknt, knumb, &knot, &kn, *iordr, &kstat);
      if (kstat < 0)
	goto error;

      if (incknt != SISL_NULL)  freearray(incknt);  /* PFU 18/07-94. */

      if (kt != SISL_NULL)
	freearray (kt);
      kt = knot;
      knot = SISL_NULL;  /* PFU 18/07-94 */
    }

  /* The knot vector produced might contain knots originating
   * from different curves that are located very close to each
   * other. We will move internal knots when they are closer
   * to each other than tepsco.
   *
   * Find first knot bigger than kt[*iordr-1] and last
   * knot less than kt[kn].					*/

  kend = kn - 1;
  tbgn = kt[*iordr - 1];
  tend = kt[kn];

  /* Find first knot kt[kbgn] > tbgn	*/

  for (kbgn = *iordr; (tbgn >= kt[kbgn]) && (kbgn < kend); kbgn++) ;

  /* Find last knot kt[kend] < tend	*/

  for (; (tend <= kt[kend]) && (kbgn < kend); kend--) ;

  /* Set greatest distance for two knots to be regarded
   * as the same knot					*/

  tepsco = (double) 0.00000000000001 * (stop - start);

  /* Loop runing through all the knots of the union knot
   * vector and moving knots when possible		*/

  for (ki = kbgn; ki < kend; ki = kl + 1)
    {
      /* Find knots closer to knot ki than tepsco in
       * positive direction				*/

      kl = ki - 1;
      do
	{
	  kl++;
	  dum = MAX (fabs(kt[kl + 1]), fabs(kt[ki]))/stop;
	  if (dum == (double) 0.0)
	    dum = (double) 1.0;
	}
      while ((kl < kend) && ((fabs (kt[kl + 1] - kt[ki]) / dum) <= tepsco));


      /* If only one value found, then shifting is not necessary */

      if (ki == kl)
	continue;


      /* The knots close to kt[ki-1] are kt[ki]...kt[kl-1]	*/

      tmin = kt[ki];
      tmax = MIN (kt[kl] + tepsco*dum, kt[kend]);


      /* If interval is degenerate we have finished the moving	*/

      if (tmin >= tmax)
	break;


      /* For each curve, count the number of distinct knot
         values in the interval [tmin,tmax].
         Accumulate max and min knot values in the interval
         where we know that tmin is a knot.			*/

      tinmin = tmin;
      tinmax = tmax;
      kmax = 0;
      for (kj = 0; kj < inbcrv; kj++)
	{
	  curve = crvarr[kj];
	  ktell = 0;
	  tlast = (double) 0.0;

	  for (kr = *iordr - 1; kr <= curve->in; kr++)
	    if ((tmin <= curve->et[kr]) && (curve->et[kr] <= tmax))
	      {
		if (ktell == 0)
		  ktell = 1;
		else if (curve->et[kr] > tlast)
		  ktell++;
		tlast = curve->et[kr];
		tinmin = MIN (tinmin, tlast);
		tinmax = MAX (tinmax, tlast);
	      }
	  kmax = MAX (kmax, ktell);
	}

      /* kmax now contains the number of distinct knot values
       * to be kept in the interval [tinmin,tinmax]. */

      if (kmax > 1)
	tincre = (tinmax - tinmin) / (double) (kmax - 1);
      else
	tincre = (double) 0.0;
      tincr2 = tincre / (double) 2.0;

      /* The values used will be
         tinmin,tinmin+tincre,...,tinmin+(kmax-1)*tincre.
         We want to use the values closest to the actual
         knot values when possible.
         Run through each curve and move knots if tinmin<tinmax.
         If they are equal, knots should not be moved. */

      if (tinmin >= tinmax)
	continue;

      for (kj = 0; kj < inbcrv; kj++)
	{
	  curve = crvarr[kj];
	  ktell = 0;
	  for (kr = *iordr - 1; kr <= curve->in; kr++)
	    if ((tinmin <= curve->et[kr]) && (curve->et[kr] <= tinmax))
	      {
		if (ktell == 0)
		  {
		    ktell = 1;
		    kfirst = kr;
		  }
		else if (curve->et[kr] >= tlast)
		  {
		    ktell++;
		    klast = kr;
		  }
		tlast = curve->et[kr];
	      }

	  /* ktell contains the number of knots in [tinmin,tinmax]
	     on curve kj. kfirst is a pointer to first knot value,
	     klast is the pointer to last knot value		*/

	  if (ktell == 1)
	    {

	      /* Only one knot value in interval, move to closest
	       * legal value					*/

	      if ((kmax == 1) || (tincre == (double) 0.0))
		tval = tinmin;
	      else
		{
		  kdum1 = (int)((curve->et[kfirst] - tinmin + tincr2) / tincre);
		  tval = tinmin + (double)kdum1 * tincre;
		}
	      tlast = curve->et[kfirst];
	      do
		{
		  curve->et[kfirst] = tval;
		  kfirst++;
	      } while (curve->et[kfirst] == tlast);
	      continue;
	    }

	  if ((ktell <= 1) || (tincre <= (double) 0.0))
	    continue;


	  /* More than one point found. */

	  kdum1 = (int)((curve->et[kfirst] - tinmin + tincr2) / tincre);
	  kdum2 = (int)((curve->et[klast] - tinmin + tincr2) / tincre);
	  kdiff = ktell - kdum2 + kdum1 - 1;
	  if (kdiff > 0)
	    {
	      /* Change kdum1 and kdum2 to allow for ktell
	         different values. */

	      kdf = kdiff / 2;
	      kdum1 = MAX (0, kdum1 - kdf);
	      kdum2 = MIN (kmax - 1, kdum2 + kdiff - kdf);
	      kdiff = ktell - kdum2 + kdum1 - 1;

	      if (kdiff > 0)
		{
		  if (kdum1 == 0)
		    kdum2 += kdiff;
		  else
		    {
		      if (kdum2 != ktell)
			goto err170;
		      kdum1 -= kdiff;
		    }
		}
	    }

	  if ((kdum1 < 0) || (kdum2 > ktell))
	    goto err170;


	  /* Use kdum1 as start value, move knots. */

	  tval = tinmin + kdum1 * tincre;
	  kr = kfirst;
	  tlast = curve->et[kr];
	  for (kp = 0; kp < ktell; kp++, kr++)
	    if (curve->et[kr] == tlast)
	      curve->et[kr] = tval;
	    else
	      {
		if (curve->et[kr] <= tlast)
		  goto err170;
		tval += tincre;
		tlast = curve->et[kr];
		curve->et[kr] = tval;
	      }
	}
    }

  /* Make start knot vector					*/

  if (kt != SISL_NULL)  freearray(kt);  /* PFU 18/07-94 */

  kt = newarray (*iordr * 2, double);
  if (kt == SISL_NULL)
    goto err101;

  for (ki = 0; ki < *iordr; ki++)
    {
      kt[ki] = start;
      kt[*iordr + ki] = stop;
    }

  kn = *iordr;

  /* Run through all knot vectors, lift the order, map into right
     interval and make the union with the candidates already existing. */

  for (ki = 0; ki < inbcrv; ki++)
    {
      /* Increase order of basis. */

      curve = crvarr[ki];
      s1754 (curve->et, curve->in, curve->ik, *iordr, &incknt, &knumb, &kstat);
      if (kstat < 0)
	goto error;


      /* Make union of old union and new knot vector */

      s1935 (kt, kn, incknt, knumb, &knot, &kn, *iordr, &kstat);
      if (kstat < 0)
	goto error;

      if (incknt != SISL_NULL)  freearray(incknt);  /* PFU 18/07-94 */

      if (kt != SISL_NULL)
	freearray (kt);
      kt = knot;
      knot = SISL_NULL;  /* PFU 18/07-94 */
    }

  /* No errors */

  *in = kn;
  *it = kt;
  kt = SISL_NULL;  /* PFU 18/07-94 */
  goto out;


  /* Memory error */

err101:
  *jstat = -101;
  s6err ("s1933", *jstat, kpos);
  goto out;

  /* Error in input */

err112:
  *jstat = -112;
  s6err ("s1933", *jstat, kpos);
  goto out;

  /* Too few curves for spline lofted surface */

err179:
  *jstat = -179;
  s6err ("s1933", *jstat, kpos);
  goto out;

  /* Error in lower level routines */

error:
  *jstat = kstat;
  s6err ("s1933", *jstat, kpos);
  goto out;

  /* Special error in moving of knot values */

err170:
  *jstat = -170;
  s6err ("s1933", *jstat, kpos);
  goto out;

out:
  if (kt != SISL_NULL)  freearray(kt);  /* PFU 18/07-94 */
  if (knot != SISL_NULL)  freearray(knot);  /* PFU 18/07-94 */
  if (incknt != SISL_NULL)  freearray(incknt);  /* PFU 18/07-94 */
  return;
}
