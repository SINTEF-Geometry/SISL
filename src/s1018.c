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
 * $Id: s1018.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1018

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1018 (SISLCurve *pc, double epar[], int inpar,
		  SISLCurve **rcnew, int *jstat)
#else
void
s1018 (pc, epar, inpar, rcnew, jstat)
     SISLCurve *pc;
     double epar[];
     int inpar;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Insert a given set of knots into the description of
*              a B-spline curve.
* NOTE       : When the curve is periodic, the input parameter values
*              must lie in the HALFOPEN [et[kk-1], et[kn), the function
*              will automatically update the extra knots and
*              coeffisients.
*              rcnew->in is still eq to pc->in + inpar !
*
* INPUT      : pc        - SISLCurve to be refined.
*              epar      - Parameter values of knots to be inserted.
*                          The values are stored in increasing order.
*              inpar     - Number of parameter values in epar.
*
*
*
* OUTPUT     : rcnew     - The new, refined curve .
*              jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              S1701.C   - Making the knot-inserten-transformation matrix.
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* CHANGED BY : Ulf J. Krystad, SI, 92-01
*              Treatment of periodic curves.
**********************************************************************/
{
  register int ki, ki2;  	/* Control variable in loop.                */
  register int kj, kj1, kj2;	/* Control variable in loop.                */
  register double *s1;		/* Pointers used in loop.                     */
  int kstat;			/* Local status variable.                     */
  int k1, k2, k3;		/* Index to array of parameter values.        */
  int kmy;			/* An index to the knot-vector.               */
  int kpl, kfi, kla;		/* To posisjon elements in trans.-matrix.     */
  int kk = pc->ik;		/* Order of the input curve.                  */
  int kn = pc->in;		/* Number of the vertices in input curves.    */
  int kdim = pc->idim;		/* Dimensjon of the space in whice curve lies.*/
  int kn1;			/* Number of vertices in the new curve .      */
  double tstart, tend;		/* Endparameters of curve.                    */
  double tpar;			/* Parameter value of knot to insert.         */
  double *st = SISL_NULL;		/* The new knot-vector.                       */
  double *sp = SISL_NULL;		/* To use in s1701.c                          */
  double *salfa = SISL_NULL;		/* A line of the trans.-matrix.               */
  double *scoef = SISL_NULL;		/* The new vertices.                          */
  double *coef = SISL_NULL;		/* The old vertices.                          */
  SISLCurve *q1 = SISL_NULL;		/* Pointer to new curve-object.               */
  /* Periodicity treatment --------------------------                         */
  double *st2 = SISL_NULL;		/* The new knot-vector.                    */
  double *scoef2 = SISL_NULL;	/* The new vertices.                       */
  double t1;			/* Start of knot vector                    */
  double t2;			/* Start of full basis part of knot vector */
  double t3;			/* End  of full basis part of knot vector  */
  double t4;			/* End of knot vector                      */
  double tmod;			/* t3 - t2, period size                    */
  double mod_neg;		/* parametervalue shifted tmod left        */
  double mod_pos;		/* parametervalue shifted tmod right       */
  int no_neg, no_pos;		/* Nmb of Xtra parval. to insert (lft,rght)*/
  double *negpar = SISL_NULL;	/* Array of parameter values shift lft  */
  double *pospar = SISL_NULL;	/* Array of parameter values shift rght */
  double *pararr = SISL_NULL;	/* Pointer to parameter array treated   */
  double *periodarr = SISL_NULL;	/* Total parameter array, periodic case*/
  /* --------------------------------------------------------- */

  /* Check that we have a curve to treat. */

  if (!pc)
    goto err150;
  if (pc->ikind == 2 || pc->ikind == 4)
    {
      kdim++;
      coef = pc->rcoef;
    }
  else
    {
      coef = pc->ecoef;
    }


  /* Periodicity treatment -------------------------- */
  if (pc->cuopen == SISL_CRV_PERIODIC)
    {
       /* Get values in knotvector for start, start of full basis,
	  end of full basis, end. */
      t1 = *(pc->et);
      t2 = *(pc->et + kk - 1);
      t3 = *(pc->et + kn);
      t4 = *(pc->et + kn + kk - 1);
      tmod = t3 - t2;

      /* Check that the input parameter values lie in the legal parameter
	 interval of the curve, ie in the HALFOPEN [et[kk-1], et[kn) */
      for (k1 = 0; k1 < inpar; k1++)
	if (epar[k1] < t2 ||epar[k1] >= t3)
	   goto err158;

      no_neg = 0;
      no_pos = 0;
      if ((negpar = newarray (inpar, double)) == SISL_NULL)
	goto err101;
      if ((pospar = newarray (inpar, double)) == SISL_NULL)
	goto err101;


      for (k1 = 0; k1 < inpar; k1++)
	{
	  mod_neg = epar[k1] - tmod;
	  mod_pos = epar[k1] + tmod;
	  if (mod_neg <= t2 &&
	      mod_neg > t1)
	    {
	      negpar[no_neg] = mod_neg;
	      no_neg++;
	    }
	    if (mod_pos >= t3 &&
		mod_pos <t4)
	    {
	       pospar[no_pos] = mod_pos;
	       no_pos++;
	    }
	}

      if ((periodarr = newarray (inpar + no_neg + no_pos, double)) == SISL_NULL)
	goto err101;
      memcopy (periodarr, negpar, no_neg, double);
      memcopy (periodarr + no_neg, epar, inpar, double);
      memcopy (periodarr + no_neg + inpar, pospar, no_pos, double);
      inpar += no_neg + no_pos;

      pararr = periodarr;
    }
  else
    pararr = epar;



  /* Check that the input parameter values lie in the parameter
     interval of the curve.  */

  tstart = *(pc->et);
  tend = *(pc->et + kn + kk - 1);
  for (k1 = 0; k1 < inpar; k1++)
    if (pararr[k1] < tstart || pararr[k1] > tend)
      goto err158;

  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix, and space for new knots
     to use in s1701.c */

  if ((salfa = newarray (kk, double)) == SISL_NULL)
    goto err101;
  if ((sp = newarray (kk, double)) == SISL_NULL)
    goto err101;

  /* Find the number of vertices in the new curve. */

  kn1 = kn + inpar;

  /* Allocating the new arrays to the new curve. */

  if ((st = newarray (kn1 + kk, double)) == SISL_NULL)
    goto err101;
  if ((scoef = new0array (kn1 * kdim, double)) == SISL_NULL)
    goto err101;

  /* Making the new knotvectors. */

  for (k2 = 0, k3 = 0, k1 = 0; k1 < inpar; k1++)
    {
      tpar = pararr[k1];
      for (; k2 < kn + kk; k2++)
	{
	  if (pc->et[k2] <= tpar)
	    st[k3++] = pc->et[k2];
	  else
	    break;
	}
      st[k3++] = tpar;
    }
  for (; k2 < kn + kk; k2++)
    st[k3++] = pc->et[k2];

  /* Updating the coefisientvector to the new curve.*/

  for (s1 = scoef, ki = 0, kmy = 0; ki < kn1; ki++)
    {
      /* Here we compute a new line with line number ki of
	 the knot inserten matrix. */

      while (pc->et[kmy + 1] <= st[ki])
	kmy++;
      s1701 (ki, kmy, kk, kn, &kpl, &kfi, &kla, st, pc->et, sp, salfa, &kstat);
      if (kstat)
	goto err153;

      /* Compute the kdim vertices with the same "index". */

      for (kj = 0; kj < kdim; kj++, s1++)
	for (*s1 = 0, kj1 = kfi, kj2 = kfi + kpl; kj1 <= kla; kj1++, kj2++)
	  {
	    ki2 = kj1 * kdim + kj;
	    *s1 += salfa[kj2] * coef[ki2];
	  }
    }

  /* Periodicity treatment -------------------------- */

  if (pc->cuopen == SISL_CRV_PERIODIC)
    {
      kn1 -= no_pos + no_neg;
      /* Allocating the new arrays to the new curve. */

      if ((st2 = newarray (kn1 + kk, double)) == SISL_NULL)
	goto err101;
      if ((scoef2 = new0array (kn1 * kdim, double)) == SISL_NULL)
	goto err101;
      memcopy (st2, st + no_neg, kn1 + kk, double);
      memcopy (scoef2, scoef + no_neg * kdim, kn1 * kdim, double);

      if (st)
	freearray (st);
      if (scoef)
	freearray (scoef);

      /* Allocating new curve-objects.*/

      if ((q1 = newCurve(kn1,kk,st2,scoef2,pc->ikind,pc->idim,2)) == SISL_NULL)
	goto err101;
    }
  else
    {
      /* Allocating new curve-objects.*/
      if ((q1 = newCurve (kn1, kk, st, scoef, pc->ikind, pc->idim, 2)) == SISL_NULL)
	goto err101;
    }

  q1->cuopen = pc->cuopen;

  /* Updating output. */

  *rcnew = q1;
  *jstat = 0;
  goto out;


  /* Error. Subrutine error. */

err153:*jstat = kstat;
  goto outfree;


  /* Error. No curve to treat.  */

err150:*jstat = -150;
  goto out;

  /* Error. Parameter value to insert is outside the curve. */

err158:*jstat = -158;
  goto out;

  /* Error. Allocation error, not enough memory.  */

err101:*jstat = -101;
  goto outfree;


outfree:
  if (q1)
    freeCurve (q1);
  else
    {
      if (st)
	freearray (st);
      if (scoef)
	freearray (scoef);
      if (st2)
	freearray (st2);
      if (scoef2)
	freearray (scoef2);
    }


  /* Free local used memory. */

out:if (salfa)
    freearray (salfa);
  if (sp)
    freearray (sp);
  if (periodarr)
    freearray (periodarr);
  if (pospar)
    freearray (pospar);
  if (negpar)
    freearray (negpar);

}
