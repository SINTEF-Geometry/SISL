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
 * $Id: s1017.c,v 1.4 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1017

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1017 (SISLCurve * pc, SISLCurve ** rc, double apar, int *jstat)
#else
void
s1017 (pc, rc, apar, jstat)
     int *jstat;
     double apar;
     SISLCurve *pc, **rc;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Insert a given knot into the description of
*              a B-spline curve.
* NOTE       : When the curve is periodic, the input parameter value
*              must lie in the HALFOPEN [et[kk-1], et[kn), the function
*              will automatically update the extra knots and
*              coeffisients.
*              rcnew->in is still eq to pc->in + 1!
*
*
* INPUT      : pc        - SISLCurve to be refined.
*              apar      - Parameter values of knot to be s1017ed.
*
*
*
* OUTPUT     : rc        - The new, refined curve.
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
*              S1701.C   - Making the knot-s1017en-transformation matrix.
*              s1017knots in periodic case.
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* CHANGED BY : Ulf J. Krystad, SI, 92-01
*              Treatment of periodic curves.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Undefined
*              symbol 's1017knots()' changed to 's1018()'. Closed ('cuopen'
*              flag=0) and 'ikind' flag are now handled correctly. Fixed memory
*              problem for rationals.  Cleaned up a bit.
*
**********************************************************************/
{

  int kstat;			/* Local status variable.                     */
  int kmy;			/* An index to the knot-vector.               */
  int kpl, kfi, kla;		/* To posisjon elements in trans.-matrix.     */
  int kk = pc->ik;		/* Order of the input curve.                  */
  int kn = pc->in;		/* Number of the vertices in input curve.     */
  int kdim = pc->idim;		/* Dimension of the space in whice curve lies.*/
  int kn1;			/* Number of vertices in the new curve.       */
  int kch;			/* First vertice to be changes.               */
  int knum;			/* Number of knots less and equal than
                                   the intersection point.                    */
  int ki;			/* Control variable in loop.                  */
  int kj, kj1, kj2;		/* Control variable in loop.                  */
  double *s1;           	/* Pointers used in loop.                     */
  double *st = SISL_NULL;		/* The first new knot-vector.                 */
  double *salfa = SISL_NULL;		/* A line of the trans.-matrix.               */
  double *scoef = SISL_NULL;		/* The first new vertice.                     */
  double *sp;			/* Help array for use in s1701.               */
  double *ecoef;                /* Pointer to input curve ecoef/rcoef vector. */
  SISLCurve *qc = SISL_NULL;		/* Pointer to new curve-object.               */



  /* Make sure the returned pointer is valid. */

  *rc = SISL_NULL;


  /* Check that we have a curve. */

  if (!pc)
    goto err150;


  /* Find the number of vertices in the new curve. */

  kn1 = kn + 1;

  if ( kn1 <= 0 )
    goto err150;



  /* Periodicity treatment -------------------------- */
  if (pc->cuopen == SISL_CRV_PERIODIC)
    {
      s1018(pc, &apar, 1, rc, &kstat);
      if (kstat < 0)
	goto err153;
      goto out;
    }



  /* Check that the intersection point is an interior point. */

  if (apar < *(pc->et) || apar > *(pc->et + kn + kk - 1))
    goto err158;


  /* Check if the curve is rational.  */

  if (pc->ikind == 2 || pc->ikind == 4)
  {
    ecoef = pc->rcoef;
    kdim++;
  }
  else
    ecoef = pc->ecoef;


  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix and for a help array of
     kk elements.*/

  if ((salfa = newarray (2 * kk, double)) == SISL_NULL)
    goto err101;

  sp = salfa + kk;

  /* Find the number of the knots which is smaller or like
     the intersection point.*/

  s1 = pc->et;

  if (apar > pc->et[0] && apar < pc->et[kn + kk - 1])
    {
      /* Using binary search*/
      kj1 = 0;
      kj2 = kk + kn - 1;
      knum = (kj1 + kj2) / 2;
      while (knum != kj1)
	{
	  if (s1[knum] < apar)
	    kj1 = knum;
	  else
	    kj2 = knum;
	  knum = (kj1 + kj2) / 2;
	}
      knum++;			/* The smaller knots. */

      while (s1[knum] == apar)
	/* The knots that are equal to the intersection point. */
	knum++;
    }
  else if (apar == pc->et[0])
    {
      knum = 0;
      while (s1[knum] == apar)
	/* The knots that are equal to the intersection point. */
	knum++;
    }
  else if (apar == pc->et[kn + kk - 1])
    {
      knum = 0;
      while (s1[knum] < apar)
	/* The knots that are less than or equal to the intersection point. */
	knum++;
    }



  /* Allocating the new arrays to the new curve. */

  if ((scoef = newarray (kn1 * kdim, double)) == SISL_NULL)
    goto err101;
  if ((st = newarray (kn1 + kk, double)) == SISL_NULL)
    goto err101;


  /* Copying the knotvectors, all but the intersection point from
     the old curve to the new curves */

  memcopy (st, pc->et, knum, double);
  st[knum] = apar;
  if (knum < kn + kk)
    memcopy (st + knum + 1, pc->et + knum, kn + kk - knum, double);


  /* Copying the coefisientvector to the new curve.  (Here 'ecoef' points
     to 'pc->rcoef' or 'pc->ecoef' depening on if 'pc' is rational or not). */

  kch = knum - kk + 1;

  if (kch > 0)
    memcopy (scoef, ecoef, kdim * kch, double);
  if (knum < kn1)
    memcopy (scoef + kdim * knum, ecoef + kdim * (knum - 1),
	     kdim * (kn1 - knum), double);


  /* Updating the coefisientvectors the new curve.*/

  /* Updating the first curve. */

  for (ki = max (0, kch), kmy = 0, s1 = scoef + ki * kdim; ki < min (knum + 1, kn1); ki++)
    {
      /* Initialising:
           ki = kch,        Index of the vertices we are going to
                             change. Starting with kch, but if
                             kch is negativ we start at zero.
           s1=scoef1+ki*kdim,Pointer at the first vertice to
                             change. */


      /* Using the Oslo-algorithm to make a transformation-vector
         from the old vertices to one new vertice. */

      while (kmy < kn + kk && pc->et[kmy] <= st[ki])
	kmy++;

      s1701 (ki, kmy - 1, kk, kn, &kpl, &kfi, &kla, st, pc->et, sp, salfa, &kstat);
      if (kstat)
	goto err153;


      /* Compute the kdim vertices with the same "index". */

      for (kj = 0; kj < kdim; kj++, s1++)
	for (*s1 = 0, kj1 = kfi, kj2 = kfi + kpl; kj1 <= kla; kj1++, kj2++)
	  *s1 += salfa[kj2] * ecoef[kj1 * kdim + kj];
    }



  /* Allocating new curve-objects.*/


  if ((qc = newCurve (kn1, kk, st, scoef, pc->ikind, pc->idim, 2)) == SISL_NULL)
    goto err101;


  /* Updating output. */

  qc->cuopen = pc->cuopen;  /* Ok since only internal knots can be inserted. */

  *rc = qc;

  *jstat = 0;
  goto out;


  /* Error. Error in low level routine. */

err153:
  *jstat = kstat;
  goto outfree;


  /* Error. No curve to s1017 a new knot.  */

err150:
  *jstat = -150;
  goto out;


  /* Error. The intersection-point is outside the curve.  */

err158:
  *jstat = -158;
  goto out;


  /* Error. Allocation error, not enough memory.  */

err101:
  *jstat = -101;
  goto outfree;


outfree:
  if (qc)
    freeCurve (qc);
  else
    {
      if (st)
	freearray (st);
      if (scoef)
	freearray (scoef);
    }

  /* Free local used memory. */

out:if (salfa)
    freearray (salfa);
}
