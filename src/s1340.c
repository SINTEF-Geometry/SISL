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
 * $Id: s1340.c,v 1.5 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1340

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1340(SISLCurve *oldcurve, double eps[], int startfix, int endfix,
      double epsco, int itmax, SISLCurve **newcurve, double maxerr[],
      int *stat)
#else
void
s1340(oldcurve, eps, startfix, endfix, epsco, itmax, newcurve,
      maxerr, stat)
    SISLCurve 	*oldcurve;
    double	eps[];
    int		startfix;
    int		endfix;
    double	epsco;
    int		itmax;
    SISLCurve	**newcurve;
    double	maxerr[];
    int		*stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* Purpose: To remove as many knots as possible from the spline curve
*	   "oldcurve" without perturbing the curve more than a given
*	   tolerance. The tolerance is given by "eps" and the
* 	   approximation by "newcurve".
*
*
* Input:
*          oldcurve    - pointer to the original spline curve. Note that
*                        if the polynomial order of the curve is k, then
*                        the knot vector of oldcurve is assumed to have
*                        knots of multiplicity k at the beginning and end.
*
*	   eps	       - double array giving the desired absolute accuracy
*                        of the final approximation as compared to oldcurve.
*                        If oldcurve is a spline curve in a space of
*			 dimension dim, then eps must have length dim.
*                        Note that it is not relative, but absolute accuracy
*                        that is being used. This means that the difference
*                        in component i at any parameter value, between
*                        the given curve and the approximation, is to be
*                        less than eps[i]. Note that in such comparisons
*                        the same parametrization is used for both curves.
*
*	   startfix    - the number of derivatives to be kept fixed at
*			 the beginning of the knot interval.
*                        The (0 - (startfix-1)) derivatives will be kept fixed.
*                        If startfix <0, this routine will set it to 0.
*                        If startfix< the order of the curve, this routine
*                        will set it to the order.
*
*	   endfix      - the number of derivatives to be kept fixed at
*		         the end of the knot interval.
*                        The (0 - (endfix-1)) derivatives will be kept fixed.
*                        If endfix <0, this routine will set it to 0.
*                        If endfix< the order of the curve, this routine
*                        will set it to the order.
*
*	   epsco       - real number used to check equality between real
*                        numbers in some routines. Two numbers differing
*                        by a relative amount less than epsco will in some
*                        cases be considered equal. A suitable value is just
*			 above the unit roundoff of the machine. In IEEE
*                        double precision arithmetic 1.0E-15 is a reasonable
*                        choice.
*                        NB! The computations are not guaranteed to have
*                            relative accuracy less than epsco.
*
*	   itmax       - maximum number of iterations. The routine will
*                        follow an iterative procedure trying to remove
*                        more and more knots. The process will almost always
*                        stop after less than 10 iterations and it will often
*                        stop after less than 5 iterations. A suitable
*                        value for itmax is therefore usually in the region
*                        3-10.
*
*
*
*
* Output:
*          newcurve    - the spline approximation on the reduced
*                        knot vector.
*	   maxerr      - double array containing an upper bound for the
*			 pointwise error in each of the components of the
*		         spline approximation. The two curves oldcurve and
*			 newcurve are compared at the same parameter value,
*                        i.e., if oldcurve is f and newcurve is g, then
*			               |f(t)-g(t)| <= eps
*			 in each of the components.
*
*           stat       - status message
*                           > 0      : warning
*                           = 0      : ok
*                           < 0      : error
*
*
*
*
* Method:
*     The method is described in two papers listed as references below.
*     First, the knots of the input spline, s, are ranked according to their
*     relative importance in the representation of s (s1353). Using this
*     information, knots are removed from s and an approximation g is computed
*     on the reduced knot vector (s1354) with the error guaranteed to be
*     less than the tolerance in all components. In order to remove more knots,
*     the knots of g are ranked, and then knots are removed from g, resulting
*     in a spline h with even fewer knots. The knot vector of h is then
*     clearly a subsequence of the knot vector of s and s is approximated
*     on this knot vector (sh1365). If the error is less than the tolerance
*     then h is accepted as the new best approximation and stored in g.
*     If the error is greater than the tolerance, knots are removed directly
*     from s based on the ranking of g resulting in a new g. In either case a
*     new best approximation g is obtained, and we can continue the iteration.
*     The iterations stop when no more knots can be removed or all the interior
*     knots have been removed. (The names s, g, and h are not used in the
*     code.)
*
* NB! It is assumed that the input spline has knots of multiplicity k (the
*     polynomial order) at each end. This property of the knot vector is
*     used extensively by the lower level routines.
*
*
* References:
*     1. A Data-Reduction Strategy for Splines with Applications to the
*        Approximation of Functions and data, IMA J. of Num. Anal. 8 (1988),
*        pp. 185-208.
*
*     2. Knot Removal for Parametric B-spline Curves and Surfaces,
*        CAGD 4 (1987), pp. 217-230.
*
*
* Calls: s1353, s1354, sh1365, s6err
*
* Written by: Knut Moerken, University of Oslo, November 1992, based
*             on an earlier Fortran version by the same author.
* CHANGED BY: Paal Fugelli, SINTEF, 1994-07.  Initialized pointers (to SISL_NULL)
*      in 'ranking' to avoid potential memory leak when exiting through 'out'.
*      Removed several other memory leaks.
* CHANGED BY: Per OEyvind Hvidsten, SINTEF, 1994-11. Added a freeCurve
*      call before overwriting the *newarray pointer (thus removing a
*      memory leak.
*
*********************************************************************
*/
{
  char ready, big;              /* Boolean variables used for computing
                                   the stopping criteria.                */
  int k = oldcurve->ik;         /* Unwrapped version of the given curve. */
  int dim = oldcurve->idim;
  double *d = oldcurve->ecoef;
  int itcount=0;                /* Iteration counter.                    */
  int n1 = oldcurve->in;        /* n1 and n2 are used for storing the number
                                   coefficients in consecutive spline
                                   approximations. The situation n1=n2
				   signifies that no more knots can be
				   removed.                              */
  int n2;
  int i, mini, maxi, indx;      /* Various auxiliary integer variables.  */
  double *local_err = SISL_NULL;     /* Variables used for storing local error
                                   estimates.                            */
  double *l2err = SISL_NULL;
  double *lepsco = SISL_NULL;
  double *temp_err = SISL_NULL;
  SISLCurve *tempcurve = SISL_NULL;  /* Variables that are used for storing
                                   temporary curves.                     */
  SISLCurve  *helpcurve = SISL_NULL;
  rank_info ranking;            /* Variable used for holding ranking
                                   information (from s1353).             */
  int lstat=0;                  /* Local status variable.                */
  int pos=0;			/* Parameter to s6err.                  */
  SISLCurve *qc_kreg = SISL_NULL;    /* Non-periodic version of the input curve. */


  /* Initialize ranking ptrs in case of early exit through 'out' (PFU 05/07-94) */
  ranking.prio = SISL_NULL;
  ranking.groups = SISL_NULL;

  /* Initialize maxerr to zero. */

  for (i=0; i<dim; i++) maxerr[i] = 0.;

  /* Only interior knots may be removed so if n1==k we can stop
     straight away.                                             */

  if (n1 == k)
  {
    *newcurve = newCurve(n1, k, oldcurve->et, oldcurve->ecoef,
			 oldcurve->ikind, oldcurve->idim, 1);
    if (*newcurve == SISL_NULL)  goto err101;

    *stat = 0;
    goto out;
  }

  /* Make sure that the input curve is non-periodic.  */

  if (oldcurve->cuopen == SISL_CRV_PERIODIC)
  {
    make_cv_kreg(oldcurve, &qc_kreg, &lstat);
    if (lstat < 0) goto error;

    /* The input curve is closed and periodic. Make sure that the
       endpoints of the curve is still matching by fixing the
       position and the derivative in the endpoints of the curve.
       The change made to startfix and endfix is only locallly. */

    startfix = MAX(startfix, 2);
    endfix = MAX(endfix, 2);
  }
  else
    qc_kreg = oldcurve;

  /* Allocate space for some local arrays. */

  temp_err = newarray(dim, double);
  if (temp_err == SISL_NULL) goto err101;

  lepsco = newarray(dim, double);
  if (lepsco == SISL_NULL) goto err101;

  /* ranking is of type rank_info which is a struct described in s1353. */

  ranking.prio = newarray(MAX(n1-k,1), int);
  if (ranking.prio == SISL_NULL) goto err101;

  ranking.groups = newarray(MAX(n1-k,1), int);
  if (ranking.groups == SISL_NULL) goto err101;

  /* lespco is needed in s1354. In component i we first store the l1-norm
     of the component i of the B-spline coefficients of oldcurve. */

  for (i=0; i<dim; i++) lepsco[i] = fabs(d[i]);

  indx = 0;
  for (i=1; i<n1*dim; i++)
  {
    lepsco[indx++] += fabs(d[i]);
    if (indx == dim) indx = 0;
  }

  /* We can now compute the final value of lepsco which will be used for
     checking if two numbers are almost equal. */

  for (i=0; i<dim; i++)
    lepsco[i] = MIN(lepsco[i]*epsco/n1, eps[i]);

  /* This is where the knot removal process starts. */

  /* mini and maxi are lower and upper bounds on how many knots that can
     be removed. */

  mini = 0;
  maxi = n1 - k + 1;

  /* Start by ranking the knots of oldcurve. */

  s1353(qc_kreg, eps, &ranking, &lstat);
  if (lstat < 0) goto error;

  /* Based on the computed ranking, we remove as close to maxi knots
     from oldcurve as possible, but such that we can always compute an
     approximation to oldcurve on the reduced knot vector, with error less
     than eps.  newcurve will be created with icopy==1.  */

  s1354(qc_kreg, qc_kreg, &ranking, eps, lepsco, startfix, endfix,
	mini, maxi, newcurve, maxerr, &lstat);
  if (lstat < 0) goto error;

  /* The spline stored in newcurve is now an approximation to oldcurve
     with error less than eps. We will now iterate and try to remove knots
     from newcurve. The integers n1 and n2 will be the number of knots in
     the two most recent approximations. */

  n2 = (*newcurve)->in;

  /* Start the iterations. We have already done one iteration and we have
     not had any problems with the error getting too big. */

  itcount = 1;
  big = 0;

  /* If n1=n2 we are unable to remove any knots, and if n2=k there are
     no interior knots left so we can stop. Likewise if itmax was 1. */

  ready = (n1 == n2) || (n2 == k) || (itcount == itmax);
  while (!ready)
  {

    /* We start by ranking the knots of newcurve which is the approximation
       with the fewest knots that we have found so far. */

    s1353(*newcurve, eps, &ranking, &lstat);
    if (lstat < 0) goto error;

    if (!big)
    {

      /* If we have not had any problems with the error in approximation
         becoming too big we take the risk of removing as many knots as
         we can from newcurve and in this way determining an even shorter
	 knot vector.
	 First we iterate in s1354 to see how many knots we can remove
	 and store the approximation to newcurve in helpcurve (will be
	 created with icopy==1). */

      mini = 0; maxi = (*newcurve)->in - k + 1;
      s1354(*newcurve, *newcurve, &ranking, eps, lepsco, startfix,
	    endfix, mini, maxi, &helpcurve, temp_err, &lstat);
      if (lstat < 0) goto error;

      /* Now, we have no guarantee that helpcurve is closer to oldcurve
         than the tolerance so we compute a new approximation to oldcurve
         on the knot vector of helpcurve and store this in tempcurve (with
	 icopy==1).
	 Must make sure that local_err is free'ed if allocated from last
	 iteration. */

      if (local_err != SISL_NULL) freearray(local_err);

      sh1365(qc_kreg, helpcurve->et, k, helpcurve->in,
	     startfix, endfix, &tempcurve,
	     &local_err, &l2err, &lstat);
      if (lstat < 0) goto error;

      /* Don't need l2err or helpcurve anymore. */

      freearray(l2err);
      freeCurve(helpcurve);
      helpcurve = SISL_NULL;

      /* We must now check if the new tempcurve is within the tolerance. */
      i = 0;
      while (!big && i<dim)
      {
	big = local_err[i] > eps[i];
	i++;
      }
    }
    else

      /* If we have had problems with the error becoming greater than the
	 tolerance, we simply store newcurve in tempcurve and proceed. */

      tempcurve = *newcurve;

    if (big)
    {

      /* If at some stage we have had problems with the error becoming
	 larger than the tolerance, we do not get involved in the risky
	 approach of removing knots from newcurve. Instead we throw away
         tempcurve and remove knots directly from oldcurve, but remember
	 that the only knots left to remove are the knots of newcurve for
	 which we use the ranking computed above. The result is stored
	 in helpcurve and later transferred to newcurve. */

      mini= 0; maxi = tempcurve->in - k + 1;
      s1354(qc_kreg, *newcurve, &ranking, eps, lepsco,
	    startfix, endfix, mini, maxi, &helpcurve, maxerr, &lstat);
      if (lstat < 0) goto error;

      if (*newcurve != SISL_NULL && *newcurve != tempcurve)
      {
	freeCurve(*newcurve);
	*newcurve = SISL_NULL;
      }

      if (tempcurve != SISL_NULL)
      {
	freeCurve(tempcurve);
	tempcurve = SISL_NULL;
      }
      *newcurve = helpcurve;
      helpcurve = SISL_NULL;
    }
    else
    {

      /* Since we now know that the difference between oldcurve and tempcurve
	 is within the tolerance, we just have to store the error in maxerr,
	 throw away the old newcurve, and store tempcurve in newcurve. */

      for (i=0; i<dim; i++) maxerr[i] = local_err[i];

      if (*newcurve != SISL_NULL) freeCurve(*newcurve);
      *newcurve = tempcurve;
      tempcurve = SISL_NULL;
    }

    /* Now we must check if it is time to stop. */

    n1 = n2;
    n2 = (*newcurve)->in;
    ++itcount;
    ready = (n1 == n2) || (n2 == k) || (itcount == itmax);
  }

  /* Set periodicity flag.  */

  if (oldcurve->cuopen == SISL_CRV_CLOSED ||
      oldcurve->cuopen == SISL_CRV_PERIODIC)
    (*newcurve)->cuopen = SISL_CRV_CLOSED;

  /* Success */

  *stat = 0;
  goto out;

  /* Error in allocation of memory. */

err101:
  *stat = -101;
  s6err("s1340", *stat, pos);
  goto out;

  /* Error in lower level routine. */

error:
  *stat = lstat;
  s6err("s1340", *stat, pos);
  goto out;

  /* Clear up and free memory before exit. */

out:
  if (temp_err != SISL_NULL) freearray(temp_err);
  if (local_err != SISL_NULL) freearray(local_err);
  if (l2err != SISL_NULL) freearray(l2err);
  if (lepsco != SISL_NULL) freearray(lepsco);
  if (ranking.prio != SISL_NULL) freearray(ranking.prio);
  if (ranking.groups != SISL_NULL) freearray(ranking.groups);
  if (qc_kreg != SISL_NULL && qc_kreg != oldcurve) freeCurve(qc_kreg);
  if (helpcurve != SISL_NULL && helpcurve != (*newcurve)) freeCurve(helpcurve);
  if (tempcurve != SISL_NULL && tempcurve != (*newcurve)) freeCurve(tempcurve);

  return;
}
