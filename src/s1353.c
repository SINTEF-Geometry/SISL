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
 * $Id: s1353.c,v 1.3 2006-05-02 15:06:56 sbr Exp $
 *
 */

#define S1353

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1353(SISLCurve *curve,double eps[],rank_info *ranking,int *stat)
#else
void s1353(curve, eps, ranking, stat)
     SISLCurve	*curve;
     double	eps[];
     rank_info	*ranking;
     int *stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* Purpose: To compute a ranking of the knots of curve, suitable for use
*          during knot removal.
*
*
*
* Input:
*          curve       - pointer to the curve whose knots are to be ranked.
*
*	   eps         - double array containing the tolerance in the
*			 different components to be used during knot removal.
*
*
*
* Output:  ranking     - a pointer to a rank_info object containing the
*                        result of the ranking computations. The struct
*			 rank_info contains four variables listed below
*                        (t, k, and n refer to the knot vector of curve,
*			 its polynomial order and the number of B-spline
*                        coefficients):
*	 ranking->prio - integer array of dimension (n-k) containing the
*			 ranking of the interior knots. The knots are listed
* 			 in order of increasing ranking number, cf. the
*			 second reference, with knots with the same ranking
*		         number listed in the order in which they occur
*		         in t. To differentiate between the different
*                        ranking numbers the array ranking->groups is used.
*      ranking->groups - integer array of dimension (n-k) used to partition
*	                 ranking->prio into blocks of knots with the same
*                        ranking number. ranking->groups[0] contains the
*			 number of knots with the smallest ranking number,
*                        ranking->groups[1] contains the number of knots with
*                        the two smallest ranking numbers etc.. Only the first
*                        ranking->antgr elements of ranking->groups are used.
*       ranking->antgr - the number of distinct ranking numbers among the
*                        interior knots of t.
*      ranking->antrem - an estimate of the number of knots that can be
*                        removed from curve.
*
*           stat       - status message
*                           > 0      : warning
*                           = 0      : ok
*                           < 0      : error
*
*
*
* Method:
*      The ranking strategy is described in the two references.
*      somewhat simplified, the computation of the ranking number of knot
*      number i consists of computing the error in the best approximation to
*      the input spline from a subspace where this knot has been removed,
*      and then finding in which of the intervals [0,eps/2], (eps/2,eps],
*      (2*eps,4*eps],... , the error lies. This interval number is then
*      the ranking number of the knot.
*      The best approximation is computed in a certain discrete max-norm,
*      and this leads to well conditioned linear systems.
*
*      NB! It is assumed that the first and lasty knots occur with
*      multiplicity k.
*
*      NB! The code is unreadable without close study of the two references.
*
*
* References:
*      1. A Data-Reduction Strategy for Splines with Applications to the
*         Approximation of Functions and data, IMA J. of Num. Anal. 8 (1988),
*         pp. 185-208.
*
*      2. Knot Removal for Parametric B-spline Curves and Surfaces,
*         CAGD 4 (1987), pp. 217-230.
*
*
* Written by: Knut Moerken, University of Oslo, July 1992, based on an
*              earlier Fortran version by the same author.
*
*********************************************************************
*/
{
  int k = curve->ik;            /* The basic parameters of curve will be */
  int dim = curve->idim;        /* used extensively and are therefore    */
  int n = curve->in;            /* extracted and stored in local         */
  double *et = curve->et;       /* variables.                            */
  double *ecoef = curve->ecoef;
  int *prio = ranking->prio;    /* Local versions of the ranking arrays. */
  int *groups = ranking->groups;
  int antgr;                    /* The two integers of ranking.          */
  int antrem;
  int pos=0;                    /* Integer used by s6err.                */
                                /* The routine requires a large number of
				   integers and floats, see the code.    */
  int i, ii, i2, inxt, mult, dm, j1, j2, nkt, dih, dj, dj1, hmult,
      p, dp, di2, j, count, di, iw, ph, gh;
  double h, x, pm1h, pm1, muj1, muj, mm, h1, h2, denom, ln2inv;
  double *c, *hh, *mu, *w, *e, *s;

  /* We start by allocating space for local arrays. */

  /* c will contain the right-hand-sides in the linear systems that
     are going to be solved. */

  c = newarray(dim*(k+1), double);
  if (c == SISL_NULL) goto err101;

  /* hh will contain the last column of the coefficient matrix during
     elimination. */

  hh = newarray(k+1, double);
  if (hh == SISL_NULL) goto err101;

  /* mu will contain one of the two diagonals of the coefficient matrix. */

  mu = newarray(k+1, double);
  if (mu == SISL_NULL) goto err101;

  /* w will contain the weights of the interior knots (of type double). */

  w = newarray(dim*(n-k), double);
  if (w == SISL_NULL) goto err101;

  /* e and s are auxialiary arrays. */

  e = newarray(dim, double);
  if (e == SISL_NULL) goto err101;

  s = newarray(dim, double);
  if (s == SISL_NULL) goto err101;

  /* Start the computation of the weights. We first compute real (not integer)
     weights, w. The weight of a knot is computed by removing it and computing
     a best approximation without the knot. The weight is the error in this
     approximation. The norm used is the max-norm of the B-spline coefficients
     on the original knot vector et. */

  /* The integer i points to the previous knot that was removed, while
     dih points to the location in w where the next weight is to be stored.
     nkt gives the number of interior knots and therefore the number of
     weights to be computed (the number of weights is nkt*dim).
     We clearly have to loop until i=n-1 (remember that we assumed the first
     and last knots to occur with multiplicity k). */

  i = k - 1; dih = 0; nkt = n - k;

  while (i < n-1)
    {

      /* The knot that we want to compute a weight for now is x=et[i+1].
	 inxt will point to the smallest knot distinct from x. */

      inxt = i + 1; x = et[inxt];
      while (x == et[inxt+1]) ++inxt;

      /* mult is the multiplicity of x, while dm is the dimension of
	 the linear system we are going to solve, see reference 1. */

      mult = inxt - i; dm = k - mult + 2;

      /* Set up the linear system of equations. The first dm-1 unknowns
         are the coefficients of the spline approximation that are different
	 from the coefficients of curve while the last unknown is the
	 error in the approximation. */

      /* pm1h is the first entry in last column of the coefficient matrix. */

	       pm1h = 1.0;
      if ( fmod((float) dm, (float) 2.0) == 0 )  pm1h = -1.0;

      /* The integer j is used to step through the equations during
	 elimination, while muj is used to hold the mu of equation j
	 in ref. 1. The array mu (and muj1) will be used for storing
	 the inverse of muj, while mm is the factor that the previous
	 (or next0 equation is multiplied with to eliminate one variable.
	 hh is used for keeping track of the last column of the coefficient
	 matrix during elimination, while j1 and j2 are pointers to knots
	 in et that are needed for computing muj. dj1 is a pointer to the
	 B-spline coefficient that occurs on the right-hand-side
	 of the equation (really dim right-hand-sides), while dj
         is a pointer to the right-hand-side after elimination
	 (i2 is always dim behind (or ahead of) dj.) */

      muj1 = 1.0; mu[0] = 1.0;
      j1 = inxt - k + 1; dj1 = dim*j1;
      j2 = inxt + 1;
      hh[0] = pm1 = pm1h;

      di2 = (j1-1)*dim;

      for (ii = 0; ii < dim; ++ii, ++di2) c[ii] = ecoef[di2];

      j = 1; dj = dim;
      i2 = 0; di2 = 0;

      /* Eliminate from the top until muj1 becomes larger than 2
	 (see ref. 1). */

      while (muj1 <= 2 && j < dm - 1)
	{
	  h = et[j1];

	  /* The two conditions tested here correspond to a singular system
	     which can only occur with a buggy knot vector. */

	  if (h == et[j2] || x == h) goto err112;
	  muj = (x-h) / (et[j2]-h);
	  mm = (1-muj) * muj1;
	  mu[j] = muj1 = 1 / muj;
	  pm1 = -pm1;
	  hh[j] = pm1 - mm*hh[i2];
	  for (ii=0; ii<dim; ++ii, ++dj, ++dj1, ++di2)
	    c[dj] = ecoef[dj1] - mm*c[di2];

	  i2 = j; ++j; ++j1; ++j2;
	}

      /* Now we remember how far down we managed to eliminate in p and prepare
	 to eliminate from the last equation to p+1. The variables have the
	 same interpretations as before except for the obvious changes since
	 we are eliminating upwards from the bottom and not dowmwards from
	 the top. */

      p = i2;

      /* Initialize variables as above. */

      j1 = i; dj1 = dim*(j1+1)-1;
      j2 = j1 + k;
      i2 = j1 + 1; di2 = dim*(i2+1) - 1;
      pm1 = 1;
      hh[dm-1] = 1.0;
      dj = dim*dm - 1;

      for (ii=0; ii<dim; ++ii, --dj, --di2) c[dj] = ecoef[di2];

      muj1 = mu[dm-1] = 1.0;
      i2 = dm - 1; di2 = dim*dm - 1;

      /* Eliminate. */

      for (j=dm-2; j>p; --j)
	{
	  h = et[j2];
	  if (h == x || h == et[j1]) goto err112;
	  muj = (h-x) / (h-et[j1]);
	  mm = (1-muj)*muj1;
	  muj1 = 1/muj;
	  mu[j] = muj1;
	  pm1 = -pm1;
	  hh[j] = pm1 - mm*hh[i2];
	  for (ii=0; ii<dim; ++ii, --dj, --dj1, --di2)
	    c[dj] = ecoef[dj1] - mm*c[di2];

	  i2 = j; --j1; --j2;
	}

      /* Now we are left with two equations (no. p and p+1) in two unknowns
	 (one B-spline coefficient and the error in the approximation) which
	 we solve the standard way. The weight is then obtained by taking
	 the absolute value of the error. */

      h1 = mu[i2]; h2 = mu[p];
      denom = 1/(h1*hh[i2]-h2*hh[p]);
      dp = dim*p; di2 = dim*i2;
      for (ii=0; ii<dim; ++ii, ++di2, ++dp, ++dih)
	{
	  h = e[ii] = (h1*c[di2] - h2*c[dp])*denom;
	  w[dih] = fabs(h);
	}

      /* We now determine the weights of the remaining knots at x. */

      count = 1;

      for (hmult=mult-1; hmult>0; --hmult)
	{

	  /* First do the back substitution of the system used to determine
             the previous weight. In this way the nontrivial B-spline
	     coefficients of the spline that best approximates the given
	     spline curve but without the knot x, is determined.
             The weight of the next knot at x is obtained by solving a
             linear system of the same sort as the one above, but based on
	     this new spline instead. In this way a weight for the new
	     spline is found; the weight of this knot with respect to the
	     original spline is then set equal to this weight plus the
	     weight of the previous spline (at this knot). */

	  /* Substitute to the end. */

	  i2 = i + count; dj = dim*(p+1);
	  for (j=p+1; j<dm; ++j)
	    {
	      h1 = hh[j]; h2 = mu[j];
	      for (ii=0; ii<dim; ++ii, ++dj)
		c[dj] = (c[dj] - h1*e[ii])*h2;
	    }

	  /* Substitute to the beginning. */

	  dj1 = dim*(p+1) - 1; dj = dj1 - dim;
	  for (j=p-1; j>=0; --j)
	    {
	      h1 = hh[j]; h2 = mu[j];
	      for (ii=dim-1; ii>-1; --ii, --dj, --dj1)
		{
		  c[dj1] = (c[dj] - h1*e[ii])*h2;
		}
	    }

	  /* Set up the new system of equations. The dimension will be
	     one greater than the previous one and the system is based
	     on the spline that was just computed so that the multiplicity
	     of x has been reduced by one. */

	  /* The variables have the same interpretations as above. */
	  ++count; ++dm;
	  pm1h = -pm1h;
	  pm1 = hh[0] = pm1h;
	  j1 = i + hmult - k + 1;
	  j2 = inxt + 1;
	  i2 = j1 - 1; di2 = dim*i2;
	  muj1 = 1.0;

	  for (ii=0; ii<dim; ++ii, ++di2) c[ii] = ecoef[di2];


	  j = 1; dj = dim;
	  i2 = 0; di2 = 0;

	  /* Eliminate unknowns until muj1 becomes greater than 2. */

	  while (muj1 <= 2 && j < dm - 1)
	    {
	      h = et[j1];
	      if (h == et[j2] || x == h) goto err112;
	      muj = (x-h) / (et[j2]-h);
	      mm = (1-muj) * muj1;
	      mu[j] = muj1 = 1 / muj;
	      pm1 = -pm1;
	      hh[j] = pm1 - mm*hh[i2];
	      for (ii=0; ii<dim; ++ii, ++dj, ++di2)
		c[dj] = c[dj] - mm*c[di2];

	      i2 = j; ++j; ++j1; ++j2;

	    }

	  /* Remember where we have are. */

	  p = i2;

	  /* Prepare elimination from the end up to p. */

	  j1 = i; dj1 = dim*(j1+1) - 1;
	  j2 = j1 + k + count - 1;
	  i2 = i + count; di2 = dim*(i2+1) - 1;
	  pm1 = 1;
	  hh[dm-1] = 1.0;
	  dj = dim*dm - 1;

	  for (ii=0; ii<dim; ++ii, --dj, --di2) c[dj] = ecoef[di2];


	  muj1 = mu[dm-1] = 1.0;
	  i2 = dm - 1; di2 = dim*dm - 1;

	  /* Eliminate from the end. */
	  for (j=dm-2; j>p; --j)
	    {
	      h = et[j2];
	      if (h == x || h == et[j1]) goto err112;
	      muj = (h-x) / (h-et[j1]);
	      mm = (1-muj)*muj1;
	      muj1 = 1/muj;
	      mu[j] = muj1;
	      pm1 = -pm1;
	      hh[j] = pm1 - mm*hh[i2];
	      for (ii=0; ii<dim; ++ii, --dj, --di2)
		c[dj] = c[dj] - mm*c[di2];

	      i2 = j; --j1; --j2;
	    }

	  /* Solve the two equations in two unknowns and determine w
	     by adding the absolute value of the current error to the
	     previous weight. */

	  h1 = mu[i2]; h2 = mu[p];
	  denom = 1/(h1*hh[i2]-h2*hh[p]);
	  dp = dim*p; di2 = dim*i2;
	  for (ii=0; ii<dim; ++ii, ++di2, ++dp, ++dih)
	    {
	      h = e[ii] = (h1*c[di2] - h2*c[dp])*denom;
	      w[dih] = w[dih-dim] + fabs(h);
	    }
	}

      /* Go to the next distinct knot. */

      i = inxt;

    }

  /* We now have the information that is needed for computing the
     ranking numbers, see refs. 1 and 2. */

  /* We start by computing the inverse of eps/2 (careful if eps is small). */

  for (ii=0; ii<dim; ++ii)
    s[ii] = 2.0/MAX(eps[ii],1.0e-300);

  /* We need the constant 1/logarithm(2). */

  ln2inv = 1.0 / log((double) 2.0);
  di = 0;

  /* Since all the weights are nonnegative, each of them must fall in one
     of the intervals [0,eps/2), [eps/2, eps), [eps, 2eps), [2eps, 4eps), ...
     If we number these intervals from 0, then the ranking number of a knot
     is the number of the interval that contains its weight. For a parametic
     curve we choose the ranking number to be the largest of the ranking
     numbers of the dim components. The ranking numbers are stored in
     the array groups. */

  for (i=0; i<nkt; ++i)
    {
      iw = 0;
      for (ii=0; ii<dim; ++ii, ++di)
	{
	  h = w[di];
	  if (h > 0.5*eps[ii])
	    iw = MAX(iw, (int)floor(log(h*s[ii])*ln2inv+1));
	}
      groups[i] = iw;
    }

  /* We next sort the ranking numbers in increasing order. Since each ranking
     number belongs to a specific knot, this induces an ordering of the knots.
     This ordering is stored in prio. Note that if two knots have the same
     ranking they are listed in prio in the order that they occur in et. */

  for (i=0; i<nkt; ++i) prio[i]=i;

  for (i=1; i<nkt; ++i)
    {
      gh = groups[i]; ph = prio[i];

      for (j=i-1; gh < groups[j];)
	{
	  groups[j+1] = groups[j];
	  prio[j+1] = prio[j];
	  --j;
	  if (j <0) goto out;
	}

    out:
      groups[j+1] = gh; prio[j+1] = ph;
    }

  /* We next prepare for estimating how many knots that can be removed (very
     heuristic and not important for the perfomance of the knot removal. */

  if (groups[0] == 0) antrem = 0;
    else
      antrem = -1;

  /* The array groups should end up with the cumulative number of knots in
     each interval. */

  j = 0;
  for (i=0; i<nkt-1; ++i)
    if (groups[i] < groups[i+1])
      {
	groups[j] = i+1; ++j;
      }

  antgr = ++j;
  groups[antgr-1] = n - k;

  /* More guessing about antrem. */

  if (antrem == 0)
    {
      if (antgr > 1)
	antrem = (int)(0.75*groups[1]/(antgr-1));
      else
	antrem = (int)(0.9*groups[0]);
    }
  else
    antrem = 1;

  antrem = MAX(MIN(1, nkt), antrem);

  /* Finish of ranking. */

  ranking->antgr = antgr;
  ranking->antrem = antrem;



  /* Free memory. */

  freearray(hh);
  freearray(c);
  freearray(mu);
  freearray(w);
  freearray(e);
  freearray(s);

  /* Successful computations. */

  return;

  /* Error in allocation of memory. */

 err101:
  *stat = -101;
  s6err("s1353", *stat, pos);
  goto err;

  /* Error in knot vector (singular system of equations). */

 err112:
  *stat = -112;
  s6err("s1353", *stat, pos);
  goto err;

  /* Clean up. */

 err:
  if (c != SISL_NULL) freearray(c);
  if (hh != SISL_NULL) freearray(hh);
  if (mu != SISL_NULL) freearray(mu);
  if (w != SISL_NULL) freearray(w);
  if (e != SISL_NULL) freearray(e);
  if (s != SISL_NULL) freearray(s);

  return;

}
