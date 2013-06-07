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


#define S1950

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1950(SISLCurve *oldcurve,SISLCurve *rankcurve,rank_info *ranking,
	   double eps[],double epsco[],int startfix,int endfix,int *jncont,
	   int mini,int maxi,SISLCurve **newcurve,double maxerr[],int *stat)
#else
void s1950(oldcurve, rankcurve, ranking, eps, epsco,
	   startfix, endfix, jncont, mini, maxi, newcurve, maxerr, stat)
     SISLCurve	*oldcurve;
     SISLCurve	*rankcurve;
     rank_info	*ranking;
     double	eps[];
     double	epsco[];
     int	startfix;
     int	endfix;
     int        *jncont;
     int	mini;
     int	maxi;
     SISLCurve	**newcurve;
     double	maxerr[];
     int	*stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* Purpose: To remove as many knots as possible from oldcurve
*          without perturbing this spline more than eps, using the ranking
*          information in rank_info which was obtained from the spline
*          rankcurve; and compute an approximation to oldcurve on the
*          knot vector of rankcurve. The knot vector of rankcurve must be
*          a subsequence of the knot vector of oldcurve and the knots are
*          always removed from rankcurve so that the final reduced knot
*          vector will always be a subsequence of both the knot vectors of
*          oldcurve and rankcurve. The final approximation is given in
*	   newcurve.
*
*
*
*
* Input:
*          oldcurve    - pointer to the original spline from which the knots
*			 are to be removed.
*
*          rankcurve   - pointer to the curve that the ranking information
*                        is based on. It is assumed the knot vector of
*                        rankcurve is a subsequence of the knot vector of
*			 oldcurve.
*
*	   ranking     - a pointer to a rank_info object containing the
*                        result of the ranking computations. The struct
*                        rank_info contains four variables listed below
*                        (t, k, and n refer to the knot vector of rankcurve,
*                        its polynomial order and the number of B-spline
*                        coefficients):
*        ranking->prio - integer array of dimension (n-k) containing the
*                        ranking of the interior knots. The knots are listed
*                        in order of increasing ranking number, cf. the
*                        second reference, with knots with the same ranking
*                        number listed in the order in which they occur
*                        in t. To differentiate between the different
*                        ranking numbers the array ranking->groups is used.
*      ranking->groups - integer array of dimension (n-k) used to partition
*                        ranking->prio into blocks of knots with the same
*                        ranking number. ranking->groups[0] contains the
*                        number of knots with the smallest ranking number,
*                        ranking->groups[1] contains the number of knots with
*                        the two smallest ranking numbers etc.. Only the first
*                        ranking->antgr elements of ranking->groups are used.
*       ranking->antgr - the number of distinct ranking numbers among the
*                        interior knots of t.
*      ranking->antrem - an estimate of the number of knots that can be
*                        removed from curve.
*
*          eps         - double array containing the tolerance. The
*			 approximation must at each parameter value deviate
*			 less than eps from oldcurve (in each component).
*
*	   epsco       - double array of containing a tolerance which
*                        is used when the number of coefficients in the
*                        approximating spline becomes very small. More
*                        specifically, if so many knots are removed that
*                        the corresponding spline approximation has fewer
*                        knots than startfix+endfix (the total number of
*                        constraints), then it is required that the error
*                        is less than epsco in each component.
*
*          startfix    - the number of derivatives to be kept fixed at the
*			 left end of the parameter interval.
*
*          endfix      - the number of derivatives to be kept fixed at the
*			 right end of the parameter interval.
*
*          jncont      - Number of continuity requirements across the seem.
*
*          mini        - a lower bound on how many knots that can be removed.
*
*          maxi        - an upper bound on how many knots that can be removed.
*
*
* Output:
*          jncont      - Number of continuity requirements across the seem.
*
*          newcurve    - pointer to the spline approximation on the reduced
*			 knot vector.
*          maxerr      - double array containing an upper bound for the
*                        pointwise error in each of the components of the
*                        spline approximation. The two curves oldcurve and
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
* Method:
*     The routine determines the maximum number of knots that can be removed
*     by binary search starting by removing ranking->antrem (antrem below)
*     knots. In other words, it first tries to remove antrem knots from
*     oldcurve. If this gives an error less than the tolerance, it tries
*     to remove (antrem+maxi)/2 knots, while if the error is greater than
*     the tolerance, the routine will try to remove (mini+antrem)/2 of the
*     knots. This is continued until the number of knots that can be removed
*     has been determined exactly, together with the corresponding
*     approximation. The acceptable approximation with the fewest knots is
*     kept in the variables igtau,igc,ign, and this is
*     initialized to the spline given by etprio,edprio,imprio.
*
*     During each iteration it must be determined which knots to remove,
*     given the number of knots to remove. The basis for this choice is
*     the ranking information (below groups and prio denote the arrays
*     ranking->groups and ranking->prio and t and k denote the knot vector
*     and order of rankcurve). Suppose we are to remove r knots.
*     First the smallest integer i satisying groups(i) >= r is determined.
*     In prio all the interior knots of t are listed (or rather their
*     index in t, t[k-1+5] has index 5), and since we are to remove
*     all knots in the first ranking groups, we remove the first groups[i-1]
*     knots of prio (assuming for simplicity that groups[-1]=0).
*     In the last group we reach into (which consists of knots no.
*     prio[groups[i-1]+1], prio[groups[i-1]+2],...,prio[groups[i]])
*     knots are removed uniformly on index, meaning that if we are to remove
*     half the knots in this group, we pick every other knot as they are
*     listed in prio, if we are to remove 1/3 of the knots, we remove every
*     three knots and so on. In general, it is difficult to say what
*     uniformly on subscripts should mean (how do you choose one knot from
*     six in a symmetric way) but the approach used here seems to work
*     reasonably well (see the program text below).
*     When it is known which knots to remove, the corresponding knot vector
*     can be constructed and an approximation computed.
*     the error is then checked as indicated above, and the iterations proceed.
*
*     For more information cf. the references:
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
* Calls: s1943, s6err
*
* Written by : Knut Moerken, University of Oslo, July 1992, based on an
*              earlier Fortran version.
* Changed by: Paal Fugelli, SINTEF, 1994-07.
*             Changed to fix several major memory leakage problems.
* Renamed and changed by: Vibeke Skytt, SINTEF Oslo, 01.95. Introduced
*                                                treatment of periodicity.
*
*********************************************************************
*/
{
  int k = oldcurve->ik;           /* Make some parameters more easily */
  int dim = oldcurve->idim;       /* accessible. */
  int mprio = rankcurve->in;
  int *prio = ranking->prio;
  int *group = ranking->groups;
  int antrem = ranking->antrem;
  int antgr = ranking->antgr;
  char *del_array;                /* Boolean array that is used for
                                     marking what knots to include.   */
  char big, bigco;                /* Boolean variables that are used
                                     for checking whether the error
				     is small enough.                 */
  int lstat=0;			  /* Local status variable.           */
  int pos=0;			  /* Parameter to s6err.              */
  int nlim = MAX(k, startfix+endfix); /* The minimum number of B-spline
					 coefficients that should be
					 left in newcurve.            */

				/* For the use of the other variables,
				   see the code below.              */
  int i, start, stop, indx, count, r, p, hn;
  int kncont = 0;
  SISLCurve *hcurve = SISL_NULL;
  double h;
  double *local_err = SISL_NULL, *l2_err = SISL_NULL, *ltau = SISL_NULL;
  double tend_domain = oldcurve->et[oldcurve->in];

  /* Allocate memory for the local arrays. */

  del_array = newarray(mprio, char);
  if (del_array == SISL_NULL) goto err101;
  if ((local_err = newarray(dim, DOUBLE)) == SISL_NULL) goto err101;
  if ((l2_err = newarray(dim, DOUBLE)) == SISL_NULL) goto err101;

  /* In case we do not enter the while loop at all we must give newcurve
     a value. */

  *newcurve = newCurve(mprio, k, rankcurve->et, rankcurve->ecoef, 1, dim, 1);
  if (newcurve == SISL_NULL) goto err101;

  /* Iterate by binary search until the lower and upper bound on how many
     knots to include are essentially equal. */

  while (mini+1 < maxi)
  {

    /* To start with none of the knots of rankspline are marked
       for removal. */

    for (i=0; i<mprio; i++) del_array[i] = 0;
    
    /* Initialize number of knots to remove at the seem. */
    
    kncont = 0;

    /* We then have to find out which knots to remove. We remove knots
       group by group, with start pointing (into prio) to the first knot
       of the current group and stop to the one following the last of this
       group. */

    start = 0;
    stop = group[0];
    count = 0;

    /* We remove all knots of each group and mark them in del_array
       until this would mean that we have removed too many. */

    while (stop <= antrem)
    {
      for (i=start; i<stop; i++) del_array[prio[i]] = 1;

      count++;

      if (count < antgr)
      {
	start = stop;
	stop = group[count];
      }
      else
      {

	/* start=stop signifies that we have reached the end. */

	stop = stop + 1;
	start = stop + 1;
      }
    }

    /* Now there are p more knots to remove from a group containing
       p knots. These p knots are removed "uniformly on index". */

    r = stop - start;
    p = antrem - start;

    if (p > 0)
    {
      h = (double) (r+1) / (double) p;
      for (i=0; i<p; i++)
      {
	indx = start - 1 + (int) floor( h*(i+0.5)+0.5 );
	del_array[prio[indx]] = 1;
      }
    }

    /* Gives the number of coefficients in the new spline to be computed. */

    hn = mprio - antrem ;

    /* The new knot vector is stored in ltau.  It might already be allocated
       from the last iteration, so free it first if required. */

    if (ltau != SISL_NULL) freearray(ltau);
    ltau = newarray(hn+2*k, double);
    if (ltau == SISL_NULL) goto err101;

    /* Set the first and last k knots. */

    for (i=0; i<k; i++)
    {
      ltau[i] = rankcurve->et[i];
    }

    /* Set the remaining knots as indicated by del_array. */

    for (indx=k, i=0; i<mprio; i++)
    {
       if (!del_array[i]) ltau[indx++] = rankcurve->et[i+k];
       else if (rankcurve->et[i+k] == tend_domain)
       {
	  /* Set knot at the seem that leads to a continutiy requirement
	     on the final spline.  */
	  
	  ltau[indx++] = tend_domain;
	  hn++;
	  kncont++;
       }
    }
    
    /* Compute an approximation on the new knot vector and store it
       in hcurve.
       s1943() will allocate space for hcurve (with icopy==1).
       Need local_err to store the error of the approximation until we know
       that we have an approximation with error smaller than the tolerance,
       then we can transfer the error to maxerr. */

    /* *jncont = MAX(kncont, *jncont); */
    *jncont = MIN(kncont, k-1);
    s1943(oldcurve, ltau, k, hn, startfix, endfix, *jncont,
	   &hcurve, local_err, l2_err, &lstat);
    if (lstat < 0) goto err;

    /* Check the error. big signifies that the error in one of the
       components is larger than the tolerance. bigco signifies that
       the error in one of the components is larger than the (usually
       much smaller) tolerance epsco. */
    big = 0;
    bigco = 0;

    for (i=0; i<dim; i++)
    {
      big = big || (local_err[i] > eps[i]);
      bigco = bigco || (local_err[i] > epsco[i]);
    }

    /* The error is too large if it is big or if it is bigco and the
       number of coefficients is smaller than nlim. The latter test is
       important when the sum of the number of derivatives to be kept
       fixed at the two ends is larger than k. */

    big = big || (bigco && hn < nlim);

    if (big)
    {

      /* If the error was too big, we just throw away hcurve and
	 indicate that an upper bound for the number of knots to
	 remove is antrem. */

      if (hcurve != SISL_NULL) freeCurve(hcurve);
      hcurve = SISL_NULL;
      maxi = antrem;
    }
    else
    {

      /* If the error is acceptable we know that we can remove at least
	 antrem knots and store hcurve in newcurve. We also save the
	 error in maxerr. */

      mini = antrem;
      if (*newcurve != SISL_NULL) freeCurve(*newcurve);
      *newcurve = hcurve;
      hcurve = SISL_NULL;
      memcopy(maxerr, local_err, dim, DOUBLE);
    }

    /* The number of knots to be removed next time is half way between
       mini and maxi. */

    antrem = mini + (maxi-mini)/2;
  }

  *stat = 0;

  goto out;



  /* Error in memory allocation. */

err101:
  *stat = -101;
  goto out;

  /* Error in lower level routine. */

err:
  *stat = lstat;
  s6err("s1950", *stat, pos);

  /* Clean up before exit. */

out:
  if (hcurve != SISL_NULL) freeCurve(hcurve);
  if (del_array != SISL_NULL) freearray(del_array);
  if (local_err != SISL_NULL) freearray(local_err);
  if (l2_err != SISL_NULL) freearray(l2_err);
  if (ltau != SISL_NULL) freearray(ltau);

  return;
}
