/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1345.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */

#define S1345

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void s1345(SISLSurf *oldsurf,double eps[],int edgefix[4],double edgeps[],
	   double epsco,int opt,int itmax,SISLSurf **newsurf,
	   double maxerr[],int *stat)
#else
void s1345(oldsurf, eps, edgefix, edgeps, epsco, opt, itmax, newsurf,
	    maxerr, stat)
     SISLSurf 	*oldsurf;
     double	eps[];
     int	edgefix[4];
     double	edgeps[];
     double	epsco;
     int 	opt;
     int	itmax;
     SISLSurf	**newsurf;
     double	maxerr[];
     int	*stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* Purpose: To remove as many knots as possible from the surface oldsurf
*          without perturbing the surface more than the given tolerance.
*          The tolerances are given by eps and epsfix, while the
*	   approximation is given by newsurf.
*          Mark that the error in continuity over the start and end of
*          a closed or periodic surface is only guaranteed to be within
*          edgeps.
*
*
*
* Input:       oldsurf - pointer to the original spline surface. Note
*			 if the polynomial orders of the surface are
*			 k1 and k2, then the two knot vectors are
*			 assumed to have knots of multiplicity k1 and
*			 k2 at the ends.
*
*	       eps     - double array of length dim (the number of
*			 components of the surface, typically three)
*			 giving the desired accuracy of the
*			 final approximation compared to oldcurve.
*                        Note that in such comparisons the two
*			 surfaces are not reparametrized in any way.
*
*	       edgefix - integer array of dimension (4) giving the number
*			 of derivatives to be kept fixed along each edge
*			 of the surface. The numbering of the edges is the
*			 same as for edgeps below. All the derivatives of
*			 order < nend(i)-1 will be kept fixed along
*			 edge i. Hence nend(i)=0 indicates that nothing is
*			 to be kept fixed along edge i.
*
*                NB! TO BE KEPT FIXED HERE MEANS TO HAVE ERROR LESS THAN
*		     EDGEPS. IN GENERAL, IT IS IMPOSSIBLE TO REMOVE KNOTS
*                    AND KEEP AN EDGE COMPLETELY FIXED.
*
*	       edgeps  - double array of length 4*dim ([4,dim]) (dim is
*                        the number of components of each coefficient)
*			 containing the maximum deviation which is
*			 acceptable along the edges of the surface.
*                        edgeps[0]-edgeps[dim-1] gives the tolerance along
*			 the edge corresponding to x1 (the first parameter)
* 			 having it's minimum value.
*			 edgeps[dim]-edgeps[2*dim-1] gives the tolerance
*			 along the edge corresponding to x1 (the first
*			 parameter) having it's maximum value.
*              		 edgeps[2*dim]-edgeps[3*dim-1] gives the tolerance
*			 along the edge corresponding to x2 (the second
*			 parameter) having it's minimum value.
*              		 edgeps[3*dim]-edgeps[4*dim-1] gives the tolerance
*			 along the edge corresponding to x2 (the second
*			 parameter) having it's maximum value.
*           	 NB! EDGEPS WILL ONLY HAVE ANY SIGNIFICANCE IF THE
*		     CORRESPONDING ELEMENT OF EDGEFIX IS POSITIVE.
*
*	       epsco   - real number used to check equality between real
*                        numbers in some routines. Two numbers differing
*                        by a relative amount less than epsco will in some
*                        cases be considered equal. A suitable value is just
*                        above the unit roundoff of the machine. In IEEE
*                        double precision arithmetic 1.0E-15 is a reasonable
*                        choice.
*                 NB! THE COMPUTATIONS ARE NOT GUARANTEED TO HAVE
*                     RELATIVE ACCURACY LESS THAN EPSCO.
*
*
*
*	       itmax   - maximum number of iterations. The routine will
*                        follow an iterative procedure trying to remove
*                        more and more knots, one direction at a time.
*                        The process will almost always stop after less
*                        than 10 iterations and it will often stop after
*                        less than 5 iterations. A suitable value for itmax
*                        is therefore usually in the region 3-10.
*
*
*	       opt     - integer indicating the order in which the
*	       	         knot removal is to be performed.
*                            = 1 - remove knots in parameter direction 1 only.
*                            = 2 - remove knots in parameter direction 2 only.
*                            = 3 - remove knots first in parameter direction
*			       1 and then 2.
*                            = 4 - remove knots first in parameter direction
*			       2 and then 1.
*
*
*
*
*
* Output: 
*         newsurf      - the approximating surface on the reduced knot vectors.
*
*	  maxerr       - double array of length dim containing an upper
*		         bound for the pointwise error in each of the
*		         components of the spline approximation. The two
*                        surfaces oldsurf and newsurf are compared at the
*                        same parameter vaues, i.e., if oldsurf is f and
*                        newsurf is g then
*                              |f(u,v)-g(u,v)| <= eps
*                        in each of the components.
*
*         stat         - status message  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*
*
* Method:
*     The data reduction is performed in one parameter direction at a time,
*     using a subroutine for data reduction for spline curves (s1340).
*     Consider for example knot removal in the second parameter direction
*     (suppose oldcurve has m1xm2 coefficients, each with dim components
*     and suppose also that the knot vectors are t1 and t2).
*     This amounts to removing knots from each of the m1*dim curves, each
*     of length m2 and with knot vector t2, obtained by considering the
*     coefficients d as an array of dimension [m2][dim*m1] instead
*     of [m2][m1][dim] (note that in the code the arrays are always
*     treated as one dimensional).
*     A similar approach is possible when removing knots in the first parameter
*     direction except that then d has to be transposed in the first two
*     dimensions first.
*     The remaining part of the subroutine is more or less obvious except
*     perhaps the fact that if knots are to be removed in both directions,
*     the permissible error must be halved so that the total error committed
*     is kept within the original tolerance.
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
*-
* Calls: s1340, s1352, s6chpar, s6err
*
* Written by: Knut Moerken, University of Oslo, December 1992, based
*             on an earlier Fortran version by the same author.
*
*********************************************************************
*/                                     

{
  int k1 = oldsurf->ik1;          /* Unwrap some of the most commonly */
  int k2 = oldsurf->ik2;          /* used parameters.                 */
  int m1 = oldsurf->in1;
  int m2 = oldsurf->in2;
  int dim = oldsurf->idim;
  int ledgefix[4];                /* Local parameter for fixing the edge. */
  double *loct1 = oldsurf->et1;
  double *loct2 = oldsurf->et2;

  double *local_eps = NULL;       /* Declaration of local variables, */
  double *local_edge_eps = NULL;  /* for usage, see the code.        */
  double *clocal_eps = NULL;
  double *clocal_err = NULL;
  double *harray = NULL;
  double *tempcoef = NULL;
  SISLSurf *qs_kreg = NULL;

  SISLCurve *local_curve = NULL;
  SISLCurve *newlcurve = NULL;

  int lopt, i, it, n1, n2, antit, j, jh, lstat;
  double factor;
  int pos = 0;                    /* Parameter in s6err.             */

  /* Make sure that the input surface is non-periodic.  */
 
  memcopy(ledgefix, edgefix, 4, INT);
  if (oldsurf->cuopen_1 == SISL_CRV_PERIODIC ||
      oldsurf->cuopen_2 == SISL_CRV_PERIODIC)
  {
     make_sf_kreg(oldsurf, &qs_kreg, &lstat);
     if (lstat < 0) goto error;
     
     /* The input surface is closed and periodic. Make sure that the
	endcurves of the surface is matching as good as possible
	by fixing the position and the derivative in the edges of the surface.
	The change made to startfix and endfix is only locallly. */
     
     if (oldsurf->cuopen_1 == SISL_CRV_PERIODIC)
	ledgefix[0] = ledgefix[1] = MAX(MAX(edgefix[0],edgefix[1]), 2);
     if (oldsurf->cuopen_2 == SISL_CRV_PERIODIC)
	ledgefix[2] = ledgefix[3] = MAX(MAX(edgefix[2],edgefix[3]), 2);
  }
  else 
     qs_kreg = oldsurf;
  
  /* We start by allocating space for local arrays. */

  local_eps = newarray(dim, double);
  if (local_eps == NULL) goto err101;

  local_edge_eps = newarray(4*dim, double);
  if (local_edge_eps == NULL) goto err101;

  /* If we are going to remove knots in both directions, we use half the
     tolerance each time, otherwise the whole tolerance once. This is
     signified by factor. */

  factor = (opt < 3) ? 1.0 : 0.5;

  /* Store a local version of the tolerance and initialize the error
     to zero. */

  for (i=0; i<dim; i++)
    {
      local_eps[i] = factor*eps[i];
      maxerr[i] = 0.0;
    }

  /* We also need a local version of the edge tolerance. */

  /* Originally the edgetolerance was also multiplied by fac. It turns
     out that this is unnecessary: Each edge is completely fixed in one
     of the calls to s1340. */

  for (i=0; i<4*dim; i++) local_edge_eps[i] = edgeps[i];

  /* Now we are ready to prepare for the iterations. it counts the iterations
     while lopt says which direction to start with and antit tells us how
     many iterations there should be (one or two). n1 and n2 are used for
     keeping track of the number of coefficients. */

  it = 0; lopt = (opt-1) % 2 + 1;
  n1 = m1; n2 = m2;

  antit = (opt < 3) ? 1 : 2;
  while (it < antit)
    {
      if (lopt == 1)
	{

	  /* Here we remove knots along the first parameter direction. */

          /* Since we are going to treat the surface as a collection of
	     n2 curves in the first direction each with dim components,
	     we need a tolerance array of length dim*n2 and a similar
	     array for storing the error. */

	  clocal_eps = newarray(dim*n2, double);
	  if (clocal_eps == NULL) goto err101;

	  clocal_err = newarray(dim*n2, double);
	  if (clocal_err == NULL) goto err101;

	  /* The tolerance will in general not be the same for all the
	     n2*dim curves. If edges 1 and/or 3 are to be kept fixed then
	     first and/or last tolerance must be smaller than the tolerance
	     for middle curves. This is accomplished by a call to s1352. */

	  s1352(loct2, n2, k2, local_eps, (local_edge_eps+2*dim),
		 (local_edge_eps+3*dim), dim, ledgefix[2], ledgefix[3],
		 clocal_eps, &lstat);
	  if (lstat < 0) goto error;

	  /* First we have to pick out the correct coefficients. */

	  if (opt == 4)	    
	    harray = local_curve->ecoef;
	  else
	    harray = qs_kreg->ecoef;

	  /* Then we must transpose the coefficient with respect to
	     n1 and n2 due to the way we store the coefficients of
	     surfaces. */

	  tempcoef = newarray(n1*n2*dim, double);
	  if (tempcoef == NULL) goto err101;

	  s6chpar(harray, n1, n2, dim, tempcoef);

	  /* We then create a curve in which to store the high dimensional
	     curve. */

	  local_curve = newCurve(n1, k1, qs_kreg->et1, tempcoef,
				   1, dim*n2, 0);
	  if (local_curve == NULL) goto err101;

	  /* We can now perform knot removal on this curve. The result
	     is stored in newlcurve. */

	  s1340(local_curve, clocal_eps, ledgefix[0], ledgefix[1],
		 epsco, itmax, &newlcurve, clocal_err, &lstat);
	  if(lstat<0) goto error;

	  /* Remember to update n1 and throw away local_curve. */

	  n1 = newlcurve->in;
	  freeCurve(local_curve);

	  /* Now we must transpose back to return to surface coeffiecients. */

	  harray = newlcurve->ecoef;
	  tempcoef = increasearray(tempcoef,n1*n2*dim, double);
	  if (tempcoef == NULL) goto err101;

	  s6chpar(harray, n2, n1, dim, tempcoef);

	  /* We save the surface as a curve and create the final surface
	     later. */

	  local_curve = newlcurve;
	  local_curve->ecoef = tempcoef;

	  if (local_curve->et == qs_kreg->et1)
	    {

	      /* If no knots were removed and the original knot vector
		 was returned, we want a true copy of the original,
		 but not a reference to it. */

	      loct1 = newarray(n1+k1, double);
	      if (loct1 == NULL) goto err101;

	      harray = qs_kreg->et1;
	      for (i=0; i<n1+k1; i++) loct1[i] = harray[i];
	    }
	  else

	    /* Here we can just save the new knot vector in loct1. */

	    loct1 = local_curve->et;

	  /* Calculate the error in the approximation. */

	  /* The error in component i in the approximation is simply
	     the largest error in component i in any of the curves. */

	  for (i=0; i<dim; i++) clocal_eps[i] = 0.0;
	  for (jh=0, i=0; i<n2; i++)
	    for (j=0; j<dim; j++, jh++)
	      clocal_eps[j] = MAX(clocal_eps[j], clocal_err[jh]);

	  /* Accumulate the error in maxerr and use the remaining tolerance
	     in the next direction. */

	  for (i=0; i<dim; i++)
	    {
	      maxerr[i] += clocal_eps[i];
	      local_eps[i] = eps[i] - maxerr[i];
	    }
	  freearray(clocal_eps);
	  freearray(clocal_err);
	}
      else
	{


	  /* Now we are to remove knots in the second direction. This is
	     very similar to removing knots in the first direction but
	     we do not need to transpose. */

	  /* First some local arrays. */

	  clocal_eps = newarray(dim*n1, double);
	  if (clocal_eps == NULL) goto err101;

	  clocal_err = newarray(dim*n1, double);
	  if (clocal_err == NULL) goto err101;

	  /* Compute a tolerance along the second direction. */

	  s1352(loct1, n1, k1, local_eps, local_edge_eps,
		local_edge_eps+dim, dim, ledgefix[0], ledgefix[1],
		clocal_eps, &lstat);
	  if (lstat < 0) goto error;

	  /* Generate a high dimensional curve along the second direction
	     with the right coefficients. */

	  if (opt == 3)
	    local_curve = newCurve(n2, k2, qs_kreg->et2, local_curve->ecoef,
				   1, dim*n1, 0);
	  else
	    local_curve = newCurve(n2, k2, qs_kreg->et2, qs_kreg->ecoef, 1,
				   dim*n1, 0);
	  if (local_curve == NULL) goto err101;

	  /* Remove knots and store in newlcurve. */

	  s1340(local_curve, clocal_eps, ledgefix[2], ledgefix[3],
		epsco, itmax, &newlcurve, clocal_err, &lstat);
	  if(lstat<0) goto error;

	  /* Remember to update n2. */

	  n2 = newlcurve->in;

	  /* Throw away the old local_curve, but keep newlcurve. */

	  freeCurve(local_curve);
	  local_curve = newlcurve;

	  /* The new knot vector in the second direction will be stored in
	     locat2. */

	  loct2 = newarray(n2+k2, double);

	  /* Get the loct2 right. */

	  if (local_curve->et == qs_kreg->et2)
	    {
	      harray = qs_kreg->et2;
	      for (i=0; i<n2+k2; i++) loct2[i] = harray[i];
	    }
	  else
	    loct2 = local_curve->et;

	  /* Calculate the error in the approximation as above. */

	  for (i=0; i<dim; i++) clocal_eps[i] = 0.0;
	  for (jh=0, i=0; i<n1; i++)
	    for (j=0; j<dim; j++, jh++)
	      clocal_eps[j] = MAX(clocal_eps[j], clocal_err[jh]);

	  for (i=0; i<dim; i++)
	    {
	      maxerr[i] += clocal_eps[i];
	      local_eps[i] = eps[i] - maxerr[i];
	    }
	  freearray(clocal_eps);
	  freearray(clocal_err);
 	}

      /* Prepare for the next iteration. */

      ++it;
      lopt = 3 - lopt;
      
    }

  /* It remains to create a new surface object. If we only removed knots
     in one direction, one of loct1 and loct2 will be NULL. */

  if (loct1 == NULL)
    {
      loct1 = newarray(n1+k1, double);

      harray = qs_kreg->et1;
      for (i=0; i<n1+k1; i++) loct1[i] = harray[i];
    }

  if (loct2 ==NULL)
    {
      loct2 = newarray(n2+k2, double);

      harray = qs_kreg->et2;
      for (i=0; i<n2+k2; i++) loct2[i] = harray[i];
    }

  /* Generate the new surface object. */

  *newsurf = newSurf(n1, n2, k1, k2, loct1, loct2, local_curve->ecoef,
		    1, dim, 2);

  /* Set periodicity flag. */
  
  if (oldsurf->cuopen_1 == SISL_CRV_CLOSED || 
      oldsurf->cuopen_1 == SISL_CRV_PERIODIC)
     (*newsurf)->cuopen_1 = SISL_CRV_CLOSED;
  if (oldsurf->cuopen_2 == SISL_CRV_CLOSED || 
      oldsurf->cuopen_2 == SISL_CRV_PERIODIC)
     (*newsurf)->cuopen_2 = SISL_CRV_CLOSED;
  
  *stat = 0;
  goto out;

  /* Error in memory allocation. */

  err101:
    *stat = -101;
    s6err("s1345", *stat, pos);
    goto out;

  /* Error in lower level routine. */

  error:
    *stat = lstat;  
    s6err("s1345", *stat, pos);

  /* Clean up. */

  out:
    if (local_eps != NULL) freearray(local_eps);
    if (local_edge_eps != NULL) freearray(local_edge_eps);
    if (clocal_eps != NULL) freearray(clocal_eps);
    if (qs_kreg != NULL && qs_kreg != oldsurf) freeSurf(qs_kreg);

  return;

}




