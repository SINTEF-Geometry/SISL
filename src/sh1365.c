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
 * $Id: sh1365.c,v 1.2 2001-03-19 15:59:03 afr Exp $
 *
 */

#define SH1365

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1365(SISLCurve *pcurve,double etau[],int ik,int in,int ileftfix,
	     int irightfix,SISLCurve **rnewcurve,double **gmaxerr,
	     double **gl2err,int *jstat)
#else
void sh1365(pcurve,etau,ik,in,ileftfix,
	    irightfix,rnewcurve,gmaxerr,gl2err,jstat)
     SISLCurve *pcurve;
     double etau[];
     int ik;
     int in;
     int ileftfix;
     int irightfix;
     SISLCurve **rnewcurve;
     double **gmaxerr;
     double **gl2err;
     int *jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate an approximation to a given spline pcurve, f,
*              in a spline space S0 (given by the knot vector and number
*              of coefficients of the spline pcurve), from a subspace 
*              S1 (given by the knot vector etau which is a subseqence
*              of the knotvector of pcurve) of dimension in. Both the 
*              input and the output spline is of order ik. The approximation
*              rnewcurve, g, is to be a least squares approximation to pcurve
*              in the sense that a weighted sum of the squares of the 
*              B-spline coefficients of f-g (expressed as a spline on the
*              knotvector of pcurve, i.e. f) is minimized. In addition, the
*              0 - ileftfix-1 derivatives at the left end of the curve pcurve
*              and the 0 - irightfix-1 derivatives at the right end of the 
*              curve are to be kept fixed, i.e. there are ileftfix+irightfix
*              side conditions.
*
*
*
* INPUT      : pcurve    - Curve to be approximated on a subspace of the spline
*                          space in which the curve lies.
*              etau      - Knot vector corresponding to the subspace in which
*                          the approximating curve will lie.
*              ik        - Order of approximating curve. Should be the same as 
*                          the order of the input curve.
*              in        - Number of coefficients of the approximating spline.
*                          The dimension of the subspace.
*              ileftfix  - The number of derivatives that are to be kept fixed
*                          at the left end of the spline.
*              irightfix - The number of derivatives that are to be kept fixed
*                          at the right end of the spline.
*              
*                       
*
* OUTPUT     : rnewcurve  - Curve in the given subspace approximating the
*                           input spline curve.
*              gmaxerr    - Array of dimension equal to the dimension of the 
*                           geometry space, containing the absolute value of 
*                           the largest B-spline coefficient of the error f-g
*                           (see purpose above) in each component when f-g is
*                           expressed as a spline on the knot vector of f, 
*                           i.e. pcurve.
*              gl2err     - Array of dimension equal to the dimension of the
*                           geometry space, containing a weighted L2-norm of
*                           the B-spline coefficients of the error curve f-g.
*              jstat      - status messages  
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
*
* USE        : 
*
*-
* CALLS      : sh1922, sh1923, sh1924, sh1925, sh1926, sh1927, sh1928,
*              sh1930, newCurve.
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 07.90.
* REWRITTEN BY: Vibeke Skytt, SI, 05.92 on a basis of a routine by
*               Tom Lyche and Knut Moerken, 08.86 and 12.85.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Status variable. */
  int kdim = pcurve->idim;  /* Dimension of geometry space.          */
  int kkind = 1;            /* Indicates polynomial B-spline curve.  */
  int kcopy = 1;            /* Copy input arrays when creating a
			       B-spline curve.      */
  int kn = pcurve->in;      /* Number of vertices of input curve.    */
  double *scoef = SISL_NULL;     /* Coeffecient array of approximating curve. */
  double *sa = SISL_NULL;        /* Transformation matrix.  */
  double *sb = SISL_NULL;
  int *lfirst = SISL_NULL;
  int *llast = SISL_NULL;
  int *l2sta = SISL_NULL;
  double *sc = SISL_NULL;
  int kk = pcurve->ik;
  int knh = in - ileftfix - irightfix;

  /* Test order of input curve. */

  if (kk != ik) goto err109;

  /* Allocate scratch for error estimates.  */

  if ((*gmaxerr = new0array(kdim,DOUBLE)) == SISL_NULL) goto err101;
  if ((*gl2err = new0array(kdim,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Allocate scratch for coefficients of approximating spline.  */

  if ((scoef = newarray(in*kdim,DOUBLE)) == SISL_NULL) goto err101;
  
  /* If kn = in then pcurve->et = etau, and the problem is
     trivial, the solution is the input spline itself resulting in
     a zero error.  */
  
  if (kn == in)
  {
     memcopy(scoef,pcurve->ecoef,in*kdim,DOUBLE);
     
     /* Express the approximating curve as a curve.  */
     
     if ((*rnewcurve = newCurve(in,ik,etau,scoef,kkind,kdim,kcopy)) == SISL_NULL)
	goto err101;
     (*rnewcurve)->cuopen = pcurve->cuopen;
  }     
  else
  {
     /* Allocate scratch for local arrays. 
	Allocate scratch for a transformation matrix. This is a
	kn*in matrix, but since at most  elements
	are nonzero in each row, it can be stored as a kn*ik
	matrix together with two integer arrays of length kn
	indicating the position of the first and last nonzero elements
	of each row of the matrix.   
	Also allocate scratch for the coefficient matrix sb of the
	normal equations. This is a in*in symmetric positive definite
	matrix, but since at most ik elemnts are nonzero in each row, it
	can be stored as a in*ik matrix together with an integer array
	of length kn indicating the position of the first nonzero element
	of each row of sb.    */
     
     if ((sa = newarray(kn*ik,DOUBLE)) == SISL_NULL) goto err101;
     if ((sb = newarray(in*ik,DOUBLE)) == SISL_NULL) goto err101;
     if ((lfirst = newarray(kn,INT)) == SISL_NULL) goto err101;
     if ((llast = newarray(kn,INT)) == SISL_NULL) goto err101;
     if ((l2sta = newarray(kn,INT)) == SISL_NULL) goto err101;
     
     /* Compute the refinement matrix sa from etau to pcurve->et.  */
     
     sh1922(pcurve->et,kn,ik,etau,in,sa,lfirst,llast,&kstat);
     if (kstat < 0) goto error;
		    
     /* Set up the normal equations.  */
		    
     if (ileftfix == 0 && irightfix == 0)
     {
	/* If there is no side constraints, the number of unknowns is
	   in and the first unknown is no. 1. */
	
	sh1926(etau,ik,in,kdim,pcurve->et,pcurve->ecoef,kn,
	       sa,lfirst,llast,sb,l2sta,scoef,&kstat);
	if (kstat < 0) goto error;
      }
     else
     {
	/* Make a local copy of the coefficients of the input spline.
	   This is necessary since the right-hand-side of the problem
	   will be adjusted and the original data should not be altered. */
	
	if ((sc = newarray(kn*kdim,DOUBLE)) == SISL_NULL) goto err101;
	memcopy(sc,pcurve->ecoef,kn*kdim,DOUBLE);
     
	/* Enforce the side conditions.  */
	
	sh1927(etau,ik,in,kdim,pcurve,ileftfix,irightfix,scoef,&kstat);
	if (kstat < 0) goto error;
		       
        /* There are now in-ileftfix-irightfix unknowns left, and the first
            of them is no. ileftfix+1.  */
		       
        if (knh > 0)
	{
	   /* Adjust the left squares problem.  */
	   
	   sh1930(sa,lfirst,llast,sc,scoef,ik,in,kn,kdim,ileftfix,
		  irightfix,&kstat);
	   if (kstat < 0) goto error;
			  
           /* Set up the normal equations for the modified least
	      squares problem.   */
			  
	   sh1928(etau,ik,in,kdim,pcurve->et,sc,kn,ileftfix,
		  irightfix,sa,knh,lfirst,llast,sb,scoef,l2sta,&kstat);
	   if (kstat < 0) goto error;
        }
     }
     
     if (knh > 0)
     {
	/* Calculate the Cholesky factorization of the coefficient
	   matrix of the normal equations.  */
	
	sh1923(sb,knh,ik,l2sta,&kstat);
	if (kstat < 0) goto error;
		       
        /* Solve the normal equations.  */
		       
        sh1924(sb,scoef+ileftfix*kdim,knh,ik,kdim,l2sta,&kstat);
	if (kstat < 0) goto error;
      }

     /* Express the approximating curve as a curve.  */
     
     if ((*rnewcurve = newCurve(in,ik,etau,scoef,kkind,kdim,kcopy)) == SISL_NULL)
	goto err101;
     (*rnewcurve)->cuopen = pcurve->cuopen;
     
     /* Multiply by dtau(-1/2) and estimate error.  */
     
     sh1925(pcurve,*rnewcurve,kdim,sa,lfirst,llast,*gmaxerr,*gl2err,
	    ileftfix,irightfix,&kstat);
     if (kstat < 0) goto error;
  }
  
  /* Approximation performed.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Order of input and output curve in conflict.  */

  err109 :
    *jstat = -109;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :

    /* Free scratch occupied by local arrays.  */

    if (scoef != SISL_NULL) freearray(scoef);
    if (sc != SISL_NULL) freearray(sc);
    if (sa != SISL_NULL) freearray(sa);
    if (sb != SISL_NULL) freearray(sb);
    if (lfirst != SISL_NULL) freearray(lfirst);
    if (llast != SISL_NULL) freearray(llast);
    if (l2sta != SISL_NULL) freearray(l2sta);
  
    return;
}
