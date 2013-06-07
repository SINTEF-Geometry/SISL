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


#define S1943

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1943(SISLCurve *pcurve,double etau[],int ik,int in,int ileftfix,
	     int irightfix,int incont,SISLCurve **rnewcurve,
	     double gmaxerr[],double gl2err[],int *jstat)
#else
void s1943(pcurve,etau,ik,in,ileftfix,
	    irightfix,incont,rnewcurve,gmaxerr,gl2err,jstat)
     SISLCurve *pcurve;
     double etau[];
     int ik;
     int in;
     int ileftfix;
     int irightfix;
     int incont;
     SISLCurve **rnewcurve;
     double gmaxerr[];
     double gl2err[];
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
*              side conditions, and there are incont continuity conditions
*              at the seam of the spline curve.
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
*              incont    - Number of continuity conditions over the seam.
*                          incont < ik.
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
* REFERENCES : T. Lyche and K. Moerken. A Data Reduction Strategy for
*              Splines. February 1987.
*              
*
* USE        : 
*
*-
* CALLS      : 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 07.90.
* REWRITTEN BY: Vibeke Skytt, SI, 05.92 on a basis of a routine by
*               Tom Lyche and Knut Moerken, 08.86 and 12.85.
* REWISED AND RENAMED BY : Vibeke Skytt, SINTEF Oslo, 12/94. Introduced
*                          continuity conditions originating from periodicity.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Status variable. */
  int kdim = pcurve->idim;  /* Dimension of geometry space.          */
  int kj;                   /* Counter.                              */
  int kkind = 1;            /* Indicates polynomial B-spline curve.  */
  int kcopy = 1;            /* Copy input arrays when creating a
			       B-spline curve.      */
  int kn = pcurve->in;      /* Number of vertices of input curve.    */
  int knlr;                 /* Number of rows in corner element.     */
  int knred = 0;            /* Number of equations removed because of 
			       periodicity. */
  int knormr;               /* Number of rows of corner element of normal
			       equations.                                */
  double *scoef = SISL_NULL;     /* Coefficient array of approximating curve. */
  double *sa = SISL_NULL;        /* Transformation matrix.  */
  double *sa2 = SISL_NULL;       /* Transformation matrix. Copy of last part of sa.*/
  double *sw1 = SISL_NULL;       /* Corner element due to continuity requirements. */
  double *sb = SISL_NULL;        /* Coefficient matrix of the normal equations. */
  double *sw2 = SISL_NULL;       /* Corner elements of normal equations.        */
  double *sfac = SISL_NULL;      /* Factors used to implement continuity 
			       requirements.                               */
  int *lfirst = SISL_NULL;       /* Integer array of dimension kn containing
			       pointers to the first nonzero element of each
			       row of the B-spline refinement matrix, sa. from
			       etau to pcurve->et.  */
  int *lfirst2 = SISL_NULL;      /* Copy of last part of lfirst.               */
  int *llast = SISL_NULL;        /* Pointer to the last nonzero element of sa. */
  int *llast2 = SISL_NULL;       /* Copy of last part of llast to avoid destroing 
			       the original information of llast.         */
  int *l2sta = SISL_NULL;        /* Pointer to the first non-zero element of the
			       coefficient matrix of the normal equations. */
  double *sc = SISL_NULL;        /* Copy of coefficients of the input curve.    */
  int kk = pcurve->ik;      /* Order of curves.                            */
  int knh;                  /* Size of coefficient matrix in matrix solver. */

  /* Test order of input curve. */

  if (kk != ik) goto err109;

  /* Make sure that the conditions on continuity and side constraints
     are consistent. */
  
  if (incont < 0) incont = 0;
  if (incont >= ik) incont = ik-1;
  
  if (ileftfix >= incont && irightfix >= incont) incont = 0;
  else if (irightfix >= ileftfix && irightfix >= incont)
  {
     ileftfix = incont;
     irightfix -= incont;
  }
  else if (irightfix > ileftfix)
  {
     ileftfix = irightfix;
     irightfix = 0;
  }
  else if (incont >= irightfix) irightfix = 0;
  knh = in - ileftfix - irightfix - incont;
  
  /* Allocate scratch for coefficients of approximating spline.  */

  if ((scoef = newarray(in*kdim,DOUBLE)) == SISL_NULL) goto err101;
  
  /* If kn = in then pcurve->et = etau and if no continuity 
     requirements are given, the problem is trivial. 
     The solution is the input spline itself resulting in
     a zero error. If any continuity requirements are given,
     we run the approximation in order to get an exact fit
     at the seem.  */
  
  if (kn == in && incont == 0)
  {
     memcopy(scoef,pcurve->ecoef,in*kdim,DOUBLE);
     
     /* Express the approximating curve as a curve.  */
     
     if ((*rnewcurve = newCurve(in,ik,etau,scoef,kkind,kdim,kcopy)) == SISL_NULL)
	goto err101;
     (*rnewcurve)->cuopen = pcurve->cuopen;
     
     memzero(gmaxerr,kdim,DOUBLE);
     memzero(gl2err,kdim,DOUBLE);
  }     
  else
  {
     /* Allocate scratch for local arrays. 
	Allocate scratch for a transformation matrix. This is a
	kn*in matrix, but since at most ik elements
	are nonzero in each row, it can be stored as a kn*ik
	matrix together with two integer arrays of length kn
	indicating the position of the first and last nonzero elements
	of each row of the matrix.   
	Also allocate scratch for the coefficient matrix sb of the
	normal equations. This is a in*in symmetric positive definite
	matrix, but since at most ik elemnts are nonzero in each row, it
	can be stored as a in*ik matrix together with an integer array
	of length kn indicating the position of the first nonzero element
	of each row of sb.  
	In the periodic case (continuity requirements are given) it is
	necessary to allocate scratch for storing corner elements. */
     
     if ((sa = newarray((kn+in)*ik+incont*incont,DOUBLE)) == SISL_NULL) goto err101;
     sb = sa + kn*ik;
     if ((lfirst = newarray(3*kn,INT)) == SISL_NULL) goto err101;
     llast = lfirst + kn;
     l2sta = llast + kn;
     
     if (incont > 0)
     {
	sfac = sb + in*ik;
	memzero(sfac, incont*incont, DOUBLE);
     }
     
     /* Compute the refinement matrix sa from etau to pcurve->et.  */
     
     sh1922(pcurve->et,kn,ik,etau,in,sa,lfirst,llast,&kstat);
     if (kstat < 0) goto error;
		
     /* Compute size of matrix storing corner element. */
     
     for (knormr=in, kj=kn-1; kj>=0; kj--)
     {
	if (llast[kj] < in-incont) break;
	knormr = MIN(knormr, lfirst[kj]);
     }
     knlr = kn - kj - 1;
     knormr = MIN(knh, in - knormr);
     
     
     if (knlr > 0)
     {
	
	/* Allocate scratch for corner element matrix. */
	if ((sa2 = new0array(knlr*(ik+incont)+knormr*in,DOUBLE)) == SISL_NULL) 
	   goto err101;
	sw1 = sa2 + knlr*ik;
	sw2 = sw1 + knlr*incont;
	if ((lfirst2 = newarray(2*knlr,INT)) == SISL_NULL) goto err101;
	llast2 = lfirst2 + knlr;
     
     /* Save the content of sa to be used for estimating the error.  */
     
	memcopy(sa2, sa+(kn-knlr)*ik, knlr*ik, DOUBLE);
	
	/* Make copy of lfirst and llast. */
	
	memcopy(lfirst2, lfirst+kn-knlr, knlr, INT);
	memcopy(llast2, llast+kn-knlr, knlr, INT);
     }
     
     if (incont > 0)
     {
	
	/* Set up as constraints expressing that incont derivatives of
	   the curve is to be equal across the seam. */
	
	s1947(sa, lfirst, llast, ik, kn, etau, in, incont, sw1, 
	      knlr, &knred, sfac, &kstat);
	if (kstat < 0) goto error;
     }
     
     /* Set up the normal equations.  */
		    
     if (ileftfix == 0 && irightfix == 0)
     {
	/* If there is no side constraints, the number of unknowns is
	   in and the first unknown is no. 1. */
	
	s1944(etau,ik,in-incont,kdim,pcurve->et,pcurve->ecoef,kn,incont,
	      knlr,knormr,sa,sw1,lfirst,llast,sb,sw2,l2sta,scoef,&kstat);
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
	   
	   s1946(sa,sw1,lfirst,llast,sc,scoef,ik,in,kn,kdim,ileftfix,
		  irightfix,knlr,incont,&kstat);
	   if (kstat < 0) goto error;
			  
           /* Set up the normal equations for the modified least
	      squares problem.   */
			  
	   s1945(etau,ik,in,kdim,pcurve->et,sc,kn,ileftfix,irightfix,
		 incont,knlr,knormr,sa,sw1,knh,lfirst,llast,sb,sw2,scoef,
		 l2sta,&kstat);
	   if (kstat < 0) goto error;
        }
     }
     
     if (knh > 0)
     {
	/* Calculate the Cholesky factorization of the coefficient
	   matrix of the normal equations.  */
	
	s1948(sb,sw2,knh,ik,(knlr==knred)?0:knormr,l2sta,&kstat);
	if (kstat < 0) goto error;
		       
        /* Solve the normal equations.  */
		       
        s1949(sb,sw2,scoef+ileftfix*kdim,knh,ik,(knlr==knred)?0:knormr,
	      kdim,l2sta,&kstat);
	if (kstat < 0) goto error;
      }
     
     /* Multiply the coefficients by dtau(-1/2) and express the incont 
	last coefficients as a weighted sum of the incont first. */
     
     s1951(etau, scoef, in, ik, kdim, ileftfix, irightfix, incont, sfac);

     /* Express the approximating curve as a curve.  */
     
     if ((*rnewcurve = newCurve(in,ik,etau,scoef,kkind,kdim,kcopy)) == SISL_NULL)
	goto err101;
     (*rnewcurve)->cuopen = pcurve->cuopen;
     
     /* Restore information of sa, lfirst and llast. */
     
     if (knlr > 0)
     {
	memcopy(sa+(kn-knlr)*ik, sa2, knlr*ik, DOUBLE);
	memcopy(lfirst+kn-knlr, lfirst2, knlr, INT);
	memcopy(llast+kn-knlr, llast2, knlr, INT);
     }
     
     /* Estimate error.  */
     
     s1942(pcurve,*rnewcurve,kdim,sa,lfirst,llast,gmaxerr,gl2err,&kstat);
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
    if (sa2 != SISL_NULL) freearray(sa2);
    if (lfirst != SISL_NULL) freearray(lfirst);
    if (lfirst2 != SISL_NULL) freearray(lfirst2);
  
    return;
}
