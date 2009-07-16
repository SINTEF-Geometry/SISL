/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1942

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1942(SISLCurve *pc1,SISLCurve *pc2,int idim,double ea[],int nstart[],
	     int nstop[],double emxerr[],double el2err[],int *jstat)
#else
void s1942(pc1,pc2,idim,ea,nstart,nstop,emxerr,el2err,jstat)
   SISLCurve *pc1;
   SISLCurve *pc2;
   int idim;
   double ea[];
   int nstart[];
   int nstop[];
   double emxerr[];
   double el2err[];
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To calculate the discrete max-error
*              and L2-error of the approximation.
* 
* 
* INPUT      : pc1    - The spline curve to be approximated.
*              pc2    - The curve approximating pc1 on a subset of the
*                       knotvector.
*              idim   - The dimension of the geometry space.
*              ea     - Real array of dimension (pc1->in*pc1->ik) containing 
*                       the B-spline refinement matrix from the knot vector
*                       pc2->et to the knot vector pc1->et. This matrix has
*                       dimension pc1->in*pc2->in but since at most
*                       pc1->ik entries are nonzero in each row, it can
*                       be stored in a pc1->in*pc1->ik array together
*                       with two integer arrays indicating the position
*                       of the first and last nonzero elements in each
*                       row.
*              nstart - Integer array of dimension (pc1->in) containing 
*                       pointers to the first nonzero element of each row 
*                       of the B-spline refinement matrix from pc2->et to
*                       pc1->et.
*              nstop  - Integer array of dimension (pc1->in) containing 
*                       pointers to the last nonzero element of each row 
*                       of the B-spline refinement matrix from pc2->et to
*                       pc1->et.
*
* 
* OUTPUT     : emxerr - Real array of dimension (idim) containing the 
*                       absolute value of the largest B-spline coefficient
*                       of f-g in each component when
*                       f-g is expressed as a spline on the pc1->et knot
*                       vector.
*              el2err - Real array of dimension (idim) containing a
*                       weighted L2-norm of the B-spline coefficients of
*                       f-g in each component when f-g is expressed as a
*                       spline on the knot vector pc1->et.
*              jstat     - status messages 
*                          > 0 : warning 
*                          = 0 : ok 
*                          < 0 : error 
*             
* 
* METHOD     : 
*
*
* REFERENCES : Any book on general numerical analysis or numerical
*              linear algebra.
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt, SI, 05.92, on the basis of a routine
*              written by Tom Lyche and Knut Moerken, 12.85.
* CHANGED AND RENAMED BY : Vibeke Skytt, SINTEF Oslo, 12.94. Introduced
*                                        periodicity.
*
*********************************************************************
*/
{
   int ki,kj,kr;
   int kjh;
   int kk = pc1->ik;
   int km = pc1->in;
   int kj1,kj2;
   double tkindv = (double)1.0/(double)kk;
   double thelp;
   double *st = pc1->et;
   double *sd = pc1->ecoef;
   double *sc = pc2->ecoef;
   double *stemp = SISL_NULL;
   
   /* Allocate scratch for local array.  */
   
   if ((stemp = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
  
   /* Initiate arrays to zero.  */
   
   memzero(stemp,idim,DOUBLE);
   memzero(emxerr,idim,DOUBLE);
   memzero(el2err,idim,DOUBLE);
		  
   /* Express the approximating spline as a spline on et by multiplying
      ec by ea and then calculate the error in the spline approximation. */
   
   for (ki=0; ki<km; ki++)
     {
	memzero(stemp,idim,DOUBLE);
	
	/* Express the approximation as a spline on et.  */
	
	kj1 = nstart[ki];
	kj2 = nstop[ki];
	for (kjh=kk+kj1-kj2-1, kj=kj1; kj<=kj2; kjh++, kj++)
	  {
	     thelp = ea[ki*kk+kjh];
	     for (kr=0; kr<idim; kr++)
	       stemp[kr] += thelp*sc[kj*idim+kr];
	  }
	
	/* Calculate the maxerror and the weighted L2-error of the
	   approximation. */
	
	thelp = (st[ki+kk] - st[ki])*tkindv;  
	for (kr=0; kr<idim; kr++)
	  {
	     stemp[kr] = fabs(stemp[kr] - sd[ki*idim+kr]);
	     el2err[kr] += thelp*stemp[kr]*stemp[kr];
	     if (stemp[kr] > emxerr[kr]) emxerr[kr] = stemp[kr];
	  }
     }
   for (kr=0; kr<idim; kr++)
     el2err[kr] = sqrt(el2err[kr]);
	
   
   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (stemp != SISL_NULL) freearray(stemp);
	  
      return;
}
   
