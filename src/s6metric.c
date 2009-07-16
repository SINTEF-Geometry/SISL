/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s6metric.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6METRIC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s6metric(double epoint[],int in,int idim,double emat[],int *jstat)
#else
void s6metric(epoint,in,idim,emat,jstat)
   double epoint[];
   int    in;
   int    idim;
   double emat[];
   int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Set up matrix used for distance computation in an
*              affine metric.
*
*
*
* INPUT      : epoint - Point set used to set up the metric. The metric
*                       is to be used to compute distances between
*                       points in the set. Dimension is in*idim, and the
*                       points are stored in sequence.
*              in     - Number of points in the point set.
*              idim   - Dimension of geometry space.
*              
*
* OUTPUT     : emat   - Matrix used to compute distances in the metric.
*                       The dimension is idim*idim.
*              jstat  - status messages  
*                                         = 1      : singular equation system
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : s6invert - Find the inverse of a given square matrix.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-03.
*
*********************************************************************
*/
{
   int kstat = 0;        /* Status variable.  */
   int ki,kk;            /* Counter.          */
   double tmedium;       /* Medium of a number of coordinates. */
   double tdot;          /* Scalar product.                    */
   double *smat1 = SISL_NULL; /* Input points minus weight points.  */
   double *smat2 = SISL_NULL; /* Inverse matrix to output matrix.   */
   double *s1,*s2,*s3;   /* Pointers into arrays.              */
   
   /* Allocate scratch for internal matrices.  */
   
   if ((smat1 = newarray(in*idim,DOUBLE)) == SISL_NULL) goto err101;
   if ((smat2 = newarray(idim*idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Set up a matrix consisting of the input points minus
      the weight point of the pointset. First copy the array of points
      into the matrix.  */
   
   memcopy(smat1,epoint,in*idim,DOUBLE);
   
   for (kk=0; kk<idim; kk++)
   {
      /* Compute weight point of coordinate.  */
      
      for (tmedium=DZERO, s1=smat1+kk, s2=s1+in*idim; s1<s2; s1+=idim)
	 tmedium += (*s1);
      tmedium /= in;
      
      /* Set up column in difference matrix.  */
      
      for (s1=smat1+kk; s1<s2; s1+=idim) *s1 -= tmedium;
   }
   
   /* Compute the inverse of the output matrix as the matrix product
      between the difference matrix and the transpose of that matrix. */
   
   for (ki=0; ki<idim; ki++)
      for (kk=0; kk<idim; kk++)
      {
	 for (tdot=DZERO, s1=smat1+ki, s2=smat1+kk, s3=s1+in*idim;
	      s1<s3; s1+=idim, s2+=idim)
	    tdot += (*s1)*(*s2);
	 smat2[ki*idim+kk] = tdot;
      }
   
   /* Compute the matrix corresponding to the metric.  */
   
   s6invert(smat2,idim,emat,&kstat);
   if (kstat < 0) goto error;
		  
   /* The matrix of the metric is computed.  */
		  
   *jstat = kstat;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in lower level routine.  */
   
   error : *jstat = kstat;
   goto out;
   
   out :
      /* Free scratch occupied by local arrays.  */
      
      if (smat1 != SISL_NULL) freearray(smat1);
      if (smat2 != SISL_NULL) freearray(smat2);
			 
      return;
}
      
