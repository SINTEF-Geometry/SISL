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
 * $Id: s6invert.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6INVERT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s6invert(double emat[],int im,double einvertmat[],int *jstat)
#else
void s6invert(emat,im,einvertmat,jstat)
   double emat[];
   int    im;
   double einvertmat[];
   int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the inverse of a given square matrix.
*
*
*
* INPUT      : emat   - The matrix of which to find the inverse.
*                       NB! The elements of emat will be changed!!!
*                       The dimension of the array is im*im.
*              im     - Dimension of the square matrix.
*              
*
* OUTPUT     : emat   - The contents of the input matrix is changed 
*                       during the computations. If you want to keep
*                       the original matrix, copy it before the inverse
*                       computation is performed.
*              einvertmat - The inverse of the input matrix. The
*                            dimension is im*im.
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
* CALLS      : s6lufacp - LU-factorization of matrix.
*              s6lusolp - Solve equation system given LU-factorization.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-03.
*
*********************************************************************
*/
{
   int kstat = 0;     /* Status variable.  */
   int ki;            /* Counter.          */
   int *lpiv = SISL_NULL;  /* Array indicating order of rows in the matrix
			 after pivoting.                               */
   double *scol = SISL_NULL;  /* Column of identity matrix.                 */
   double *s1,*s2,*s3;   /* Pointers into double array.                */
   
   /* Allocate scratch for local arrays.  */
   
   if ((lpiv = newarray(im,INT)) == SISL_NULL) goto err101;
   if ((scol = newarray(im,DOUBLE)) == SISL_NULL) goto err101;

   /* Perform LU-factorization of input matrix.  */
   
   s6lufacp(emat,lpiv,im,&kstat);
   if (kstat < 0) goto error;
   if (kstat == 1) goto sing;		  
		  
   /* Solve the equation systems using the LU-factorization of the
      input matrix as the left side and the columns of the identity
      matrix as the right side, to find the inverse of the input matrix. */
		  
   for (ki=0; ki<im; ki++)
   {
      /* Set up column of identity matrix.  */
      
      for (s1=scol, s2=s1+im; s1<s2; s1++) *s1 = DZERO;
      scol[ki] = (double)1.0;
      
      /* Solve equation system.  */
      
      s6lusolp(emat,scol,lpiv,im,&kstat);
      if (kstat < 0) goto error;
      if (kstat == 1) goto sing;
		     
      /* Copy contents of column to output matrix.  */
		     
      for (s1=scol, s2=s1+im, s3=einvertmat+ki;
           s1<s2; s1++, s3+=im)
	 *s3 = *s1;
   }
   
   /* Inverse of matrix computed.  */
   
   *jstat = 0;
   goto out;
   
   /* Singular equation system.  */
   
   sing : *jstat = 1;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in lower level routine.  */
   
   error : *jstat = kstat;
   goto out;
   
   out :
      /* Free space occupied by local arrays.  */
      
      if (lpiv != SISL_NULL) freearray(lpiv);
      if (scol != SISL_NULL) freearray(scol);
			
      return;
}
