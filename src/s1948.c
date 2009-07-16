/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"


#define S1948

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1948(double *ea,double *ew,int in,int ik,int inlr, 
	 int *nstart,int *jstat)
#else
void s1948(ea,ew,in,ik,inlr,nstart,jstat)
   double *ea;
   double *ew;
   int in;
   int ik;
   int inlr;
   int *nstart;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To calculate the Cholesky factorization (A=L*L(tr)) of
*              a (symmetric positive definite) in*in matrix with a band of
*              maximum of 2*ik-1 nonzero elements in each row. In addition
*              it might be a corner element of size inlr*inlr. Because 
*              of the  symmetry it is sufficient to store the nonzero
*              elements to the left of and including the diagonal, i.e.
*              at most ik nonzero elements, and the lower inlr rows. 
*              Hence, the input matrix is given by the three arrays
*              ea, ew and nstart of dimension (in*ik), (inlr*in)
*              and in, respectively, with nstart indicating the position
*              of the first nonzero elements of each row of ea. If the
*              factorization is successfull, it is given as output in ea
*              ew, and nstart using the same representation as the input.
* 
* 
* INPUT      : ea     - Real array of dimension (in*ik) containing the 
*                       nonzero band elements of the matrix factorized. This
*                       is possible since it is assumed that the input 
*                       matrix is symmetric with at most 2*ik nonzero
*                       elements in each row.
*              ew     - Corner element originating from periodic geometry.
*                       The size of ew is inlr*in. It is expected 
*                       that the last inlr rows of ea is copied into ew.
*              in     - The dimension of the input matrix, i.e. the number
*                       of rows in ea.
*              ik     - The maximum number of different nonzero elements
*                       in each row of the input matrix and therefore the
*                       number of columns of ea.
*              inlr   - The number of rows of the corner element
*              nstart - Integer array of dimension (in) containing pointers
*                       to the first nonzero element of each row of ea.
*              
*
* 
* OUTPUT     : ea     - Real array of dimension (in*ik) containing the
*                       nonzero elements of the Cholesky factorization.
*                       This is possible because of the special structure 
*                       of the input matrix and the fact that Cholesky
*                       factorization preserves this structure.
*              ew    -  The corner element and the last inlr rows of ea
*                       after the factorization. Due to fill in all inlr*in
*                       elements are expected to be non-zero.
*              jstat  - status messages 
*                          > 0 : warning 
*                          = 0 : ok 
*                          < 0 : error 
*             
* 
* METHOD     : The Cholesky factorization is computed in the form
*              A = L*L(tr). (L(tr) denotes the transpose of L) one
*              row at a time.
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
* REWRITTEN AND RENAMED BY : Vibeke Skytt, SINTEF Oslo, 12.94. Introduced
*                            a corner element originating from periodicity.
*
*********************************************************************
*/
{
   int ki,kj,kr;           /* Counters.              */
   int kjs,kjh,krhs,krh;   /* Pointers into ea matrix.  */
   int ki2, ki3;           /* Pointers into ew.      */
   int kik1 = ik-1;        /* Order minus one.       */
   double tsum;            /* Help variable.         */
   double thelp;           /* Help variable.         */

   
   /* Go through the rows of the matrix one by one excluding the rows
      storing the corner element.  */
   
   for (ki=0; ki<in-inlr; ki++)
     {
	/* Calculate the off diagonal elemnts of row no ki of the
	   Cholesky factorization. The first nonzero element of this
	   row is given by nstart[ki]. The integer kjh gives the
	   position of this element in the underlying in*in matrix.  */
	
	for (kjs=nstart[ki], kjh=ki+kjs-ik+1, krhs=kik1, kj=kjs;
	 kj<kik1; kjh++, krhs--, kj++)
	  {
	     tsum = (double)0.0;
	     krh = krhs;
	     for (kr=kjs; kr<kj; krh++,kr++)
	       tsum += ea[ki*ik+kr]*ea[kjh*ik+krh];
	     ea[ki*ik+kj] -= tsum;
	     ea[ki*ik+kj] /= ea[kjh*ik+kik1];
	  }
	
	/* Finish off row ki by calculating the diagonal element of the
	   Cholesky factorization.   */
	
	tsum = (double)0.0;
	for (kr=kjs; kr<kik1; kr++)
	  {
	     thelp = ea[ki*ik+kr];
	     tsum += thelp*thelp;
	  }
	
	/* Check if nonpositive ea[ki*ik+ik-1], i.e. singular non positive
	   definite matrix.  */
	
	tsum = ea[ki*ik+kik1] - tsum;
	if (tsum <= DZERO) goto err106;
        ea[ki*ik+kik1] = sqrt(tsum);
     }
   
   for (ki2=0; ki<in; ki++, ki2++)
     {
	/* Calculate the off diagonal elemnts of row no ki of the
	   Cholesky factorization. The current rows are stored in ew.
	   The first nonzero element of the band part of this
	   row (stored in ea) is given by nstart[ki]. The integer kjh gives 
	   the position of this element in the underlying in*in matrix. 
	   This loop is entered only if continuity requirements are
	   given across the seem.  */
	
	for (kjs=nstart[ki], kjh=ki+kjs-ik+1, krhs=kik1, kj=0;
	 kj<MIN(ki,in-inlr); kjh++, krhs--, kj++)
	  {
	     tsum = (double)0.0;
	     krh = nstart[kj];
	     for (kr=kj-kik1+nstart[kj]; kr<kj; krh++,kr++)
	       tsum += ew[ki2*in+kr]*ea[kj*ik+krh];
	     ew[ki2*in+kj] -= tsum;
	     ew[ki2*in+kj] /= ea[kj*ik+kik1];
	  }
	
	/* Compute the factorization of the last inlr elements. */
	
	for (ki3=0; kj<ki; kj++, ki3++)
	{
	     tsum = (double)0.0;
	     for (kr=0; kr<kj; kr++)
		tsum += ew[ki2*in+kr]*ew[ki3*in+kr];
	     ew[ki2*in+kj] -= tsum;
	     ew[ki2*in+kj] /= ew[ki3*in+kj];
	}
	
	/* Finish off row ki by calculating the diagonal element of the
	   Cholesky factorization.   */
	
	tsum = (double)0.0;
	for (kr=0; kr<ki; kr++)
	  {
	     thelp = ew[ki2*in+kr];
	     tsum += thelp*thelp;
	  }
	
	/* Check if nonpositive ea[ki*ik+ik-1], i.e. singular non positive
	   definite matrix.  */
	
	tsum = ew[ki2*in+ki] - tsum;
	if (tsum <= DZERO) goto err106;
        ew[ki2*in+ki] = sqrt(tsum);
     }
   
   /* Cholesky factorization computed.  */
   
   *jstat = 0;
   goto out;
   
   /* Singular matrix.  */
   
   err106: *jstat = -106;
   goto out;
   
   out:
      return;
}
   
