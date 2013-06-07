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
 * $Id: sh1923.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */

#define SH1923

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1923(double *ea,int in,int ik,int *nstart,int *jstat)
#else
void sh1923(ea,in,ik,nstart,jstat)
   double *ea;
   int in;
   int ik;
   int *nstart;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To caluclate the Cholesky factorization (A=L*L(tr)) of
*              a (symmetric positive definite) in*in matrix with a
*              maximum of 2*ik-1 nonzero elements in each row. Because 
*              of the  symmetry it is sufficient to store the nonzero
*              elements to the left of and including the diagonal, i.e.
*              at most ik nonzero elements. Hence, the input matrix is
*              given by the two arrays ea and nstart of dimension (in*ik)
*              and in, respectively, with nstart indicating the position
*              of the first nonzero elements of each row of ea. If the
*              factorization is successfull, it is given as output in ea
*              a nstart using the same representation as the input.
* 
* 
* INPUT      : ea     - Real array of dimension (in*ik) containing the 
*                       nonzero elements of the matrix factorized. This
*                       is possible since it is assumed that the input 
*                       matrix is symmetric with at most 2*ik nonzero
*                       elements in each row.
*              in     - The dimension of the input matrix, i.e. the number
*                       of rows in ea.
*              ik     - The maximum number of different nonzero elements
*                       in each row of the input matrix and therefore the
*                       number of columns of ea.
*              nstart - Integer array of dimension (in) containing pointers
*                       to the first nonzero element of each row of ea.
*              
*
* 
* OUTPUT     : ea     - Real array of dimension (in*ik) containing the
*                       nonzero elements of the Cholesky factorization.
*                       This is possible because fo the special structure 
*                       of the input matrix and the fact that Cholesky
*                       factorization preserves this structure.
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
*
*********************************************************************
*/
{
   int ki,kj,kr;           /* Counters.              */
   int kjs,kjh,krhs,krh;   /* Pointers into matrix.  */
   int kik1 = ik-1;        /* Order minus one.       */
   double tsum;            /* Help variable.         */
   double thelp;           /* Help variable.         */

   
   /* Go through the rows of the matrix one by one.  */
   
   for (ki=0; ki<in; ki++)
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
   
   /* Cholesky factorization computed.  */
   
   *jstat = 0;
   goto out;
   
   /* Singular matrix.  */
   
   err106: *jstat = -106;
   goto out;
   
   out:
      return;
}
   
