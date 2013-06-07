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
 * $Id: sh1924.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */

#define SH1924

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1924(double *ea,double *eb,int in,int ik,int idim,
	     int *nstart,int *jstat)
#else
void sh1924(ea,eb,in,ik,idim,nstart,jstat)
   double *ea;
   double *eb;
   int in;
   int ik;
   int idim;
   int *nstart;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To solve idim linear systems of equations A*X=eb, given
*              the right-hand-side eb and the Cholesky factorization of
*              the matrix A (i.e. A is assumed to be symmetric and
*              positive definite). It is also assumed that A has at most
*              2*ik-1 nonzero elements in each row, and since this
*              structure is not destroyed by the Cholesky factorization,
*              it is sufficient to store a maximum of ik elements per
*              row, see the subroutine sh1923. The right-hand-side eb has
*              idim components each of length in. The output of the 
*              routine is the solution fo the linear system which 
*              overwrites eb.
* 
* 
* INPUT      : ea     - Real array of dimension (in*ik) containing the 
*                       nonzero elements of the Cholesky factorization.
*              eb     - Real array of dimension (in*idim) containing
*                       the right-hand-side(s) of the linear systems.
*              in     - The dimension of the linear systems, i.e. the number
*                       of rows in ea.
*              ik     - The maximum number of different nonzero elements
*                       in each row of the Cholesky factorization and 
*                       therefore the number of columns of ea.
*              idim   - The number of different right-hand-sides.
*              nstart - Integer array of dimension (in) containing pointers
*                       to the first nonzero element of each row of ea.
*
* 
* OUTPUT     : eb     - The solution of the idim linear systems.
*              jstat  - status messages 
*                          > 0 : warning 
*                          = 0 : ok 
*                          < 0 : error 
*             
* 
* METHOD     : The linear system is solved in the usual way, first a 
*              forward substition and the a back substitution.
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
   int ki,kj,kr;         /* Counters.    */
   int ki1,kih,kim;      /* Counters.    */
   int kjs,kjh; /* Pointers into matrix.  */
   int kik1 = ik-1;      /* Order minus one.       */
   double thelp;         /* Help variable.         */
   double *ssum=SISL_NULL;    /* Help array.            */

   /* Allocate scratch for help array.  */
   
   if ((ssum = new0array(idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Forward substitution.  */
   
   for (ki=0; ki<in; ki++)
     {
	ki1 = ki - 1;
	for (kr=0; kr<idim; kr++) ssum[kr] = (double)0.0;
	
	/* kjs points to the firs nonzero element of row ki of ea, and
	   kjh points to the corresponding element in the underlying
	   matrix.       */
	
	for (kjs=nstart[ki], kjh=ki+kjs-ik+1, kj=kjs;
	 kj<kik1; kjh++, kj++)
	  {
	     thelp = ea[ki*ik+kj];
	     for (kr=0; kr<idim; kr++)
	       ssum[kr] += thelp*eb[kjh*idim+kr];
	  }
	
	/* Check if the linear system is singular.  */
	
	if (DEQUAL(ea[ki*ik+kik1],DZERO)) goto err106;
	
	thelp = (double)1.0/ea[ki*ik+kik1];
	for (kr=0; kr<idim; kr++)
	  eb[ki*idim+kr] = (eb[ki*idim+kr] - ssum[kr])*thelp;
     }
   
   /* Back substitution.  */
   
   for (ki=in-1, ki1=in-1, kih=0;
    kih<in; ki--, kih++)
     {
	/* Compute the index of the last nonzero element in row ki of the
	   transpose of the Cholesky factorization. The integer ik-ki1+ki
	   gives the position of element no (ki1*ik+ki) in row ki1 of ea.
	   ki1 is reduced until the first nonzero element or row ki1 is
	   ik-ki1+ki.   */
	
	for (;;ki1--)
	  if (nstart[ki1] < ik-ki1+ki) break;
	
	for (kr=0; kr<idim; kr++) ssum[kr] = (double)0.0;
	
	/* Calculate eb[.*ik+kik1].  */
	
	for (kjs=ki+1, kim=ik-kjs+ki-1, kj=kjs;
	 kj<=ki1; kim--, kj++)
	  {
	     thelp = ea[kj*ik+kim];
	     for (kr=0; kr<idim; kr++)
	       ssum[kr] += thelp*eb[kj*idim+kr];
	  }
	
	thelp = (double)1.0/ea[ki*ik+ik-1];
	for (kr=0; kr<idim; kr++)
	  eb[ki*idim+kr] = (eb[ki*idim+kr] - ssum[kr])*thelp;
     }
   
   
   /* The linear system is solved.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;
   
   /* Singular matrix.  */
   
   err106: *jstat = -106;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (ssum != SISL_NULL) freearray(ssum);
	  
      return;
}
   
