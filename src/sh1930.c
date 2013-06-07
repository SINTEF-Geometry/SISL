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
 * $Id: sh1930.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */

#define SH1930

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1930(double ea[],int nfirst[],int nlast[],double ed[],double ec[],
	     int ik,int in,int im,int idim,int ilend,int irend,int *jstat)
#else
void sh1930(ea,nfirst,nlast,ed,ec,ik,in,im,idim,ilend,irend,jstat)
   double ea[];
   int nfirst[];
   int nlast[];
   double ed[];
   double ec[];
   int ik;
   int in;
   int im;
   int idim;
   int ilend;
   int irend;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To modify a least squares problem to take into account 
*              equality constraints at the beginning and end of the vector
*              of unknowns.
* 
* INPUT      : ea     - Real array of dimension (im*ik) containing 
*                       the coefficient matrix of the problem. This matrix has
*                       dimension im*in but since at most
*                       ik entries are nonzero in each row, it can
*                       be stored in a im*ik array together
*                       with two integer arrays indicating the position
*                       of the first and last nonzero elements in each
*                       row. In addition it is known that the last nonzero
*                       element of row ki is stored in ea[ki*ik+ik-1].
*              nfirst - Integer array of dimension (im) containing 
*                       pointers to the first nonzero element of each row 
*                       of the original matrix.
*              nlast  - Integer array of dimension (im) containing 
*                       pointers to the last nonzero element of each row 
*                       of the original matrix.
*              ed     - Real array of dimension (in*idim) containing the 
*                       'right hand side' of the least squares problem.
*              ec     - Real array of dimension (in*idim) which is to contain
*                       the solution of the least squares problem. Upon entry
*                       of this routine ec[0],...ec[ilend] and 
*                       ec[in-irend],...,ec[in-1] (assuming idim=1) contain
*                       the values wich are forced by the constraints
*                       (the first ilend derivatives are to be kept fixed
*                       at the beginning of the curve, and the last irend at
*                       the end of the curve).
*	       ik     - The order of the spline space in the underlying least
*                       squares problem.
*              in     - The number of unknowns in the original least squares
*                       problem. This is reduced by (ilend+irend) due to the
*                       constraints.
*              im     - The number of constraints in the least squares problem
*                       excluding the (ilend+irend) side constraints.
*              idim   - The dimension of the geometry space.
*              ilend  - The number of coefficients at the beginning of ec
*                       which have a value.
*              irend  - The number of coefficients at the end of ec
*                       which have a value.
*              
*
* 
* OUTPUT     : ed     - Real array of dimension (im*idim) containing the new
*                       'right hand side' of the problem, which takes into
*                       consideration the side constraints.
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
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt, SI, 05.92, on the basis of a routine
*              written by Tom Lyche and Knut Moerken, 12.85.
*
*********************************************************************
*/
{ 
   int ki,kj,kr;
   int kjst;
   int knk1;
   double th;
   double *shelp = SISL_NULL;

   /* Check input.  */
   
   if (ilend + irend >= in) goto err203;

   /* Allocate scratch for help array and set to zero. */
			    
   if ((shelp = new0array(idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Make adjustments for the fact that ec[0],...,ec[ilend-1] are known
      This is done by going through the rows of the coefficient matrix until
      the first nonzero element is in column ilend or further to the right. */
   
   for (ki=0; ; ki++)
     {
	/* The maximum column is im-1.  */
	
	if (ki >= im) break;
	
	/* Check if the first nonzero element is far to the right. */
	
	if (nfirst[ki] > ilend) break;
	
	/* kjst points to the last column among the first ilend that has 
	   a nonzero element in row ki.     */
	
	kjst = MIN(ilend,nlast[ki]);
	
	/* Compute the necessary adjustment for row no. ki by multiplying
	   together the appropriate elements fo ec and ea and
	   accumulate.  */
	
	for (kj=nfirst[ki]; kj<=kjst; kj++)
	  {
	     th = ea[ki*ik+ik-nlast[ki]+kj-1];
	     for (kr=0; kr<idim; kr++)
	       shelp[kr] += th*ec[kj*idim+kr];
	  }
	
	/* Adjust ed[ki].  */
	
	for (kr=0; kr<idim; kr++)
	  {
	     ed[ki*idim+kr] -= shelp[kr];
	     shelp[kr] = (double)0.0;
	  }
     }
   
   /* Make adjustments for the fact that ec[in-irend],...,ec[in-1] are known.
      This is done by going through the rows of the coefficient matrix from
      the right until the last nonzero element is in column in-irend-1 or
      further to the left.  */
   
   for (knk1=in-irend, ki=im-1; ; ki--)
     {
	/* In any case we can stop when we get passed column 0.  */
	
	if (ki < 0) break;
	
	/* Check if the last nonzero element is sufficiently far to the
	   right. */
	
	if (nlast[ki] < knk1) break;
	
	/* kjst points to the first element in this row that is of
	   interest. */
	
	kjst = MAX(knk1,nfirst[ki]);
	
	/* Multiply together the appropriate elements fo ea and ec and
	   accumulate.  */
	
	for (kj=kjst; kj<=nlast[ki]; kj++)
	  {
	     th = ea[ki*ik+ik-nlast[ki]+kj-1];
	     for (kr=0; kr<idim; kr++)
	       shelp[kr] += th*ec[kj*idim+kr];
	  }
	
	/* Adjust ed[ki].  */
	
	for (kr=0; kr<idim; kr++)
	  {
	     ed[ki*idim+kr] -= shelp[kr];
	     shelp[kr] = (double)0.0;
	  }
     }
   
   /* Modification performed.  */
	      
   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;
   
   /* Error in input.  */
   
   err203: *jstat = -203;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (shelp != SISL_NULL) freearray(shelp);
	  
      return;
}
   
