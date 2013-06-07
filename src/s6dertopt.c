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
 * $Id: s6dertopt.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DERTOPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s6dertopt(double eder[],int ntype[],int inpt,int idim,
	       double epoint[],int *jstat)
#else
void s6dertopt(eder,ntype,inpt,idim,epoint,jstat)
   double eder[];
   int    ntype[];
   int    inpt;
   int    idim;
   double epoint[];
   int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Translate derivateve conditions in an interpolation
*              problem to position conditions.
*
*
*
* INPUT      : eder   - Array of interpolation conditions. The last
*                       condition must be a point. Dimension is inpt*idim.
*              ntype  - Array containing kind of condition. Dimension 
*                       is inpt.
*                       =  0 : A point is given.
*                       =  d : The d'th derivatative condition to the
*                              previous point is given.
*                       = -d : The d'th derivatative condition to the
*                              next point is given.
*              inpt   - Number of interpolation conditions.
*              idim   - Dimension of geometry space.
*              
*
* OUTPUT     : epoint - Positional interpolation conditions.
*                       The dimension is inpt*idim.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : s6lufacp  -  LU-factorization of matrix.
*              s6lusolp  -  Solve equation system.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-04.
*
*********************************************************************
*/
{
   int kstat = 0;        /* Status variable.                        */
   int ki,kj,kk,kh,kl;   /* Counters.                               */
   int kord;             /* Order of Bezier segment.                */
   int kder;             /* Order of derivative condition.          */
   int ksgn;             /* Indicates endpoint of Bezier segment.   */
   int *lpiv = SISL_NULL;     /* Pivot array.                            */
   double *sc = SISL_NULL;    /* Matrix of equation system.              */
   double *sd = SISL_NULL;    /* Right side of equation system.          */
   double *s1;           /* Pointer into matrix of equation system. */
   double *spt;          /* Pointer to interpolation condition.     */
   double *sdum1 = SISL_NULL; /* Help array.                             */
   double *sdum2 = SISL_NULL; /* Help array.                             */
   
   /* Test if the last condition is a point. */
   
   if (ntype[inpt-1] != 0) goto err151;
   
   /* Copy interpolation conditions to output array. */
   
   memcopy(epoint,eder,inpt*idim,DOUBLE);
   
   /* Allocate scratch for equation system. Make sure that the arrays
      are large enough.  */
   
   if ((sc = newarray(inpt*inpt,DOUBLE)) == SISL_NULL) goto err101;
   if ((sd = newarray(inpt,DOUBLE)) == SISL_NULL) goto err101;
   if ((lpiv = newarray(inpt,INT)) == SISL_NULL) goto err101;
   
   /* Allocate scratch for local help arrays. */
   
   if ((sdum1 = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   if ((sdum2 = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Traverse interpolation conditions. */
   
   for (ki=0; ki<inpt; ki=kj)
   {
      for (kj=ki+1; kj<inpt && ntype[kj]!=0; kj++);
      
      if (kj-ki > 1)
      {
	 /* A derivative condition is found. Express the segment between
	    two positional conditions as a Bezier segment in [0,1]. The 
	    order is given by the number of derivative conditions. First
	    put up matrix of equation system. */
	 
	kord = kj - ki + 1;  /* Order of Bezier segment.  */
	
	/* Traverse conditions.  */
	
	for (kk=0, s1=sc; kk<kord; kk++, s1+=kord)
	{
	   /* Set line of matrix to zero.  */
	   
	   for (kh=0; kh<kord; kh++) s1[kh] = DZERO;
	   
	   /* Fetch order of differentiation and endpoint of Bezier
	      segment.  */
	   
	   kder = abs(ntype[ki+kk]);
	   ksgn = (ntype[ki+kk] > 0) ? 1 : -1;
	   if (kk == 0) ksgn = 1;
	   
	   if (ksgn == 1)
	   {
	      /* Condition in startpoint of segment. Initiate to
		 positional condition.  */
	      
	      s1[0] = (double)1.0;
	      
	      /* Compute coefficients to derivative of segment recursively. */
	      
	      for (kh=0; kh<kder; kh++)
	      {
		 for (kl=kder; kl>0; kl--)
		    s1[kl] = (double)(kord-1)*(s1[kl-1] - s1[kl]);
		 s1[0] *= -(double)(kord-1);
	      }
	   }
	   else
	   {
	      /* Condition in endpoint of segment. */
	      
	      s1[kord-1] = (double)1.0;
	      for (kh=0; kh<kder; kh++)
	      {
		 for (kl=kord-kder-1; kl<kord-1; kl++)
		    s1[kl] = (double)(kord-1)*(s1[kl]-s1[kl+1]);
		 s1[kord-1] *= (double)(kord-1);
	      }
	   }
	}
	
	/* Perform LU-factorization of matrix.  */
	
	s6lufacp(sc,lpiv,kord,&kstat);
	if (kstat < 0) goto error;
		       
        for (kh=0; kh<idim; kh++)
	{
	   /* Set up right side of equation system.  */
	   
	   for (kk=0; kk<kord; kk++) sd[kk] = epoint[(ki+kk)*idim+kh];
	   
	   /* Solve equation system. */
	   
	   s6lusolp(sc,sd,lpiv,kord,&kstat);
	   if (kstat < 0) goto error;
			  
  	   /* Move result into condition array.  */
	
	   for (kk=1; kk<kord-1; kk++) epoint[(ki+kk)*idim+kh] = sd[kk];
	}
	
	/* Now the output array contains the coefficients of a Bezier
	   segment at the current position. Estimate positional
	   conditions that would give about the same Bezier segment. */
	
	memcopy(sdum1,epoint+ki*idim,idim,DOUBLE);
	for (kk=ki+1, spt=epoint+kk*idim; kk<kj; kk++, spt+=idim)
	{
	   /* Traverse inner coefficients of segment. */
	   
	   memcopy(sdum2,spt,idim,DOUBLE);
	   for (kh=0; kh<idim; kh++)
	      spt[kh] = (double)0.25*(sdum1[kh]+spt[idim+kh])
		 + (double)0.5*spt[kh];
	   memcopy(sdum1,sdum2,idim,DOUBLE);
	}
      }
   }
   
   /* Derivative conditions translated.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Illegal interpolation point.  */
   
   err151 : *jstat = -151;
   goto out;
   
   /* Error in lower level routine.  */
   
   error : *jstat = kstat;
   goto out;
   
   out :
      /* Free scratch occupied by local arrays. */
      
      if (sc != SISL_NULL) freearray(sc);
      if (sd != SISL_NULL) freearray(sd);
      if (lpiv != SISL_NULL) freearray(lpiv);
      if (sdum1 != SISL_NULL) freearray(sdum1);			
      if (sdum2 != SISL_NULL) freearray(sdum2);
			 
      return;
}
