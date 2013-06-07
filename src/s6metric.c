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
      
