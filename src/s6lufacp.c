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
 * $Id: s6lufacp.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6LUFACP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6lufacp(double ea[],int nl[],int im,int *jstat)
#else
void s6lufacp(ea,nl,im,jstat)
     double ea[];
     int    nl[];
     int    im;
     int    *jstat;
#endif
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : LU-factorisation of a full matrix. The matrix is stored 
*             in a one-dimensional array. 
*
*
*   INPUT   : im   - The number of equations.
*
*
*
*   INPUT/OUTPUT : eb - the array in which the matrix is stored.
*                       The elements are stored row-wise.
*
*   
*   OUTPUT  : nl    - Array that indicates order of rows after pivoting.
*             jstat - Status variable.
*                       < 0 : error
*                       = 0 : ok
*                       > 0 : warning
*
*
*   METHOD  : Gauss-elimination with partial pivoting.
*
*
*   REFERENCES : Cheney & Kincaid : Numerical Mathematics and
*                                   Computing.
*                                  
*-
*   CALLS      :
*
*   WRITTEN BY : Vibeke Skytt, SI, 86-12.
*
************************************************************************
*/
{
  int kpos = 0;   /* Position of error.                               */
  int ki,kj,kk;   /* Counters.                                        */
  int kmax = 0;   /* Number of row with maximum greates element.      */
  int kchange;    /* Help variable to change order of two rows.       */
  double tmult;   /* Factor with which to multiply a row before it is 
		     added to another.                                */
  double t1;      /* Help variabel to find maximum of a number of elements.*/
  double tmax;    /* Maximum of a number of elements.                 */
  double tdiv;    /* Dividend in expression.                          */
  double *smax = SISL_NULL; /* Maximum elements in the rows of ea.         */
  
  /* Allocate space for local array.  */
  
  if ((smax = new0array(im,double)) == SISL_NULL) goto err101;
  
  /* Find largest element in each row.  */
  
  for (ki=0; ki<im; ki++)
    {                      
      nl[ki] = ki;
      for (kj=0; kj<im; kj++) 
	smax[ki] = MAX(smax[ki],fabs(ea[ki*im+kj]));
    }
  
  for (ki=0; ki<im-1; ki++)
    {
      
      /* Find row with maximum greates element to treat now.  */
      
      tmax = DZERO;  
      for (kj=ki; kj<im; kj++)
	{
	  tdiv = smax[nl[kj]];
	  if (DEQUAL(tdiv,DZERO)) goto warn1;
	  t1 = fabs(ea[nl[kj]*im+ki]/tdiv);
	  if (t1 > tmax)
	    {
	      tmax = t1; 
	      kmax = kj;
	    }
	}
      kchange  = nl[kmax];
      nl[kmax] = nl[ki];
      nl[ki]   = kchange;
      
      /* Add a multiplum of current row to all later rows in order to
	 get an upper triangular matrix ea.                           */
      
      for (kj=ki+1; kj<im; kj++)
	{
	  tdiv = ea[ki+kchange*im];
	  if (DEQUAL(tdiv,DZERO)) goto warn1;
	  tmult = ea[ki+nl[kj]*im]/tdiv;
	  ea[ki+nl[kj]*im] = tmult;
	  
	  for (kk=ki+1; kk<im; kk++)
	    ea[kk+nl[kj]*im] -= ea[kk+kchange*im]*tmult;
	}                                     
    }
  
  /* LU-factorizing performed.  */
  
  *jstat = 0;
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6lufacp",*jstat,kpos);
        goto out;

out:

/* Free space occupied by local array.  */

if (smax != SISL_NULL) freearray(smax);

return;
}

                                    

                                        
