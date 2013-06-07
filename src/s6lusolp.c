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
 * $Id: s6lusolp.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6LUSOLP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6lusolp(double ea[],double eb[],int nl[],int im,int *jstat)
#else
void s6lusolp(ea,eb,nl,im,jstat)
     double ea[];
     double eb[];
     int    nl[];
     int    im;
     int    *jstat;
#endif
/*
************************************************************************
*
***********************************************************************
*
*   PURPOSE : Solve the equationsystem LU=eb by forth- and back-
*             substitution. L and U are stored in ea.
*
*
*   INPUT   : ea   - The factorized coeffecient-matrix.
*             im   - The number of equations.
*             nl   - Ordering of lines in the matrix.
*
*
*   INPUT/OUTPUT : eb - the right side of the equationsystem and
*                       the found values for the unknowns.
*
*   OUTPUT  : jstat  - Status variable.
*                        < 0 : Error
*                        = 0 : ok
*                        > 0 : Warning
*
*                                                                       
*   METHOD  : Solve on the equation-system LU=eb.
*
*
*   REFERENCES : Cheney & Kincaid : Numerical Mathematics and
*                                   Computing.
*
*-
*   CALLS      :
*
*   WRITTEN BY : Vibeke Skytt, SI, 86-10.
*
************************************************************************
*/
{
  int kpos = 0;      /* Position of error.                             */
  int ki,kj;         /* Counters.                                      */
  double *sx = SISL_NULL; /* Array used to keep solution of equation system
			internally.                                    */
  double tdiv;       /* Dividend in expression.                        */
  
  /* Allocate space for local array.  */
  
  if ((sx = newarray(im,double)) == SISL_NULL) goto err101;
  
  for (ki=0; ki<im-1; ki++)
    {
      /*  Gauss on right side of equation  */
      
      for (kj=ki+1; kj<im; kj++)      
	eb[nl[kj]] -= eb[nl[ki]]*ea[ki+nl[kj]*im];
    }
  
  tdiv = ea[im-1+nl[im-1]*im];
  if (DEQUAL(tdiv,DZERO)) goto warn1;
  sx[im-1] = eb[nl[im-1]]/tdiv;
  
  for (ki=im-2; ki>=0; ki--)
    {
      /*  Backwards substitution.   */
      
      for (kj=ki+1; kj<im; kj++)
	eb[nl[ki]] -= sx[kj]*ea[kj+nl[ki]*im];
      
      tdiv = ea[ki+nl[ki]*im];
      if (DEQUAL(tdiv,DZERO)) goto warn1;
      sx[ki] = eb[nl[ki]]/tdiv;
    }   
  for (ki=0; ki<im; ki++) eb[ki] = sx[ki];
  
  /* Equation system solved.  */
  
  *jstat = 0; 
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6lusolp",*jstat,kpos);
        goto out;

out:

/* Free space occupied by local array.  */

if (sx != SISL_NULL) freearray(sx);

return;
}
