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
 * $Id: s6takunion.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6TAKEUNION

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s6takeunion(double evec1[],int ielem1,double evec2[],int ielem2,
		 double **gunion,int *jnmbelem,int *jstat)
#else	 
void s6takeunion(evec1,ielem1,evec2,ielem2,gunion,jnmbelem,jstat)
     int ielem1,ielem2,*jnmbelem,*jstat;
     double evec1[],evec2[],**gunion;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Take the union between two sorted double vectors. Identical 
*              elements in the two vectors are only represented once. 
*              If several elements in one array are equal, all are
*              represented.
*
*
* INPUT      : evec1    - First vector.
*              ielem1   - Number of elements of evec1.
*              evec2    - Second vector.
*              ielem2   - Number of elements of evec2.
*                       
*
* OUTPUT     : gunion   - Union vector.
*              jnmbelem - Number of elements of gunion.
*              jstat    - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int knelem;
  int knunion;
  double *sunion = SISL_NULL;
  double *s1,*s1stop;
  double *s2,*s2stop;
  
  /* Make local array to store the union of the vectors.  */

  knelem = ielem1 + ielem2;
  if ((sunion = newarray(knelem,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Produce union vector. */

  for (s1=evec1,s1stop=s1+ielem1,s2=evec2,s2stop=s2+ielem2,knunion=0; 
       s1<s1stop && s2<s2stop;)
      {
	if (*s1 < *s2)
	  sunion[knunion++] = *s1++;
	else if (*s2 < *s1)
	  sunion[knunion++] = *s2++;
	else
	  {
	    sunion[knunion++] = *s1++;
	    s2++;
	  }
      }

  for (; s1<s1stop; s1++,knunion++)
    sunion[knunion] = *s1;

  for (; s2<s2stop; s2++,knunion++)
    sunion[knunion] = *s2;
  
  /* Allocate scratch for output union vector.  */

  *gunion = SISL_NULL;
  if ((*gunion = newarray(knunion,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Copy union vector to output vector.  */

  memcopy(*gunion,sunion,knunion,DOUBLE);
  *jnmbelem = knunion;
  
  /* Union found.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  out :
    /* Free scratch occupied by local array.  */

    if (sunion != SISL_NULL) freearray(sunion);
  
  return;
}
