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
 * $Id: s6decomp.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DECOMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6decomp(double ea[],double gx[],double eb1[],double eb2[],double eb3[],int *jstat)
#else
void s6decomp(ea,gx,eb1,eb2,eb3,jstat)
     double ea[];
     double gx[];
     double eb1[];
     double eb2[];
     double eb3[];
     int    *jstat;
#endif
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : In dimention tree we change basis from a euclides bases
*             to eb1, eb2, eb3.
*
*
*   INPUT   : ea   - The orginale vector.
*             eb1  - the first bases vector.
*             eb2  - the first bases vector.
*             eb3  - the first bases vector.
*
*
*   
*   OUTPUT  : gx    - The new vector
*             jstat - Status variable.
*                       < 0 : error
*                       = 0 : ok
*                       > 0 : warning
*
*
*   METHOD  : 
*
*
*   REFERENCES : 
*-
*   CALLS      :
*
*   WRITTEN BY : Arne Laksaa, SI, 89-07.
*
************************************************************************
*/
{
  int kstat =0;       /* Local status variable.    */
  int ki;             /* Counter.                  */
  int n1[3];          /* Array for use in lufac.   */
  double sc[9],se[3]; /* Matrix and help vector.   */
  
  
  /* Copy new bases into local matrix.  */
  
  memcopy(sc,eb1,3,double);
  memcopy(sc+3,eb2,3,double);
  memcopy(sc+6,eb3,3,double);
  
  
  s6lufacp(sc,n1,3,&kstat);
  if (kstat < 0) goto error;
  else if (kstat > 0) goto warn1;
  
  for (ki=0; ki<3; ki++)
    {                      
      se[0] = se[1] = se[2] = DZERO;
      se[ki] = (double)1;
      
      s6lusolp(sc,se,n1,3,&kstat);
      if (kstat < 0) goto error;
      else if (kstat > 0) goto warn1;
      
      gx[ki] = s6scpr(ea,se,3);
    }
  
  /* Change of bases performed.  */
  
  *jstat = 0;
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in subrutines.  */

error: *jstat = kstat;
        s6err("s6decomp",*jstat,0);
        goto out;

 out: ;
}

                                    

                                        
