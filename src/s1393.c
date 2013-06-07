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
 * $Id: s1393.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1393

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1393(int n1,SISLCurve *pc1[],SISLCurve *sc1[],SISLCurve *ec1[],int *jstat)
#else
void s1393(n1,pc1,sc1,ec1,jstat)
     int   n1;
     SISLCurve *pc1[];
     SISLCurve *sc1[];
     SISLCurve *ec1[];
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* Purpose :  Split curve at midpoint, turn second curve, 
*            and normalize parameterinterval.
*
* Input     : pc1       - Pointers to first boundary curves
*             n1        - Number of curves
*
* Output    : sc1       - Pointers to first part of curve.
*             ec1       - Pointers to second part of curve.
*
*             jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                      
*-
* Calls      : s6err - error messages.
*              s1710 - split curve at given parameter value.
*              s1706 - Turn parametrization of curve.
*              s1399 - Normalize parameterinterval.
*
* Written by : Mortend Daehlen, SI, Aug. 88.
*
*********************************************************************
*/                                     
{
  int kpos = 0;
  int ki;
  int kstat = 0;
  double ax,astart,astop;
  SISLCurve *h1,*h2;
  
  astart = DZERO;
  astop  = (double)1.0;
  
  /* For each curve in pc1 split/turn and normalize. */
  for (ki=0;ki<n1;ki++)
    {
      
      /* Split */
      
      ax=(pc1[ki]->et[pc1[ki]->in]-(pc1[ki]->et[(pc1[ki]->ik)-1]))/(double)2.0;
      s1710(pc1[ki],ax,&h1,&h2,&kstat); 
      if (kstat < 0) goto error;
      
      /* Turn */
      
      s1706(h2);
      if (kstat < 0) goto error;
      
      /* Normalize */
      
      s1399(h1,astart,astop);
      if (kstat < 0) goto error;
      s1399(h2,astart,astop);
      if (kstat < 0) goto error;
      sc1[ki] = h1;
      ec1[ki] = h2;
    }
  
  *jstat=0;
  goto out;
  
  /* Error in lower level routine.   */
  
 error: *jstat = kstat;
  s6err("s1393",*jstat,kpos);
  goto out;
  
 out: return;
}
