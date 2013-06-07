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
 * $Id: s1602.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1602

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1602(double estapt[],double endpt[],int ik,int idim,double astpar,
	   double *cendpar,SISLCurve **rc,int *jstat)
#else
void s1602(estapt,endpt,ik,idim,astpar,cendpar,rc,jstat)
     double estapt[];
     double endpt[];
     int    ik;
     int    idim;
     double astpar;
     double *cendpar;
     SISLCurve  **rc;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Convert a straight line to a B-spline described curve.
*             
*
* INPUT      : estapt - start point of the straight line 
*              endpt  - end point of the straight line
*              ik     - the order of the B-spline curve to be found
*              idim   - The dimension of the space
*              astpar - start value of parameterization of the curve
*             
* OUTPUT
*            : cendpar - parameter used at the end of the curve
*              rc      - Pointer to the found curve
*              jstat  - status messages
  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : First the knots are found with ik knots in the start
*              and ik knots in the end of the curve. Then a modified
*              Marsdens identity is used to find ik vertices equally 
*              spaced on the straight line.
*
* EXAMPLE OF USE:
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6dist,newCurve,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 10. Nov 1988
*
*********************************************************************
*/
{
  int kit;            /* Loop control                                    */
  int kit2;           /* Loop contero                                    */
  int kvert;          /* Counter for position in vertex array            */
  int kpos=0;         /* Position of error                               */
  
  double *st=SISL_NULL;    /* Pointer to the first element of the knot vector
			 of the curve.                                   */
  double *scoef=SISL_NULL; /* Pointer to the first element of the curve's
			 B-spline coefficients.                          */
  double tdist;       /* Distance                                        */
  double tdel;        /* Delta x, y , ....                               */
  
  /* Check input */          
  
  if (idim <  1) goto err102;
  if (ik   <  2) goto err109;
  
  /* Find distance between start nd end point */
  tdist = s6dist(estapt,endpt,idim);
  
  
  /* Make knots. First allocate space */
  
  st = newarray(ik*2,DOUBLE);
  if (st == SISL_NULL) goto err101;
  
  for (kit=0; kit<ik; kit++) 
    {
      st[kit]    = astpar;
      st[kit+ik] = astpar + tdist;
    }
  
  /* calculate the vertices. First allocate space */ 
  
  /* First allocate space for vertices */ 
  
  scoef = newarray(ik*idim,DOUBLE);
  if (scoef == SISL_NULL) goto err101;
  
  /* Find first and last vertex. */ 
  
  kvert = (ik-1) * idim;
  for (kit=0; kit<idim; kit++,kvert++) 
    {
      scoef[kit]   = estapt[kit];
      scoef[kvert] = endpt[kit];
    }
  
  /* Find other vertices */ 
  
  for (kit=0; kit<idim; kit++)
    {   
      tdel = (endpt[kit] - estapt[kit])/(ik - 1);
      for (kit2=2; kit2<ik; kit2++)
	scoef[(kit2-1)*idim + kit] = scoef[(kit2-2)*idim + kit] + tdel; 
    }
  
  /* Make the curve */
  
  *rc = SISL_NULL;              
  *rc = newCurve(ik,ik,st,scoef,1,idim,1);
  if (*rc == SISL_NULL) goto err101;                
  
  *cendpar = st[ik];
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
 err101: 
  *jstat = -101;
  s6err("s1602",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension less than 1 */
  
 err102: 
  *jstat = -102;
  s6err("s1602",*jstat,kpos);
  goto out;                          
  /* Error in input. Order less than 2 */
  
 err109: 
  *jstat = -109;
  s6err("s1602",*jstat,kpos);
  goto out;                          
    
 out:
  if (st     != SISL_NULL) freearray(st);
  if (scoef  != SISL_NULL) freearray(scoef);
  return;
}          
                    
