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
 * $Id: s1364.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1364

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1364(SISLCurve *pc,double aepsge,int *jstat)
#else
void s1364(pc,aepsge,jstat)
     SISLCurve  *pc;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To decide if a B-spline curve is closed within a
*              tolerance
*
* INPUT      : pc     - The B-spline curve.   
*              aepsge - Geometric tolerance
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         = 1      : SISLCurve closed
*                                         = 0      : SISLCurve open
*                                         < 0      : error
*
* METHOD     : 
*
*
* REFERENCES :
*
*-                                                 
* CALLS      : s1221, s6dist, s6err
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdim;           /* Dimension of space                              */
  int kleft=0;        /* Pointer to knots                                */
  int kder=0;         /* Derivatives to be calculated                    */
  int kstat;          /* Local status variable                           */
  int kpos=0;         /* Position of error                               */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  double sdum1[3];    /* Arrays for calculation of points                */ 
  double sdum2[3];    /* Arrays for calculation of points                */
  double *sder1 = SISL_NULL; /* Pointers to points                            */
  double *sder2 = SISL_NULL; /* Pointers to points                            */
  double tdist;       /* Distance between points                         */
  
  /* Check if curve is correct */
  
  s1707(pc,&kstat);
  if (kstat<0) goto error;
  
  
  /* Copy curve attributes to local parameters.  */
  
  kn = pc -> in;
  kk = pc -> ik;
  kdim = pc -> idim;
  st = pc -> et;
  
  if (kdim>3)
    {
      sder1 = newarray(kdim,DOUBLE);
      sder2 = newarray(kdim,DOUBLE);
    }
  else
    {
      sder1 = sdum1;
      sder2 = sdum2;
    }
  
  /* Calculate start point of curve */
  
  s1221(pc,kder,st[kk-1],&kleft,sder1,&kstat);
  if (kstat<0) goto error;
  
  /* Calculate end point of curve */
  
  s1221(pc,kder,st[kn],&kleft,sder2,&kstat);
  if (kstat<0) goto error;
  
  tdist = s6dist(sder1,sder2,kdim);
  
  if (tdist>aepsge)
    *jstat = 0;
  else
    *jstat = 1;
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1364",*jstat,kpos);
  goto out;
 out:
  
  if (kdim>3)
    {
      if (sder1 != SISL_NULL) freearray(sder1);
      if (sder2 != SISL_NULL) freearray(sder2);
    }
  
  return;
}          
