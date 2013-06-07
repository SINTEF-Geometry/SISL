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
 * $Id: s1240.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1240

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1240(SISLCurve *pcurve,double aepsge,double *clength,int *jstat)
#else
void s1240(pcurve,aepsge,clength,jstat)
     SISLCurve  *pcurve;
     double aepsge;
     double *clength;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Calculate the length of a B-spline curve. The length
*              calculated will not deviate more than (aepsge/the
*              length calculated) from the real length of the curve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : clength - The length of the curve.
*              jstat   - status messages  
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
*-
* CALLS      : s6dist,s1251,make_cv_kreg.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;  /* Local status variable.                          */
  int kpos = 0;   /* Position of error.                              */
  int ki;         /* Counter.                                        */
  int kdim;       /* Dimension of the space in which the curve lies. */
  int kn;         /* Number of vertices of curve.                    */
  int kcalc;      /* Indicates if correct length of curve is found.  */
  double tlength; /* Length of curve.                                */
  double tprev;   /* Previous length of curve calculated.            */
  double teps;    /* Local tolerance.                                */
  double *s1;     /* Pointer used to traverse real array.            */
  SISLCurve *qc=SISL_NULL;  /* k-regular local curve.                     */
  
  if (pcurve->cuopen == SISL_CRV_PERIODIC)
    {
       /* Make curve k-regular. */
       
       make_cv_kreg(pcurve,&qc,&kstat);
       if (kstat < 0) goto error;
    }
  else qc = pcurve;
       
  /* Copy curve information to local parameters. */
  
  kdim = qc -> idim;
  kn   = qc -> in;
  
  /* Calculate length of control polygon.  */
  
  tlength = 0;
  for (ki=1,s1=qc->ecoef+kdim; ki<kn; ki++,s1+=kdim)
    tlength += s6dist(s1-kdim,s1,kdim);
  
  /* Set up local tolerance.  */
  
  teps = aepsge*100;
  
  kcalc = 0;
  while (kcalc == 0)
    {
      teps = teps/2.0;
      tprev = tlength;
      
      /* Compute length of curve.  */
      
      s1251(qc,teps,&tlength,&kstat);
      if (kstat < 0) goto error;
      
      /* Test if the error is within the tolerance. */
      
      if (fabs(tprev-tlength)/MAX(tprev,tlength) < aepsge) kcalc = 1;
      
    }
  
  /* Length of curve calculated. */
  
  *clength = tlength;
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine. */
  
  error : *jstat = kstat;
  s6err("s1240",*jstat,kpos);
  goto out;
  
 out: 
    
    /* Free local curve.  */
    
    if (qc != SISL_NULL && qc != pcurve) freeCurve(qc);
    
    return;
}

