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
 * $Id: s1840.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1840

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1840(SISLCurve *pcurve,double *cdist,int *jstat)
#else
void s1840(pcurve,cdist,jstat)
     SISLCurve  *pcurve;
     double *cdist;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the distance of the control polygon of a
*              B-spline curve to the control polygon of a straight line
*              from the first to the last vertex of the curve described
*              by the same control polygon.
*
*
* INPUT      : pcurve  - Pointer to surface object
*
*
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              *cdist - Deviation of curve polygon fromstraight line
*
* METHOD     : The knot vector is normalized to the interval [0,1]
*              Then Marsdens identity is used for representing the
*              function f(x)=x by the knot vector.
*
* REFERENCES : Theorem 4.21. (Marsdens identity). page 125. in
*              Larry L. Schumaker: Basic spline theory. 
*              John Wiley & Sons. 1981. ISBN 0-471-76475-2
*
*
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 10. Sept. 1988
*
*********************************************************************
*/
{
  int kstop,kj,ki,kp;        /* Loop control variable                     */
  int kpos = 0;              /* Position of errors                        */
  int kkm1;                  /* Orders minus 1                            */
  int kn,kk,kdim;            /* Variables in surface description          */
  double *st,*scoef;         /* Pointers to arrays describing surface     */
  double tstart,tlength;     /* Temproary variables                       */
  double tdiff;              /* Temproary variables                       */
  double *sm1=SISL_NULL;          /* Arrays for coefficients from Marsdens idnt*/
  double *start,*send;       /* Pointers to first and last vertex         */
  double tval0,tval1;        /* s and 1-s used in Marsdens identity       */
  double tsum;               /* Dummy variables                           */
  
  /* Initiate distance */
  
  *cdist = (double)0.0;         
  
  /* Make local versions of the curve object */
  
  st    = pcurve->et;
  scoef = pcurve->ecoef;
  kn    = pcurve->in;
  kk    = pcurve->ik;
  kdim  = pcurve->idim;
  
  /* Allocate space for arrays to be used by coefficients in Marsdens
   * identity, average of corners, diagonals, normal vector and
   * projection of corners into the plane defined by average of corners
   * and normalvector */
  
  sm1 = newarray(kn,double);
  if (sm1 == SISL_NULL) goto err101;
  
  /* Make representation of coefficients from Marsdens identity for the
     function f(t) = t, with the knot vector scaled to [0,1]. */
  
  tstart = st[kk-1];
  tlength = st[kn] - tstart;
  kkm1 = kk - 1;
  for (ki=0;ki<kn;ki++)
    {
      tsum = 0;
      kstop = ki+kk;
      for (kj=ki+1;kj<kstop;kj++)
        tsum +=st[kj];
      
      sm1[ki] = (tsum/kkm1-tstart)/tlength;
    }
  
  /* Calculate distances between the control polygon of the curve
     and the control polygon of the straight line */
  
  start = scoef;
  send  = scoef + kdim*(kn-1);
  kp    = 0;
  
  for (ki=0;ki<kn;ki++)
    {
      tval1 = sm1[ki];
      tval0 = (double)1.0 - tval1;
      tsum = (double)0.0; 
      
      for (kj=0;kj<kdim;kj++)
        {
	  
	  /* Accumulate distance between curve vertex and line vertex */
	  
	  tdiff = scoef[kp] - (tval0*start[kj] + tval1*send[kj]);
	  tsum += tdiff*tdiff;
	  kp++;
        }
      *cdist = MAX(*cdist,tsum);
    }
  
  *cdist = sqrt(*cdist);
  
  /* Everything ok */

  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
  err101: *jstat = -101;      
    s6err("s1840",*jstat,kpos);
    goto out;
  
  /* Free allocated space */
  
  out:
    if (sm1 != SISL_NULL) freearray(sm1);
    return;
}
