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

#define SH1831

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     sh1831(SISLCurve *pc1, SISLCurve *pc2, int isign, double epoint[], 
	    double enorm[], double aepsge, int *jstat)
#else
void sh1831(pc1, pc2, isign, epoint, enorm, aepsge, jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     int isign;
     double epoint[];
     double enorm[];
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Interception test between two curves. Test if it is
*              possible to split the curves by a given plane.
*
*
*
* INPUT      : pc1    - First curve.
*              pc2    - Second curve.
*              epoint - Point in plane.
*              enorm  - Plane normal.
*              aepsge - Geometric tolerance.
*
*
*
* OUTPUT     : jstat  - status messages  
*                              = 1      : Possibility of inner intersections.
*                              = 0      : No possibility of inner intersections.
*                              < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6scpr - Scalar product between two vectors.
*              s6diff - Difference vector between two vectors.
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 94-02.
*
*********************************************************************
*/
{
  int kpos = 0;          /* Position of error.               */
  int ki;
  int kdim;              /* Dimension of space.              */
  int kbez1, kbez2;      /* Indicates if the curves are of type Bezier. */
  int ksignprev = 0;     /* Sign of distance between previous curve and plane.*/
  int ksign1 = 0;        /* Sign of distance between curve and plane.*/
  int ksign2 = 0;        /* Sign of distance between curve and plane.*/
  double tdist;          /* Distance between coefficient and plane.     */
  double *s1;            /* Pointer to coefficient of curve. */
  double sdiff[3];       /* Difference vector.               */
  
  /* Test dimension of geometry space. */
  
  kdim = pc1->idim;
  if (kdim != 2 && kdim != 3) goto err105;
  if (kdim != pc2->idim) goto err106;
  
  /* Test if the curves are Bezier curves. */
  
  kbez1 = (pc1->ik == pc1->in) ? 1 : 0;
  kbez2 = (pc2->ik == pc2->in) ? 1 : 0;
  
  /* For each curve, compute the distance between the coefficients of the
     curve and the given plane.    */
  
  for (s1=pc1->ecoef, ki=0; ki<pc1->in; ki++, s1+=kdim)
  {
     s6diff(epoint, s1, kdim, sdiff);
     tdist = s6scpr(sdiff, enorm, kdim);
     
     if (fabs(tdist) <= aepsge && !kbez1 && !(ki==0 || ki==pc1->in-1)) break;
     ksign2 = (DEQUAL(tdist,DZERO)) ? 0 : ((tdist > 0) ? 1 : -1);
     if (ksign1*ksign2 < 0) break;
     ksign1 = ksign2;
  }
  if (ki < pc1->in)
  {
     *jstat = 1;
     goto out;
  }

  ksignprev = isign*ksign1;
  ksign1 = 0;
  for (s1=pc2->ecoef, ki=0; ki<pc2->in; ki++, s1+=kdim)
  {
     s6diff(epoint, s1, kdim, sdiff);
     tdist = s6scpr(sdiff, enorm, kdim);
     
     if (fabs(tdist) <= aepsge && !kbez2 && !(ki==0 || ki==pc2->in-1)) break;
     ksign2 = (DEQUAL(tdist,DZERO)) ? 0 : ((tdist > 0) ? 1 : -1);
     if (ksign1*ksign2 < 0) break;
     if (ksignprev*ksign1 > 0) break;
     ksign1 = ksign2;
  }
  if (ki < pc2->in)
  {
     *jstat = 1;
     goto out;
  }
  
  goto out;
  
  /* Error in input. Dimension not equal to 2 or 3.  */
  
 err105: *jstat = -105;
  s6err("sh1831",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("sh1831",*jstat,kpos);
  goto out;
  
 out:
  
  return;
}
