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
 * $Id: s1606.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1606

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1606 (SISLCurve * pc1, SISLCurve * pc2, double aepsge, double epoint1[], double epoint2[],
       int itype, int idim, int ik, SISLCurve ** rc, int *jstat)
#else
void
s1606 (pc1, pc2, aepsge, epoint1, epoint2, itype, idim, ik, rc, jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     double aepsge;
     double epoint1[];
     double epoint2[];
     int itype;
     int idim;
     int ik;
     SISLCurve **rc;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate a blending curve between two curves.
*              Points indicate between which ends the blend is to be
*              produced. the blending curve produced is a conic section
*              if possible, else it is a quadratic polynomial spline curve.
*
* INPUT      : pc1    - The first input B-spline curve.
*              pc2    - The second input B-spline curve.
*              aepsge - Geometry resolution.
*              epnt1  - SISLPoint near the end of curve 1 where the blend starts.
*              epnt2  - SISLPoint near the end of curve 2 where the blend starts.
*              itype  - Type of blending:
*                     = 1  - circle, interpolating tangent on first curve,
*                       not on curve 2.
*                     = 2  - conic if possible
*                     else - polynomial segment
*              idim   - Dimension of space.
*              ik     - order of blending curve.
*
* OUTPUT     : rc     - Blending curve produced
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : The input points are used to find between which ends
*              the blending curve is to be calculated. Then the
*              parameter values of these ends are used for calculating
*              the end points and tangents. The next step is to interpolate
*              these points and tangents with the curve type specified
*              by the blend type indicator. If a conic blend is specified
*              and can not be produced, then a quadratic blend is calculated
*              instead.
*
*
* USE        :
*
* REFERENCES :
*
* CALLS      : s1707,s1607,s6err.
*
* WRITTEN BY : Christophe Birkeland, SI, 1991-07
* REVISED BY : Johannes Kaasa, SI, Aug. 92 (Removed call to s1707, this is
*              called from s1607).
*
*********************************************************************
*/
{
  int ki;			/* Loop control variable			   */
  int kstat = 0;		/* Status variable                                 */
  int kpos = 0;			/* Position of error                               */
  int in1;			/* The number of B-splines, i.e., the dimension of
				   the spline space associated with the knot
				   vector. Curve 1 and 2.			   */
  int in2;
  int ik1;			/* The polynomial order of the curves.             */
  int ik2;
  double fil1, end1;		/* Used to store start and end of curve 1	   */
  double fil2, end2;		/* Used to store start and end of curve 2	   */
  double sum11, sum12;		/* Variables used to find end od curves		   */
  double sum21, sum22;
  double dum;			/* Dummy variable		                   */

  *jstat = 0;


  /* Check if conflicting dimensions. */

  if (pc1->idim != idim || pc2->idim != idim)
    goto err106;

  /* Find ends of curve. */

  sum11 = (double) 0.0;
  sum12 = (double) 0.0;
  sum21 = (double) 0.0;
  sum22 = (double) 0.0;

  in1 = pc1->in -1;
  in2 = pc2->in -1;

  for (ki = 0; ki < idim; ki++)
    {
      dum = epoint1[ki] - pc1->ecoef[ki];
      sum11 += dum * dum;
      dum = epoint1[ki] - pc1->ecoef[in1 * idim + ki];
      sum12 += dum * dum;
      dum = epoint2[ki] - pc2->ecoef[ki];
      sum21 += dum * dum;
      dum = epoint2[ki] - pc2->ecoef[in2 * idim + ki];
      sum22 += dum * dum;
    }
  in1++;
  in2++;

  ik1 = pc1->ik - 1;
  ik2 = pc2->ik - 1;
  if (sum11 < sum12)
    {
      /* Start of curve 1 closer than end of curve 1. */

      fil1 = pc1->et[ik1];
      end1 = pc1->et[in1];
    }
  else
    {
      /* End of curve 1 closer than start of curve 1. */

      end1 = pc1->et[ik1];
      fil1 = pc1->et[in1];
    }
  if (sum21 < sum22)
    {
      /* Start of curve 2 closer than end of curve 2. */

      fil2 = pc2->et[ik2];
      end2 = pc2->et[in2];
    }
  else
    {
      /* End of curve 2 closer than start of curve 2. */

      end2 = pc2->et[ik2];
      fil2 = pc2->et[in2];
    }

  /* Make blend */

  s1607 (pc1, pc2, aepsge, end1, fil1, end2, fil2, itype, idim, ik,
	 rc, &kstat);
  if (kstat < 0)
    goto error;

  goto out;


  /* Conflicting dimensions */

err106:
  *jstat = -106;
  s6err ("s1606", *jstat, kpos);
  goto out;

  /* Error in lower level function */

error:
  *jstat = kstat;
  s6err ("s1606", *jstat, kpos);
  goto out;

out:
  return;
}
