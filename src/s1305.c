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
 * $Id: s1305.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1305

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1305(double epar1[],double epar2[],double eval1[],double eval2[],
	   int *jbound,double gpar[],int *jstat)
#else
void s1305(epar1,epar2,eval1,eval2,jbound,gpar,jstat)
     double epar1[];
     double epar2[];
     double eval1[];
     double eval2[];
     int    *jbound;
     double gpar[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To find if there is an intersection between epar1 and
*              epar2 with the 2-D SISLbox desribed by eval1[0]:eval1[1]
*              in the first parameter direction and eval2[0]:eval2[1]
*              in the second parameter direction. If there is an
*              intersection find the intersection closest to the point that is
*              outside the area.
*
* INPUT      : epar1  - First parameter pair
*              epar2  - Second parameter pair
*              eval1  - Interval in first parameter direction
*              eval2  - Interval in second parameter direction
*
*
* OUTPUT     : gpar   - Parameter pair of intersection
*              jbound - Indicator telling along which boundary
*                       we have an intersection
*                       = 0      :  no intersection
*                       = 1      : intersection along u=eval1[0]
*                       = 2      : intersection along v=eval2[1]
*                       = 3      : intersection along u=eval1[1]
*                       = 4      : intersection along v=eval2[0]
*              jstat  - status messages  
*                       = 0      : Line outside no intersection
*                       = 1      : Line inside  no intersection
*                       = 2      : epar2 outside, epar1 inside, step out
*                       = 3      : epar1 outside, epar2 inside, step in
*                       = 4      : epar2 outside, epar1 on boundary
*                       < 0      : error         
*
*                                  
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 4-August-1988
* Revised by : Tor Dokken, SI, Oslo, Norway, Mars-1989
*              Improved clipping of long steps               
* Revised by : UJK, SI, Oslo, Norway, Oct-1992
*              Noise treatment.
*********************************************************************
*/
{
  int    kstat;          /* Local status variable                            */
  int    kins1;          /* epar1 inside/outside SISLbox                         */
  int    kins2;          /* epar2 inside/outside SISLbox                         */
  double tdom;           /* Denominator of last intersection point           */
  double tfak1,tfak2;    /* Distance to straight line                        */
  double tcdist=HUGE;    /* Current distance                                 */
  double tpdist=HUGE;    /* Previous distance                                */
  double simpli[5];      /* Corners of the parameter area pu into the
			    implicit equation of the line through epar1 and
			    epar2.                                           */
  double snorm[2];       /* Normal vector                                    */
  double stemp[2];       /* Candidate point                                  */
  double *outside=epar1; /* Pointer to the point that is outside;            */
  *jbound = 0;
  
  /* Test if both ends are inside */
  
  kins1 = kins2 = 0; 
  
  if (eval1[0] <= epar1[0] && epar1[0] <= eval1[1] &&
      eval2[0] <= epar1[1] && epar1[1] <= eval2[1]) kins1 = 1;
  
  if( eval1[0] <= epar2[0] && epar2[0] <= eval1[1] &&
     eval2[0] <= epar2[1] && epar2[1] <= eval2[1]) kins2 = 1;
  
  if (kins1==1 && kins2==1) goto war01;
  
  if (kins1) outside = epar2;
  
  /* Test if we step from the boundary and out */
  
  if ((eval1[0] == epar1[0] && epar2[0] < eval1[0]) ||
      (epar1[0] == eval1[1] && eval1[1] < epar2[0]) ||
      (eval2[0] == epar1[1] && epar2[1] < eval2[0]) ||
      (epar1[1] == eval2[1] && eval2[1] < epar2[1])    ) goto war04;
  
  /* Test if both ends are to the left, right, below or above */
  
  if ( (epar1[0] < eval1[0] && epar2[0] < eval1[0]) ||
      (eval1[1] < epar1[0] && eval1[1] < epar2[0]) ||
      (epar1[1] < eval2[0] && epar2[1] < eval2[0]) ||
      (eval2[1] < epar1[1] && eval2[1] < epar2[1])   ) goto war00;
  
  /* Make normal vector of line though epar1 and epar2 */
  
  snorm[0] = -(epar2[1] - epar1[1]);
  snorm[1] =   epar2[0] - epar1[0];
  
  (void)s6norm(snorm,2,snorm,&kstat);
  
  /* Put corners of parameter area into the implicit equation of the straight
     line */
  
  simpli[0] = (eval1[0]-epar1[0])*snorm[0] + (eval2[0]-epar1[1])*snorm[1];  
  simpli[1] = (eval1[0]-epar1[0])*snorm[0] + (eval2[1]-epar1[1])*snorm[1];  
  simpli[2] = (eval1[1]-epar1[0])*snorm[0] + (eval2[1]-epar1[1])*snorm[1];  
  simpli[3] = (eval1[1]-epar1[0])*snorm[0] + (eval2[0]-epar1[1])*snorm[1];  
  simpli[4] = simpli[0];
  
  /* If simpli[0:3] all have the same sign, the straight line is outside */
  
  if ((simpli[0]>(double)0.0 && simpli[1]>(double)0.0 && 
       simpli[2]>(double)0.0 && simpli[3]>(double)0.0) ||
      (simpli[0]<(double)0.0 && simpli[1]<(double)0.0 && 
       simpli[2]<(double)0.0 && simpli[3]<(double)0.0)   ) goto war00;
  
  /* Treate intersections with left boundary */
  
  if (simpli[0]*simpli[1] <= (double)0.0 && epar1[0] != eval1[0])
    {
      /*  Intersection along left boundary */
      
      tfak1 = fabs(simpli[0]);
      tfak2 = fabs(simpli[1]);
      tdom  = tfak1 + tfak2;
      
      if (DNEQUAL(tdom,(double)0.0))
        {
	  /* The straight line and the left boundary does not coinside */
	  
	  stemp[0] = eval1[0];
	  stemp[1] = (tfak2*eval2[0] + tfak1*eval2[1])/tdom;
	  tcdist   = s6dist(stemp, outside,2);
	  if (*jbound == 0 || tcdist < tpdist)
            {
	      /* New point closer than previous intersection point */
	      
	      gpar[0] = stemp[0];
	      gpar[1] = stemp[1];
	      *jbound  = 1;
	      tpdist  = tcdist;
            }
        }
    }
  
  /* Treate intersections with upper boundary */                       
  
  if (simpli[1]*simpli[2] <= (double)0.0 && epar1[1] != eval2[1])
    {
      /* Intersection along upper boundary */
      
      tfak1 = fabs(simpli[1]);
      tfak2 = fabs(simpli[2]);
      tdom  = tfak1 + tfak2;
      if (DNEQUAL(tdom,(double)0.0))
        {
	  /* The straight line and the upper boundary does not coinside */
	  
	  stemp[0] = (tfak2*eval1[0] + tfak1*eval1[1])/tdom;
	  stemp[1] = eval2[1];
	  tcdist   = s6dist(stemp, outside,2);
	  
	  if (*jbound == 0 || tcdist < tpdist)
            {
	      /* New point closer than previous intersection point */   
	      
	      gpar[0] = stemp[0];
	      gpar[1] = stemp[1];
	      *jbound  = 2;
	      tpdist  = tcdist;
            }
        }
    }
  
  /* Treate intersections with right boundary */
  
  if (simpli[2]*simpli[3] <= (double)0.0 && epar1[0] != eval1[1])
    {
      /* Intersection along right boundary */
      
      tfak1 = fabs(simpli[2]);
      tfak2 = fabs(simpli[3]);
      tdom  = tfak1 + tfak2;
      if (DNEQUAL(tdom,(double)0.0))
        {
	  /* The straight line and the right boundary does not coinside */
	  
	  stemp[0] = eval1[1];
	  stemp[1] = (tfak2*eval2[1] + tfak1*eval2[0])/tdom;
	  tcdist   = s6dist(stemp, outside,2);
	  
	  if (*jbound == 0 || tcdist < tpdist)
            {
	      /* New point closer than previous intersection point */
	      
	      gpar[0] = stemp[0];
	      gpar[1] = stemp[1];
	      *jbound  = 3;
	      tpdist  = tcdist;
            }
        }
    }
  
  /* Treate intersections with lower boundary */
  
  if (simpli[3]*simpli[4] <= (double)0.0 && epar1[1] != eval2[0])
    {
      /* Intersection along lower boundary */
      
      tfak1 = fabs(simpli[3]);
      tfak2 = fabs(simpli[4]);
      tdom  = tfak1 + tfak2;
      if (DNEQUAL(tdom,(double)0.0))
        {
	  /* The straight line and the lower boundary does not coinside */
	  
	  stemp[0] = (tfak2*eval1[1] + tfak1*eval1[0])/tdom;
	  stemp[1] = eval2[0];
	  tcdist   = s6dist(stemp, outside,2);
	  
	  if (*jbound == 0 || tcdist < tpdist)
            {
	      /* New point closer than previous intersection point */   
	      
	      gpar[0] = stemp[0];
	      gpar[1] = stemp[1];
	      *jbound = 4;
	      tpdist  = tcdist;
            }
        }
    }
  
  if (kins1 == 1)
    goto war02; 
  
  if (kins2 == 1 || *jbound != 0) goto war03; 
  
  goto war05;
  
  /* Line outside */
  
 war00:
  *jstat = 0;
  goto out;                                                                      
  
  /* Line inside */                          
  
 war01:
  *jstat = 1;
  goto out;
  
  /* epar1 inside epar2 outside */
 war02:
  *jstat = 2;
  goto out;
  
  /* epar2 inside epar1 outside */
 war03:
  *jstat = 3;
  goto out;
  
  /* epar1 on boundary, epar2 outside */
 war04:
  *jstat = 4;
  goto out;
  
  /* Special error */
 war05:
  *jstat = 5;
  goto out;
  
 out:
  return;
}
