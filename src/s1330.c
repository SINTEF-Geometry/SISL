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
 * $Id: s1330.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1330

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1330(double epar11[],double epar12[],double epar21[],double epar22[],
	   double eval11[],double eval12[],double eval21[],double eval22[],
	   int *jbound,double gpar1[],double gpar2[],int *jstat)
#else
void s1330(epar11,epar12,epar21,epar22,eval11,eval12,eval21,eval22,
           jbound,gpar1,gpar2,jstat)
     double epar11[];
     double epar12[];
     double epar21[];
     double epar22[];
     double eval11[];
     double eval12[];
     double eval21[];
     double eval22[];
     int    *jbound;
     double gpar1[];
     double gpar2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To find if there is an intersection between epar1 and
*              epar2 with the 4-D SISLbox desribed by eval11[0]:eval11[1]
*              in the first parameter direction and eval12[0]:eval12[1]
*              in the second parameter direction in the first patch, and
*              eval21[0]:eval21[1] in the first parameter direction and
*              eval22[0]:eval22[1] in the second parameter direction in
*              the second patch. If there is an intersection find the
*              intersection closest to epar1.
*
* INPUT      : epar11 - Start of line in first surface
*              epar12 - Start of line in second surface
*              epar21 - End of line in first surface
*              epar22 - End of line in second surface
*              eval11 - Interval in first parameter direction in patch 1
*              eval12 - Interval in second parameter direction in patch 1
*              eval21 - Interval in first parameter direction in patch 2
*              eval22 - Interval in second parameter direction in patch 2
*
*
* OUTPUT     : gpar1  - Parameter pair of intersection in first surface
*            : gpar2  - Parameter pair of intersection in second surface
*              jbound - Indicator telling along which boundary
*                       we have an intersection
*                       = 0      :  no intersection
*                       = 1      : intersection along u=eval11[0]
*                       = 2      : intersection along v=eval12[1]
*                       = 3      : intersection along u=eval11[1]
*                       = 4      : intersection along v=eval12[0]
*                       = 5      : intersection along s=eval21[0]
*                       = 6      : intersection along t=eval22[1]
*                       = 7      : intersection along s=eval21[1]
*                       = 8      : intersection along t=eval22[0]
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
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 16-August-1988
* Revised by : Tor Dokken, SI, Oslo, Norway, 3-april-1989
*              Correction of long steps
*
*********************************************************************
*/
{
  int    kstat1,kstat2;  /* Local status variable                            */
  int    kstat;          /* Local status variable                            */
  int    kpos=0;         /* Position of error                                */
  int    kins1;          /* epar1 inside/outside SISLbox                         */
  int    kins2;          /* epar2 inside/outside SISLbox                         */
  int    kbound1;        /*Intersection indicator along boundary of surface 1*/
  int    kbound2;        /*Intersection indicator along boundary of surface 2*/
  
  double tdom;           /* Denominator of last intersection point           */
  double tfak1,tfak2;    /* Distance to straight line                        */
  double spar11[2];      /*Candidate intersection point two first coordinates*/
  double spar12[2];      /*Candidate intersection point two last coordinates */
  double spar21[2];      /*Candidate intersection point two first coordinates*/
  double spar22[2];      /*Candidate intersection point two last coordinates */
  
  *jbound = 0;
  
  /* Test if both ends are inside */
  
  kins1 = kins2 = 0; 
  
  if (eval11[0] <= epar11[0] && epar11[0] <= eval11[1] &&
      eval12[0] <= epar11[1] && epar11[1] <= eval12[1] &&
      eval21[0] <= epar12[0] && epar12[0] <= eval21[1] &&
      eval22[0] <= epar12[1] && epar12[1] <= eval22[1]) kins1 = 1;
  
  if (eval11[0] <= epar21[0] && epar21[0] <= eval11[1] &&
      eval12[0] <= epar21[1] && epar21[1] <= eval12[1] &&
      eval21[0] <= epar22[0] && epar22[0] <= eval21[1] &&
      eval22[0] <= epar22[1] && epar22[1] <= eval22[1]) kins2 = 1;
  
  
  /* Test if we step from the boundary and out */
  
  if ((eval11[0] == epar11[0] && epar21[0] < eval11[0]) ||
      (epar11[0] == eval11[1] && eval11[1] < epar21[0]) ||
      (eval12[0] == epar11[1] && epar21[1] < eval12[0]) ||
      (epar11[1] == eval12[1] && eval12[1] < epar21[1]) ||
      (eval21[0] == epar12[0] && epar22[0] < eval21[0]) ||
      (epar12[0] == eval21[1] && eval21[1] < epar22[0]) ||
      (eval22[0] == epar12[1] && epar22[1] < eval22[0]) ||
      (epar12[1] == eval22[1] && eval22[1] < epar22[1])) goto war04;
  
  if (kins1==1 && kins2==1) goto war01;
  
  /* Test if both ends are to the left, right, below or above */
  
  if ( (epar11[0]  < eval11[0] && epar21[0]  < eval11[0]) ||
      (eval11[1] < epar11[0]  && eval11[1] < epar21[0] ) ||
      (epar11[1]  < eval12[0] && epar21[1]  < eval12[0]) ||
      (eval12[1] < epar11[1]  && eval12[1] < epar21[1] ) ||
      (epar12[0]  < eval21[0] && epar22[0]  < eval21[0]) ||
      (eval21[1] < epar12[0]  && eval21[1] < epar22[0] ) ||
      (epar12[1]  < eval22[0] && epar22[1]  < eval22[0]) ||
      (eval22[1] < epar12[1]  && eval22[1] < epar22[1] )   ) goto war00;
  
  
  
  /* Check if intersection in first two dimensions */                            
  
  s1305(epar11,epar21,eval11,eval12,&kbound1,spar11,&kstat);
  
  if (kstat<0) goto error;
  kstat1 = kstat;
  if (kstat1==0) goto war00;
  
  /* Calculate two last coefficients */
  
  if (kstat1==2 || kstat1==3)
    {
      tfak1 = fabs(spar11[0]-epar11[0]) + fabs(spar11[1]-epar11[1]);
      tfak2 = fabs(epar21[0]-spar11[0]) + fabs(epar21[1]-spar11[1]);
      tdom = tfak1 + tfak2;
      if (DNEQUAL(tdom,DZERO))
        {
	  spar12[0] = (tfak2*epar12[0] + tfak1*epar22[0])/tdom;
	  spar12[1] = (tfak2*epar12[1] + tfak1*epar22[1])/tdom;
	  
	  /* If the two last coefficients are zero, then then this intersection
	     point must be discarded */
	  
	  if (spar12[0]<eval21[0] || eval21[1]<spar12[0] ||
	      spar12[1]<eval22[0] || eval22[1]<spar12[1])
            {
	      /* Intersection point outside */
	      
	      kbound1 = 0;
            }
        }
      else
        {
	  /* epar1, spar and epar2 has equal first coordinates, since there
	     is an intersection all must lie on the boundary, thus all are
	     inside */
	  
	  kbound1 = 0;
        }
    }
  else if (kstat1==4 && kins1==1)
    {
      /* On boundary and stepping out */
      goto war04;
    }
  
  /* Check if intersection in last two dimensions */
  
  s1305(epar12,epar22,eval21,eval22,&kbound2,spar22,&kstat);
  
  if (kstat<0) goto error;
  kstat2 = kstat;
  if (kstat2==0) goto war00;
  
  if (kstat1==1 && kstat2==1) goto war01;
  
  
  /* Calculate two last coefficients */
  
  if (kstat2==2 || kstat2==3)
    {
      tfak1 = fabs(spar22[0]-epar12[0]) + fabs(spar22[1]-epar12[1]);
      tfak2 = fabs(epar22[0]-spar22[0]) + fabs(epar22[1]-spar22[1]);
      tdom = tfak1 + tfak2;
      if (DNEQUAL(tdom,DZERO))
        {
	  spar21[0] = (tfak2*epar11[0] + tfak1*epar21[0])/tdom;
	  spar21[1] = (tfak2*epar11[1] + tfak1*epar21[1])/tdom;
	  
	  /* If the two last coefficients are zero, then then this intersection
	     point must be discarded */
	  
	  if (spar21[0]<eval11[0] || eval11[1]<spar21[0] ||
	      spar21[1]<eval12[0] || eval12[1]<spar21[1])
            {
	      /*          Intersection point outside */
	      
	      kbound2 = 0;
            }
        }
      else
        {
	  /* epar1, spar and epar2 has equal last coordinates, since there
	     is an intersection all must lie on the boundary, thus all are
	     inside */
	  
	  kbound2 = 0;
        }
    }
  else if (kstat2==4 && kins1==1)
    {
      /*  On boundary and stepping out */
      goto war04;
    }
  
  /* kbound1 and kbound2 tells if we have got and intersection with the
     boundary */
  
  
  /* If intersections along both boundaries then find which intersection
     is closest to epar1 */
  
  if (kbound1!=0 && kbound2!=0)
    {
/*guen      int t1,t2,t3,t4;*/ /* temporary varuiables */
/*guen changed to           */
      double t1,t2,t3,t4; /* temporary variables */

      t1 = s6dist(spar11,epar11,2);
      t2 = s6dist(spar12,epar12,2);
      t3 = s6dist(spar21,epar11,2);
      t4 = s6dist(spar22,epar12,2);
      
      if ((t1*t1+t2*t2) < (t3*t3+t4*t4) )
        kbound2 = 0;
      else
        kbound1 = 0;
    }
  
  if (kbound1==0 && kbound2 ==0)
    {
      /*  No intersection */
      goto war00;
    }
  else if (kbound1!=0 && kbound2==0)
    {
      /*  Intersection with boundary of first patch */
      memcopy(gpar1,spar11,2,DOUBLE);
      memcopy(gpar2,spar12,2,DOUBLE);
      *jbound = kbound1;
    }
  else if (kbound1==0 && kbound2!=0)
    {
      /*  Intersection with boundary of second patch */
      memcopy(gpar1,spar21,2,DOUBLE);
      memcopy(gpar2,spar22,2,DOUBLE);
      *jbound = kbound2+4;
    }
  
  if (kins1 == 1)
    {
      if (eval11[0] == epar11[0] || epar11[0] == eval11[1] ||
	  eval12[0] == epar11[1] || epar11[1] == eval12[1] ||
	  eval21[0] == epar12[0] || epar12[0] == eval21[1] ||
	  eval22[0] == epar12[1] || epar12[1] == eval22[1])
        {
	  goto war04;
	}
      else
	{
	  goto war02; 
	}
    }
  
  if (kins2 == 1) goto war03; 
  
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
  
  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1330",*jstat,kpos);
  goto out;
  
  
 out:
  return;
}
