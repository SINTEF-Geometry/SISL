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
 * $Id: s9clipimp.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S9CLIPIMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9clipimp(double epar1[],double epar2[],SISLSurf *psurf1,double eimpli[],
	       int ideg,double euval[],double evval[],double aepsge,
	       double gpnt1[],double gpar1[],int *jstat)
#else
void s9clipimp(epar1,epar2,psurf1,eimpli,ideg,euval,evval,
               aepsge,gpnt1,gpar1,jstat)
     double epar1[];
     double epar2[];
     SISLSurf   *psurf1;
     double eimpli[];
     int    ideg;
     double euval[];
     double evval[];
     double aepsge;
     double gpnt1[];
     double gpar1[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To clip the intersection curve between epar1 and epar2
*              against the parameter boundary of the patch defined
*              by euval and evval.
*
*
* INPUT      : epar1  - Parameter pair of start point of curve branch
*              epar2  - Parameter pair of end   point of curve branch
*              psurf1 - Description of B-spline surface              
*              eimpli - Description of implicit surface
*              ideg   - Degree of implicit surface
*                        ideg=1:    Plane              
*                        ideg=2;    Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              euval  - Parameter interval in first parameter direction
*              evval  - Parameter interval in second parameter direction
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in B-spline surface
*              gpar1  - Parameter pair of result of iteration in B-spline
*                       surface
*              jstat  - status messages  
*                       = 2      : Iteration diverged or to many iterations
*                       = 1      : ok, The line cross the boundary, point
*                                  found
*                       = 0      : ok, The line does not cross the boundary
*                       < 0      : error
*
*
* METHOD     :
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 4-July-1988
* Revised by : Tor Dokken, Si, Oslo, Norway, March 1989
*              Taken into consideration that when we step we might step
*              from outside to outside of patch
* Revised by : Mike Floater, SI, 1991-01
*                   Add perspective and circular silhouettes (ideg=1004,ideg=1005)
*
*********************************************************************
*/
{
  int kpos=0;                   /* Position of error                           */
  int klfu=0;                   /* Variable used as pointer in knot vector     */
  int klfv=0;                   /* Variable used as pointer in knot vector     */
  int kder=2;                   /* Number of derivatives to be calculated      */
  int kstat;                    /* Local status variable                       */
  int kdir;                     /* Parameter direction of tpar                 */
  int kcross;                   /* Control variable in while loop              */
  int knbit;                    /* Number of iterations                        */
  int krem;                     /* Remember status                             */
  int kbound;                   /* Numbering of boundary                       */
  int ksize;                    /* Number of doubles for storage of derivateves
				   and normal vector */
  int ksizem3;                  /* ksize - 3                                   */
  double tpar;                  /* Constant parameter value                    */
  double spnt1[33];             /* Coordinates of point                        */
  double spar1[2],spar2[2];     /* Local parameter values                      */
  double spar3[2];              /* Local parameter values                      */
  
  
  
  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles are to be calculated */
  
  if (ideg==1003 || ideg==1004 || ideg==1005)
    {
      kder = 3;
      ksize = 33;
    }
  else
    {
      ksize = 21;
      kder =2;
    }
  ksizem3 = ksize -3;
  
  
  /* Make local copy of parameters */
  
  memcopy(spar1,epar1,2,DOUBLE);                        
  memcopy(spar2,epar2,2,DOUBLE);
  
  kcross = 1;
  knbit  = 0;
  
  while(kcross && knbit<4)
    {
      
      /*  Find intersection between boundary of parameter area and patch */
      
      s1305(spar1,spar2,euval,evval,&kbound,spar3,&kstat);
      if (kstat<0) goto error;
      
      /*  Remember status so that the line can be updated properly */
      krem = kstat;
      
      if (kstat<2 || kbound == 0)
        {
	  /*      No intersection */
	  kcross = 0;
        }
      else
        {
	  
	  /*      Calculate start point */
	  
	  s1421(psurf1,kder,spar3,&klfu,&klfv,spnt1,spnt1+ksizem3,&kstat); 
	  if (kstat<0) goto error;
	  
	  if (kbound==1)
            {
	      kdir = 1;
	      tpar = euval[0];
            }
	  else if (kbound==2)
            {
	      kdir = 2;
	      tpar = evval[1];
            }
	  else if (kbound==3)
            {
	      kdir = 1;                                   
	      tpar = euval[1];
            }
	  else if (kbound==4)
            {
	      kdir = 2;
	      tpar = evval[0];
            }
	  
	  
	  /* Iterate to boundary intersection */
	  
	  s9boundimp(spnt1,gpar1,psurf1,eimpli,ideg,tpar,kdir,aepsge,
		     gpnt1,gpar1,&kstat);
	  if (kstat<0) goto error;
	  if (kstat==2) goto war02;
	  
	  /* Iteration converged, copy output if new loop necessary */
	  
	  
	  if (krem == 2)
            {
	      /* spar1 was inside, update spar2 */
	      memcopy(spar2,gpar1,2,double);
            }
	  else
            {
	      /* spar2 was inside, update spar1 */
	      memcopy(spar1,gpar1,2,double);
            }
	  
	  /* Update number of iterations */
	  
	  knbit++;
        }
    }
  
  
  
  /* Problem solved if kcross==0. In this case we might have two cases:
     - iteration not used: Then knbit=0
     - iteration used    : Then knbit>0
     
     if kcross==1, then we have stopped on the condition knbit>3, and we
     have no sucess. */
  
  if (kcross==0 && knbit ==0)
    {
      /*  Iteration not used because boundary not crossed */
      *jstat = 0;
    }
  else if (kcross==0 && knbit > 0)
    {
      /*  Boundary crossed, point found, more than one intersection point
	  is possible, check which to return */
      
      if (spar1[0] == epar1[0] && spar1[1] == epar1[1])
        {
	  memcopy(gpar1,spar2,2,double);
        }
      else
        {
	  memcopy(gpar1,spar1,2,double);
	  
	  /*  Calculate crossing point, only necessary when we step into
	      the patch since we already could have stepped out and this
	      point is recorded in gpnt1 */
	  
	  s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+ksizem3,&kstat); 
	  if (kstat<0) goto error;
        }
      
      *jstat = 1;
    }
  else
    {
      /*  To many tries */
      goto war02;
    }
  goto out;
  
  /* Iteration without sucess */
  
 war02: *jstat = 2;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s9clipimp",*jstat,kpos);
  goto out;
  
 out:
  return;                                                                        
  
}

