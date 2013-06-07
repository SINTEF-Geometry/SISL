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
 * $Id: s9clipit.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S9CLIPIT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9clipit(double epar11[],double epar12[],double epar21[],double epar22[],
	      SISLSurf *psurf1,SISLSurf *psurf2,double euval[],double evval[],double esval[],
	      double etval[],double aepsge,double gpnt1[],double gpnt2[],
	      double gpar1[],double gpar2[],int *jstat)
#else
void s9clipit(epar11,epar12,epar21,epar22,psurf1,psurf2,
              euval,evval,esval,etval,
              aepsge,gpnt1,gpnt2,gpar1,gpar2,jstat)
     double epar11[];
     double epar12[];
     double epar21[];
     double epar22[];
     SISLSurf   *psurf1;
     SISLSurf   *psurf2;
     double euval[];
     double evval[];
     double esval[];
     double etval[];
     double aepsge;
     double gpnt1[];
     double gpnt2[];
     double gpar1[];
     double gpar2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To clip the intersection curve between epar1 and epar2
*              against the parameter boundary of the patch 1 defined
*              by euval and evval, and the parameter boundary of the
*              patch 2 defined by esval and etval.
*
*
* INPUT      : epar11  - Parameter pair of start point in first surface
*              epar12  - Parameter pair of start point in second surface
*              epar21  - Parameter pair of end point in first surface
*              epar22  - Parameter pair of end point in second surface
*              psurf1 - Description of B-spline surface 1
*              psurf2 - Description of B-spline surface 2            
*              euval  - Parameter interval in first parameter direction
*              evval  - Parameter interval in second parameter direction
*              esval  - Parameter interval in third parameter direction
*              etval  - Parameter interval in fourth parameter direction
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in B-spline surface 1
*              gpnt2  - 0-2 Derivatives + normal of result of iteration
*                       in B-spline surface 2
*              gpar1  - Parameter pair of result of iteration in B-spline
*                       surface 1
*              gpar2  - Parameter pair of result of iteration in B-spline
*                       surface 2
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
* Revised by : Tor Dokken, SI, Oslo, Norway, 4-APril-1989
*              Correction of steps over two boundaries
*
*********************************************************************
*/
{
  int kpos=0;                   /* Position of error                       */
  int klfu=0;                   /* Variable used as pointer in knot vector */
  int klfv=0;                   /* Variable used as pointer in knot vector */
  int klfs=0;                   /* Variable used as pointer in knot vector */
  int klft=0;                   /* Variable used as pointer in knot vector */
  int kder=2;                   /* Number of derivatives to be calculated  */
  int kstat;                    /* Local status variable                   */
  int kdir;                     /* Parameter direction of tpar             */
  int kcross;                   /* Control variable in while loop          */
  int knbit;                    /* Number of iterations                    */
  int krem;                     /* Remember status                         */
  int kbound;                   /* Numbering of boundary                   */
  double tpar;                  /* Constant parameter value                */
  double spnt1[21];             /* Coordinates of point                    */
  double spnt2[21];             /* Coordinates of point                    */
  double spar11[2],spar12[2];   /* Local parameter values                  */
  double spar21[2],spar22[2];   /* Local parameter values                  */
  double spar31[2],spar32[2];   /* Local parameter values                  */
  
  /* Make local copy of parameters */
  
  memcopy(spar11,epar11,2,DOUBLE);
  memcopy(spar12,epar12,2,DOUBLE);
  memcopy(spar21,epar21,2,DOUBLE);
  memcopy(spar22,epar22,2,DOUBLE);
  
  kcross = 1;
  knbit  = 0;
  
  while(kcross && knbit<8)
    {
      
      /* Find intersection between boundary of parameter area and patch */
      
      s1330(spar11,spar12,spar21,spar22,euval,evval,esval,etval,
	    &kbound,spar31,spar32,&kstat);
      if (kstat<0) goto error;
      
      /* Remember status so that the line can be updated properly */

      krem = kstat;
      
      if (kstat<2 || kbound == 0)
        {
	  /* No intersection */

	  kcross = 0;
        }
      else
        {
	  
	  /* Calculate start point for iteration */
	  
	  s1421(psurf1,kder,spar31,&klfu,&klfv,spnt1,spnt1+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  s1421(psurf2,kder,spar32,&klfs,&klft,spnt2,spnt2+18,&kstat); 
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
	  else if (kbound==5)
            {
	      kdir = 3;
	      tpar = esval[0];
            }
	  else if (kbound==6)
            {
	      kdir = 4;
	      tpar = etval[1];
            }
	  else if (kbound==7)
            {
	      kdir = 3;                                   
	      tpar = esval[1];
            }
	  else if (kbound==8)
            {
	      kdir = 4;
	      tpar = etval[0];
            }
	  
	  
	  /* Iterate to boundary intersection */
	  
	  s9boundit(spnt1,spnt2,spar31,spar32,psurf1,psurf2,tpar,kdir,aepsge,
		    gpnt1,gpnt2,gpar1,gpar2,&kstat);
	  if (kstat<0) goto error;
	  if (kstat==2) goto war02;
	  
	  /* Iteration converged, copy output if new loop necessary */
	  
	  if (krem == 2)
            {
	      /* spar1 was inside, update spar2 */

	      memcopy(spar21,gpar1,2,double);
	      memcopy(spar22,gpar2,2,double);
            }
	  else
            {
	      /* spar2 was inside, update spar1 */

	      memcopy(spar11,gpar1,2,double);
	      memcopy(spar12,gpar2,2,double);
            }
	  
	  /* Update number of iterations */
	  
	  knbit++;
        }
    }
  
  
  
  /* Problem solved if kcross==0. In this case we might have two cases:
     - iteration not used: Then knbit=0
     - iteration used    : Then knbit>0
     
     if kcross==1, then we have stopped on the condition knbit>7, and we
     have no success. */
  
  if (kcross==0 && knbit ==0)
    {
      /* Iteration not used because boundary not crossed */

      *jstat = 0;
    }
  else if (kcross==0 && knbit > 0)
    {
      /* Boundary crossed, point found, more than one intersection point
	 is possible, check which to return */
      
      if (spar11[0] == epar11[0] && spar11[1] == epar11[1] &&
	  spar12[0] == epar12[0] && spar12[1] == epar12[1] )
        {
	  memcopy(gpar1,spar21,2,double);
	  memcopy(gpar2,spar22,2,double);
        }
      else
        {
	  memcopy(gpar1,spar11,2,double);
	  memcopy(gpar2,spar12,2,double);
	  
	  /* Calculate crossing point, only necessary when we step into
	     the patch since we already could have stepped out and this
	     point is recorded in gpnt1 */
	  
	  s1421(psurf1,kder,gpar1,&klfu,&klfv,gpnt1,gpnt1+18,&kstat); 
	  if (kstat<0) goto error;
	  
	  s1421(psurf2,kder,gpar2,&klfu,&klfv,gpnt2,gpnt2+18,&kstat); 
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
  
  /* Iteration without success */
  
  war02: 
    *jstat = 2;
    goto out;
  
  /* Error in lower level routine.  */
  
  error : 
    *jstat = kstat;
    s6err("s9clipit",*jstat,kpos);
    goto out;
  
  out:
    return;
}
