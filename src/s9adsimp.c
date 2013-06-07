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
 * $Id: s9adsimp.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S9ADSIMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s9adsimp(double epnt1[],double epar1[],double eimpli[],int ideg,double egd1[],
		double epgd1[],double etang[],double eptan[],double astep,int *jstat)
#else
double s9adsimp(epnt1,epar1,eimpli,ideg,egd1,epgd1,etang,eptan,astep,jstat)
     double epnt1[];
     double epar1[];
     double eimpli[];
     int    ideg;
     double egd1[];
     double epgd1[];
     double etang[];
     double eptan[];
     double astep;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To decide if we should step directly on to a guide point
*
* INPUT      : epnt1  - start point, in first surface of iteration step
*              epar1  - Parameter value of epnt1
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1: Plane
*                        ideg=2; Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              egd1   - Guide point in first surface
*                       For ideg=1,2 and 1001 the sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) derivative
*                       and normal. (21 numbers)
*                       For ideg=1003,1004,1005 the second derivatives are followed
*                       by the third derivatives and the normal (33 numbers)
*                       Compatible with output of s1421
*              epgd1  - Parameter value of egd1
*              etang  - Tangent in step direction at epoint
*              eptan  - Tangent in parameter plane at current point
*              idim   - Dimension of space the vectors lie in
*              astep  - Current step length
*
*
* OUTPUT     : s9adsimp - Distance to guide point
*              jstat  - Status  variable
*                        = 0  :  Don't step onto the guide point
*                        = 1  :  Step through the guide point
*
* METHOD     : The step is changed if
*               - The guide point lies in the tangent direction
*               - The distance between the guide point and start point
*                 of iteration step is less than or equal to astep
*
* USE        : This routine is only working in 3-D
*
*
* REFERENCES :
*
*-
* CALLS      : s6diff,s6scpr,s6length,s6crss
*
* WRITTEN BY : Tor Dokken, SI, 1988-june-11
* Revised by : Tor Dokken, SI, 1989-03
*              Corrected stepping through singular points.
* Revised by : Mike Floater, SI, 1991-01
*                   Add perspective and circular silhouettes (ideg=1004,ideg=1005)
*
*********************************************************************
*/
{
  int kpos=1;              /* Position indicator for errors           */
  int kstat;               /* Dummy status variable                   */
  int kdim=3;              /* This routine is only working in 3-D     */
  int k2dim=2;             /* Dimension of parameter plane            */
  int ksize;               /* Number of doubles for storage of derivateves
			      and normal vector */
  int ksizem3;             /* ksize - 3                               */
  double tdum;             /* Variable for storage of reals           */
  double tdist=(double)0.0;/* Distance between guide point and point  */
  double sdiff[3];         /* Vector for difference between two points*/
  double scr1[3],scr2[3];  /* Normal vectors                          */
  double snorm[3];         /* Normal vector                           */
  
  
  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles are to be calculated */
  
  if (ideg==1003 || ideg==1004 || ideg==1005)
    {
      ksize = 33;
    }
  else
    {
      ksize = 21;
    }
  ksizem3 = ksize -3;
  
  /* First see that we are not turning direction in the parameter plane */
  
  s6diff(epgd1,epar1,k2dim,sdiff);
  if (s6scpr(sdiff,eptan,k2dim) < DZERO) goto dontstepthrough;
  
  s6diff(egd1,epnt1,kdim,sdiff);
  tdum  = s6scpr(sdiff,etang,kdim);
  tdist = s6length(sdiff,kdim,&kstat);
  
  /* Step onto point if it is within 2.0*aepsge */
  
  if (tdum > DZERO)
    {
      if (DZERO < tdist && tdist <= (double)2.0*astep)
        {
	  /* Guide point lies within step length and in step direction, test
	     if cross products of normal vectors at current point and guide point
	     point in the same direction, this is not possible to test on for
	   silhouette curves */
	  
	  if (ideg < 1003)
            {
	      
	      /* Make normal to implicit surface at epnt1 */
	      s1308(epnt1,3,eimpli,ideg,snorm,&kstat);
	      if (kstat<0) goto error;
	      
	      /* Make cross product of normals in start point */
	      s6crss(epnt1+ksizem3,snorm,scr1);
	      
	      
	      /* Make normal to implicit surface at epnt1 */
	      s1308(epnt1,3,eimpli,ideg,snorm,&kstat);
	      if (kstat<0) goto error;
	      
	      /* Make cross product of normals in guide point */
	      s6crss(egd1+ksizem3,snorm,scr2);
	      
	      
	      /* Make scalar product of these two vectors     */
	      tdum = s6scpr(scr1,scr2,kdim);
	      if (kstat<0) goto error;
	      
	      /* If positive scalar product the curve at the two points point in
		 the same direction, step through point */
	      if (tdum > DZERO) goto stepthrough;
	      else if (tdum == DZERO)
                {
		  
		  double tl1,tl2;
		  
		  /* Vectors orthogonal or at least one has length zero */
		  
		  tl1 = s6length(scr1,kdim,&kstat);
		  tl2 = s6length(scr2,kdim,&kstat);
		  
		  if (tl1 != DZERO && tl2 != DZERO)
		    goto dontstepthrough;
		  else if (tl1 == DZERO && tl2 == DZERO)
		    goto stepthrough;
		  else if (tl2 == DZERO)
		    goto stepthrough;
		  else if (tl2 != DZERO)
		    {
		      /* Test if scr2 points in the direction from start */
		      
		      tl1 = s6scpr(sdiff,scr2,kdim);
		      
		      if (tl1 < DZERO)
			goto dontstepthrough;
		      else 
			goto stepthrough;
		    }
                }
            }
	  else
            goto stepthrough;
        }
    }
  
 dontstepthrough:
  
  *jstat = 0;
  goto out;
  
 stepthrough:
  *jstat = 1;
  goto out;
  
  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s9adsimp",*jstat,kpos);
  goto out;
  
 out:
  return(tdist);
}
