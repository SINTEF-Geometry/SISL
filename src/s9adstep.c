/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s9adstep.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define S9ADSTEP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s9adstep(double epnt1[],double epar1[],double epnt2[],double epar2[],
		double egd1[],double epgd1[],double egd2[],double epgd2[],
		double etang[],double eptan1[],double eptan2[],double astep,int *jstat)
#else
double s9adstep(epnt1,epar1,epnt2,epar2,egd1,epgd1,egd2,epgd2,etang,
                eptan1,eptan2,astep,jstat)
     double epnt1[];
     double epar1[];
     double epnt2[];
     double epar2[];
     double egd1[];
     double epgd1[];
     double egd2[];
     double epgd2[];
     double etang[];
     double eptan1[];
     double eptan2[];
     double astep;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To decide if we should step directly on to a guide point
*
* INPUT      : epnt1  - start point, in first surface of iteration step
*              epar1  - parameter pair of epnt1
*              epnt2  - start point, in second surface of iteration step
*              epar2  - parameter pair of epnt2
*              egd1   - Guide point in first surface
*              epgd1  - parameter pair of egd1
*              egd2   - Guide point in second surface.
*                       For all these points we have stored the following
*                       values: Position, (1,0)-der, (0,1)-der, normal,
*                       (2,0)-der, (1,1)-der and (0,2)-der. Thus the
*                       normal is stored in (edg2+18).
*              epgd2  - parameter pair of egd2
*              etang  - Tangent in step direction at epnt1 and epnt2
*              eptan1 - Tangent direction in parameter pair at epnt1
*              eptan2 - Tangent direction in parameter pair at epnt2
*              idim   - Dimension of space the vectors lie in
*              astep  - Current step length
*
*
* OUTPUT     : s9adstep - Distance to guide point
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
* Revised by : Tor Dokken, SI, 1989-april-04
*              Correction of stepping onto singular points
*
*********************************************************************
*/
{
  int kstat;               /* Dummy status variable                   */
  int kdim=3;              /* This routine is only working in 3-D     */
  int k2dim=2;             /* Dimension of parameter plane            */
  double tdum;             /* Variable for storage of reals           */
  double tdist=DZERO;      /* Distance between guide point and point  */
  double sdiff[3];         /* Vector for difference between two points*/
  double scr1[3],scr2[3];  /* Normal vectors                          */
  
  
  
  /* First see that we are not turning direction in the parameter plane 1
   */
  
  s6diff(epgd1,epar1,k2dim,sdiff);
  if (s6scpr(sdiff,eptan1,k2dim) <= DZERO) goto dontstepthrough;
  
  /* Then see that we are not turning direction in the parameter plane 2
   */
  
  s6diff(epgd2,epar2,k2dim,sdiff);
  if (s6scpr(sdiff,eptan2,k2dim) <= DZERO) goto dontstepthrough;
  
  
  s6diff(egd1,epnt1,kdim,sdiff);
  tdum  = s6scpr(sdiff,etang,kdim);
  tdist = s6length(sdiff,kdim,&kstat);
  *jstat = 0;
  if (tdum > DZERO)
    {
      
      /* Step onto point if it is within 2.0*astep */
      
      if (DZERO < tdist && tdist <= (double)2.0*astep)
        {
	  /* Guide point lies within step length and in step direction, test
	     if cross products of normal vectors at current point and guide point
	     point in the same direction */
	  
	  /* Make cross product of normals in start point */
	  s6crss(epnt1+18,epnt2+18,scr1);
	  
	  /* Make cross product of normals in guide point */
	  s6crss(egd1+18,egd2+18,scr2);
	  
	  /* Make scalar product of these two vectors     */
	  tdum = s6scpr(scr1,scr2,kdim);
	  
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
    }
  
 dontstepthrough:
  
  *jstat = 0;
  goto out;
  
 stepthrough:
  *jstat = 1;
  goto out;
  
 out:
  return(tdist);
}


