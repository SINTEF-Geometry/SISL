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
 * $Id: sh6ptouchp.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6PUTTOUCH
#include "sislP.h"                            

#if defined (SISLNEEDPROTOTYPES)
void
      sh6puttouch( SISLIntpt *psource, SISLIntpt *pdest, int seq)
#else
	 
	 void sh6puttouch(psource,pdest, seq)
	    
	    SISLIntpt *psource,*pdest;
	    int seq;
#endif

/*
*********************************************************************
*                                                                   
* PURPOSE    : To set right direction in a touch point (Singular, but no branch).
*              The invariant to hold is that the direction on the curve 
*              in the parameter space is parallel to the delta vector between
*              the to points.
*
*
*
* INPUT      : psource - The previous point on the curves.
*              pdest   - The point that is to be fit in.
*                        Posisiton must be ok.
*              seq     - The sequncing of the points
*                        +1 - psource comes before pdest.
*                        -1 - psource comes after pdest.
* OUTPUT     : pdest   - Geometric values, ie tangent is changed.

*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : s6scpr, s6diff, s6norm
*
* WRITTEN BY : UJK
*
*********************************************************************
*/
{
  int kdim2=2;
  int kdim3=3;
  double diffv[3];
  int ki;
  double dot;
  /* ______________________ */
  
  if (psource->iinter == SI_ORD)
    sh6putsing(psource, pdest);
  else
    {
       kdim2 = 2;
       s6diff(pdest->geo_track_2d_1,psource->geo_track_2d_1,kdim2,diffv);
       dot = s6scpr(pdest->geo_track_2d_1+kdim2, diffv, kdim2);
       if (dot * seq < 0)
	 {
	    /* Turn direction */
	    for (ki=0;ki<kdim2;ki++) 
	      {
		 pdest->geo_track_2d_1[kdim2+ki] = -pdest->geo_track_2d_1[kdim2+ki];
		 pdest->geo_track_2d_2[kdim2+ki] = -pdest->geo_track_2d_2[kdim2+ki];
	      }  
	    
	    for (ki=0;ki<kdim3;ki++) 
	      pdest->geo_track_3d[kdim3+ki] = -pdest->geo_track_3d[kdim3+ki];
	 }
       
    }
}
