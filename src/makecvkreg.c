/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: makecvkreg.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define MAKE_CV_KREG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
    make_cv_kreg (SISLCurve * pc, SISLCurve ** rcnew, int *jstat)
#else
void 
   make_cv_kreg (pc, rcnew, jstat)
     SISLCurve *pc;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a curve to a k-regular basis.
*
*
*
* INPUT      : pc	- Curve to be made k-regular.
*
*
*
* OUTPUT     : rcnew	- The new curve on a k-regular basis.
*              jstat	- status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
*
* WRITTEN BY : Ulf J. Krystad, SI, 04.92.
*
**********************************************************************/
{
   int kn=pc->in;	/* Number of vertices in 1. par. dir.  */
   int kk=pc->ik;	/* Order in 1. par. dir.               */
   /* --------------------------------------------------------- */
   /* Pick part of curve */
   s1712 (pc, pc->et[kk-1], pc->et[kn], rcnew, jstat);
   
}
