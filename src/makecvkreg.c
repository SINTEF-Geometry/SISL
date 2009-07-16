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
 * $Id: makecvkreg.c,v 1.7 1994-11-30 14:37:16 pfu Exp $
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
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-08. Added error propagation.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov.1994. Set cuopen flag to
*              closed when changed from periodic.
**********************************************************************/
{
   int kn=pc->in;	/* Number of vertices in 1. par. dir.  */
   int kk=pc->ik;	/* Order in 1. par. dir.               */
   /* --------------------------------------------------------- */
   /* Pick part of curve */
   s1712 (pc, pc->et[kk-1], pc->et[kn], rcnew, jstat);
  if (*jstat < 0)  goto error;

   if (pc->cuopen == SISL_CRV_PERIODIC )
     (*rcnew)->cuopen = SISL_CRV_CLOSED;

  goto out;

  /* Error in lower level routine */
error:
  s6err ("make_cv_kreg", *jstat, 0);

out:;

}
