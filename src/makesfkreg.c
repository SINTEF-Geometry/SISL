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
 * $Id: makesfkreg.c,v 1.6 1994-11-30 12:53:02 pfu Exp $
 *
 */


#define MAKE_SF_KREG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void make_sf_kreg (SISLSurf * ps, SISLSurf ** rsnew, int *jstat)
#else
void
   make_sf_kreg (ps, rsnew, jstat)
     SISLSurf *ps;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a surface to a k-regular basis.
*
*
*
* INPUT      : ps	- Surface to be made k-regular.
*
*
*
* OUTPUT     : rsnew	- The new surface on a k-regular basis.
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
*
**********************************************************************/
{
  int kn1=ps->in1;	/* Number of vertices in 1. par. dir.  */
  int kn2=ps->in2;	/* Number of vertices in 2. par. dir.  */
  int kk1=ps->ik1;	/* Order in 1. par. dir.               */
  int kk2=ps->ik2;	/* Order in 2. par. dir.               */
  /* --------------------------------------------------------- */

  s1001 (ps, ps->et1[kk1-1], ps->et2[kk2-1],
		ps->et1[kn1], ps->et2[kn2], rsnew, jstat);
  if (*jstat < 0)  goto error;

  goto out;

  /* Error in lower level routine */
error:
  s6err ("make_sf_kreg", *jstat, 0);

out:;

}
