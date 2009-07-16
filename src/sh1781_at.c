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
 * $Id: sh1781_at.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1781_AT

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
sh1781_at (SISLObject * po1, SISLObject * po2,
	   SISLIntpt * pintpt,
	   int *jstat)
#else
void
sh1781_at (po1, po2, pintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pintpt;
     int *jstat;
#endif
/*
******************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data AT in 1-dimensional curve-point
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point after updating pre-topology.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh6gettop - Get topological info
*              sh6settop - Set topological info
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.      */
  double *st;			/* Pointer to knot vector of curve.        */
  double *sptpar = pintpt->epar;/* Pointer to parameter array of int.pt.   */
  double tref;			/* Referance value in equality test.       */
  SISLCurve *qc;		/* Pointer to current curve.               */
  /* --------------------------------------------------------------------- */


  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Set pointers into the arrays storing pre-topology information. */

  if (po1->iobj == SISLCURVE)
    {
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  /* Get pre-topology information. */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);
  if (kstat < 0)
    goto error;

  /* Test dimension of geometry space. */
  if (po1->iobj == SISLCURVE)
    qc = po1->c1;
  else
    qc = po2->c1;

  /* Store curve information in local parameters. */
  kn = qc->in;
  kk = qc->ik;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  /* Test if the intersection point lies at an endpoint of
     the curve. */

  if (DEQUAL (sptpar[0] + tref, st[kn] + tref))
    *lr1 = SI_AT;
  if (DEQUAL (sptpar[0] + tref, st[kk - 1] + tref))
    *ll1 = SI_AT;

  /* Update pretopology of intersection point.  */

  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);
  if (kstat < 0)
    goto error;

  *jstat = 0;
  goto out;


  /* Error lower level routine.  */
error:*jstat = kstat;
  goto out;

out:
  return;
}
