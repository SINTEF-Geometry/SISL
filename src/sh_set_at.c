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
 * $Id: sh_set_at.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH_SET_AT

#include "sislP.h"




#if defined(SISLNEEDPROTOTYPES)
void
sh_set_at (SISLObject * po1, SISLObject * po2, 
	SISLIntdat * pintdat, int *jstat)
#else
void
sh_set_at (po1, po2, pintdat,  jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set SI_AT topology part.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintdat  - Intersection data of object-object intersection.
*
* OUTPUT       : jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh1781_AT          - Find pre-topology of 1D curve-point.
*              sh1780_AT          - Find pre-topology of curve-curve.
*              sh1779_AT          - Find pre-topology of 3D curve-surface.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 09.91
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int ki;			/* Counter.                                */
  int kdim;			/* Dimension of geometry space.            */
  SISLIntpt *qpt = SISL_NULL;	/* Pointer to intersection point.          */
  /* --------------------------------------------------------------------- */

  /* Init */
  *jstat = 0;

  /* Test if an intersection data structure exist.  */
  if (pintdat == SISL_NULL)
    goto out;


  /* Fetch dimension of geometry space. */
  if (po1->iobj == SISLPOINT)

    kdim = po1->p1->idim;
  else if (po1->iobj == SISLCURVE)
    kdim = po1->c1->idim;
  else
    kdim = po1->s1->idim;

  /* Treat only cases:
     crv vs pt 1D
     crv vs crv
     crv vs sf
     (?sf vs pt 2D)
     */

  if (!(((po1->iobj == SISLCURVE && po2->iobj >= SISLCURVE) ||
	 (po2->iobj == SISLCURVE && po1->iobj >= SISLCURVE)) ||
	(kdim == 1 && (po1->iobj + po2->iobj) == (SISLPOINT + SISLCURVE)) ||
	(kdim == 2 && (po1->iobj + po2->iobj) == (SISLPOINT + SISLSURFACE))))
    goto out;


  for (ki = 0; ki < (pintdat)->ipoint; ki++)
  {
     qpt = (pintdat)->vpoint[ki];

      /* Browse on the dimension of geometry space and the type of
         the input objects.     */

      if (kdim == 1 && ((po1->iobj == SISLCURVE && po2->iobj == SISLPOINT)
		     || (po2->iobj == SISLCURVE && po1->iobj == SISLPOINT)))
	{
	  /* Compute pre-topology in one-dimensional curve-level value
             intersection.            */

	  sh1781_at (po1, po2,qpt, &kstat);
	  if (kstat < 0)
	    goto error;
	}
      else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE)
	{
	  /* curve-curve intersection.  */
	  sh1780_at (po1, po2, qpt, &kstat);
	  if (kstat < 0)
	    goto error;
	}
      else if (kdim == 3 &&
	       ((po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE) ||
		(po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE)))
	{
	  /* Surface-curve intersection in 3-dimensional geometry space. */

	  sh1779_at (po1, po2, qpt, &kstat);
	  if (kstat < 0)
	    goto error;
	}
    }

  /* Task performed.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
   return;
}
