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
 * $Id: refine_all.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define REFINE_ALL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
refine_all (SISLIntdat ** pintdat,
	    SISLObject * po1,
	    SISLObject * po2,
	    double eimpli[],
	    int ideg,
	    double aepsge,
	    int *jstat)

#else
void
refine_all (pintdat,
	    po1,
	    po2,
	    eimpli,
	    ideg,
	    aepsge,
	    jstat)

     SISLIntdat **pintdat;
     SISLObject *po1;
     SISLObject *po2;
     double eimpli[];
     int ideg;
     double aepsge;
     int *jstat;

#endif
/*
*********************************************************************
*
* PURPOSE    : An empty function, (to ensure similarity with other
*              versions).
*
*
* INPUT      : pintdat     - Pointer to pointer to the SISLIntdat data.
*              po1         - Pointer surface object.
*              po2         - Pointer surface object.
*              eimpli      - Array containing descr. of implicit surf
*	       ideg        - Type of impl surf.
              ang_tol     - Angle control tolerance ie ??
*              aepsge      - Absolute tolerance
*
*
* OUTPUT     :  jstat  - status messages
*                       = ?      : ?
*                       = 0      : ok
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway, July-1990
*
*********************************************************************
*/
{
  *jstat = 0;
}
